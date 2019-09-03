
# ## Case study 3: Model-parameter estimation in the ESDL
# ### Example of the temperature sensitivity of ecosystem respiration
#
# #### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)
#
# * Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).
#
# * The NB is written based on Julia 1.1
#
# * Normal text are explanations referring to notation and equations in the paper
#
# * `# comments in the code are itended explain specific aspects of the coding`
#
# * ### New steps in workflows are introduced with bold headers
#
# Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

## for plotting later on (need to be loaded first, to avoid conflicts)
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using PyCall, PyPlot, PlotUtils

## for operating the Earth system data lab
using ESDL, ESDLPlots

## for parallel computing
using Distributed

## other relevant packages
using Dates, Statistics, DataFrames, MultivariateStats
#----------------------------------------------------------------------------

# ### Select and prepare (subset/gapfill) an Earth system data cube
#
# We need to choose a cube and here select a 8-dayily, 0.25° resolution global cube. The cube name suggests it is chunked such that we have one time chunk and 720x1440 spatial chunks

cube_handle = Cube("/scratch/DataCube/v2.0.0/esdc-8d-0.25deg-1x720x1440-2.0.0.zarr/")

## one cube for the tair respiration
world_rec = getCubeData(c, variable = ["air_temperature_2m","terrestrial_ecosystem_respiration"],
    region = "Europe",
    time = (Date("2000-01-01"), Date("2016-02-26")))
world_gpp = getCubeData(c, variable = ["gross_primary_productivity"],
    region = "Europe",
    time = (Date("2000-01-01"), Date("2016-02-26")))
world_h2o = getCubeData(c, variable = ["surface_moisture"],
    region = "Europe",
    time = (Date("2000-01-01"), Date("2016-02-26")))

world_rec = gapFillMSC(world_rec)
world_gpp = gapFillMSC(world_gpp)
world_h2o = gapFillMSC(world_h2o)
#----------------------------------------------------------------------------

# ### Function we need: moving average for decomposition

## Moving Average decomposes a singal into fast and slow oscillations
## by caluclating a moving average over a window of points.
## This creates a smoothed curve (slow osc.) which can be subtracted from the original singla,
## to obtain fast oscillations separately.
function movingAverage(xin, window)
    Z     = length(xin)
    ## calculate moving average over window
    ## truncating windows for data points at beginning and end
    movAv = map(1:Z) do i
        if i <= window
            mean(xin[1:i+window])
        elseif i >= (Z-window)
            mean(xin[i-window:end])
        else
            mean(xin[i-window:i+window])
        end
    end
    ## return slow oscillations in col 1 and fast oscillations in col 2
    xout = hcat(movAv, xin .- movAv)
    return xout
end
#----------------------------------------------------------------------------

# ### Function we need: SCAPE Q10

function Q10decomposed(xout_p, xout_t, xin)
    ## xin: array with columns "air_temperature_2m","terrestrial_ecosystem_respiration"

    ## simple Q10
    tair = xin[:,1]
    resp = xin[:,2]
    tau  = (tair .- (273.15+15))./10
    rho  = log.(resp)

    ## solve the regression
    b    = cor(tau, rho)*std(rho)/std(tau)
    a    = mean(rho) - b*mean(tau)

    Q10  = exp(b)
    Rb   = exp(a)

    corQ10 = cor(tau,rho)^2

    ## signal decomposition; extraction of fast band
    tau_deco = movingAverage(tau, 4)
    rho_deco = movingAverage(rho, 4)

    tau_slow = tau_deco[:, 1]
    tau_fast = tau_deco[:, 2]
    rho_slow = rho_deco[:, 1]
    rho_fast = rho_deco[:, 2]

    ## EQ S5 - reference to the orginal SCAPE Q10 paper Mahecha et al. (2010)
    ## Q10 calculated on fast oscillations only
    d    = cor(tau_fast, rho_fast)*std(rho_fast)/std(tau_fast)
    c    = mean(rho_fast) - b*mean(tau_fast)
    Q10_fast = exp(d)

    ## EQ S6: Influence of low frequency temperature on Rb
    rho_sc = (tau_slow .+ mean(tau)) .* d

    ## EQ S7: Time varying estimate for Rb
    rho_b  = rho_slow .+ mean(rho) .- rho_sc
    Rb_sc  = exp.(rho_b)

    Q10_delta = Q10_fast/Q10
    Q10_diff  = Q10_fast-Q10
    corQ10_fast = cor(tau_fast,rho_fast)^2

    xout_p[1] = Rb
    xout_p[2] = Q10
    xout_p[3] = Q10_fast
    xout_p[4] = Q10_delta
    xout_p[5] = Q10_diff
    xout_p[6] = corQ10
    xout_p[7] = corQ10_fast

    xout_t[:] = Rb_sc
end
#----------------------------------------------------------------------------

# ### Phase lag function to understand the relation with low freq Rb

function phaselag(xout, xin)
    NpY = 46
    Ntot = 365
    ##any(ismissing,xin) && return xout = [missing,missing,missing,missing]
    try
        data = xin
        data = map(i->ismissing(i) ? NaN : i,data)
        x  = data[:, 3] .+ 1.5 .* abs.(minimum(data[:, 3]))
        y  = data[:, 1] .+ 1.5 .* abs.(minimum(data[:, 1]))
        dx = [diff(x);x[1]-x[end]]
        X  = hcat(x,dx)
        (b, c, interc) = llsq(X, y)
        yp = x .* b .+ c .* dx .+ interc
        R2 = var(yp)/var(y)
        xout .= [b,c,atan(-c * 2 * pi / NpY/b) * Ntot / (2*pi),R2]
    catch
        xout .= [missing,missing,missing,missing]
    end
end
#----------------------------------------------------------------------------

# ### We need two output cubes of different length

outAxis1 = OutDims(CategoricalAxis("Parameter", ["Rb", "Q10_allOsc", "Q10_fastOsc", "deltaQ10", "diffQ10", "R2_tau_rho", "R2_tau_fast_rho_fast"]))
outAxis2 = OutDims("Time")
#----------------------------------------------------------------------------

q10,rb  = mapCube(Q10decomposed,
             world_rec,
             indims  = InDims("Time", "Variable", filter = ESDL.DAT.AnyMissing()),
             outdims = (outAxis1,outAxis2))
#----------------------------------------------------------------------------

saveCube(q10, "/Net/Groups/BGI/scratch/mmahecha/test/q10europe/")
saveCube(rb, "/Net/Groups/BGI/scratch/mmahecha/test/rbeurope/")
#----------------------------------------------------------------------------

q10 = loadCube("/Net/Groups/BGI/scratch/mmahecha/test/q10/")
#----------------------------------------------------------------------------

cd("/Net/Groups/BGI/scratch/mmahecha/tempcubes/")
exportcube(q10, "q10b.nc")
#----------------------------------------------------------------------------

plotMAP(q10, var = "diffQ10")
#----------------------------------------------------------------------------

q10.axes[1].values
#----------------------------------------------------------------------------

# ### Now we leave the data cube context and plot

data = q10[3, :, :]'
data = replace(data, missing => NaN)
data = replace(data, 1 => NaN)
data = collect(data)
q10_sc = data

##Q10_norm = copy(data)
Q10_scape = copy(data)
##Q10_diff = copy(data)

lons = collect(q10.axes[2].values)
lats = collect(q10.axes[3].values)

using PyPlot, PyCall

@pyimport mpl_toolkits.basemap as basemap
@pyimport numpy as np
@pyimport matplotlib.colors as mcolors
@pyimport matplotlib.pyplot as plter

## Make room for the ridiculously large title.
##plt.subplots_adjust(top=0.8)

##f   = figure("pyplottest", figsize=(10, 10))
##ax1 = f[:add_subplot](3,1,1)
mymap = basemap.Basemap(
    llcrnrlon = minimum(lons)-0.25/2, llcrnrlat = maximum(lons)+0.25/2,
    urcrnrlon = minimum(lats)-0.25/2, urcrnrlat = maximum(lats)+0.25/2,
    resolution="l")
##mymap[:drawcoastlines](color="0.6", linewidth=0.5)
##mymap[:drawcountries](color="0.6", linewidth=0.5)
##mymap[:drawmapboundary](fill_color="aqua")
##mymap[:fillcontinents](color="0.8", lake_color="white",zorder=0)
xx, yy = np.meshgrid(lons, lats)
xi, yi = mymap(xx, yy)
##mymap[:pcolormesh](xi, yi, Q10_diff, cmap = cmap, vmin=1, vmax=2)

cmap = get_cmap("viridis_r")
##cmap = get_cmap("RdYlBu")
##mdata  = maskoceans(xi, yi, data, resolution = 'h', grid = 0.25, inlands=true)
mymap[:imshow](reverse(Q10_scape, dims=(1)),cmap=cmap, vmin=0, vmax=3)

mymap[:drawmapboundary](fill_color=[0.2, 0.2, 0.2])
cbar = mymap[:colorbar](location = "bottom", extend = "max")
cbar[:ax][:tick_params](labelsize = 10)
cbar[:set_label]("Q10scape")
##ax1[:get_figure]()
##ax1[:set_title]("TREND: Covariance of precipitation vs temperature")
##plter.savefig("/Net/Groups/BGI/scratch/mmahecha/test.png", dpi=1000)
#----------------------------------------------------------------------------

data = q10[3, :, :]'
data = replace(data, missing => NaN)
data = replace(data, 1 => NaN)
data = collect(data)
q10_sc = data

##Q10_norm = copy(data)
Q10_scape = copy(data)
##Q10_diff = copy(data)

lons = collect(q10.axes[2].values)
lats = collect(q10.axes[3].values)

#----------------------------------------------------------------------------

## compute the mean seasonal cycles
msc_rb  = getMSC(rb)
msc_gpp = getMSC(world_gpp)
msc_h2o = getMSC(world_h2o)
#----------------------------------------------------------------------------

new_var_axis = CategoricalAxis("new_var_axis", ["msc_rb", "msc_gpp", "msc_h2o"])

## we generate a new cube that has both results
merged_cube  = concatenateCubes([msc_rb, msc_gpp, msc_h2o], new_var_axis)
#----------------------------------------------------------------------------

Rb_lag = mapCube(phaselag, merged_cube,
    indims  = InDims("MSC", "new_var", filter = ESDL.DAT.AnyMissing()),
    outdims = OutDims(CategoricalAxis("Parameter", ["b", "c", "lag", "R2"])))
#----------------------------------------------------------------------------

Rb_lag_h2o = mapCube(phaselag, merged_cube,
    indims  = InDims("MSC", "new_var", filter = ESDL.DAT.AnyMissing()),
    outdims = OutDims(CategoricalAxis("Parameter", ["b", "c", "lag", "R2"])))
#----------------------------------------------------------------------------

GPP_lag_h2o = mapCube(phaselag, merged_cube,
    indims  = InDims("MSC", "new_var", filter = ESDL.DAT.AnyMissing()),
    outdims = OutDims(CategoricalAxis("Parameter", ["b", "c", "lag", "R2"])))
#----------------------------------------------------------------------------

test = merged_cube[:, 90, 70, :]
plot(test[:, 1], test[:, 2])
#----------------------------------------------------------------------------

data = Rb_lag[3, :, :]'
data = replace(data, missing => NaN)
data = replace(data, 1 => NaN)
data = collect(data)
timelag = copy(data)

lons = collect(q10.axes[2].values)
lats = collect(q10.axes[3].values)


using PyPlot, PyCall

@pyimport mpl_toolkits.basemap as basemap
@pyimport numpy as np
@pyimport matplotlib.colors as mcolors
@pyimport matplotlib.pyplot as plter


## Set up orthographic map projection with perspective of satellite looking down at 45N, 100W.
## Use low resolution coastlines.

##f   = figure("pyplot_subplot_touching", figsize=(10, 10))
##ax1 = f[:add_subplot](3,1,1)
mymap = basemap.Basemap(
    llcrnrlon = minimum(lons)-0.25/2, llcrnrlat = maximum(lons)+0.25/2,
    urcrnrlon = minimum(lats)-0.25/2, urcrnrlat = maximum(lats)+0.25/2,
    resolution="l")
##mymap[:drawcoastlines](color="0.6", linewidth=0.5)
##mymap[:drawcountries](color="0.6", linewidth=0.5)
##mymap[:drawmapboundary](fill_color="aqua")
##mymap[:fillcontinents](color="0.8", lake_color="white",zorder=0)
xx, yy = np.meshgrid(lons, lats)
xi, yi = mymap(xx, yy)
##mymap[:pcolormesh](xi, yi, Q10_diff, cmap = cmap, vmin=1, vmax=2)

##cmap = get_cmap("inferno_r")
cmap = get_cmap("Spectral")
##mdata  = maskoceans(xi, yi, data, resolution = 'h', grid = 0.25, inlands=true)
mymap[:imshow](reverse(timelag, dims=(1)),cmap=cmap, vmin = -10, vmax = 10)

mymap[:drawmapboundary](fill_color=[0.2, 0.2, 0.2])
cbar = mymap[:colorbar](location = "bottom", extend = "both")
cbar[:ax][:tick_params](labelsize = 10)
cbar[:set_label]("Time lag in days")
##ax1[:get_figure]()
##ax1[:set_title]("TREND: Covariance of precipitation vs temperature")
##plter.savefig("/Net/Groups/BGI/scratch/mmahecha/test.png", dpi=1000)
#----------------------------------------------------------------------------

# ### Q10 - Example Europe

plotMAP
#----------------------------------------------------------------------------

using Plots; gr()
#----------------------------------------------------------------------------

o[1, :, :] |> Plots.heatmap
#----------------------------------------------------------------------------

plotMAP(q10_eu, symmetric = TRUE)
#----------------------------------------------------------------------------

# ### Q10 - Example World


using ProgressMeter
@loadOrGenerate q10_world=>"WorldQ10" begin
    q10_world = mapCube(Q10decomposed, world,indims=InDims("Time","Variable",miss=ESDL.NaNMissing()), outdims=OutDims(outAxisQ10,miss=ESDL.NaNMissing()))
end
plotMAP(q10_world, dmax=5)
#----------------------------------------------------------------------------

rm("worldq10param.nc")
exportcube(q10_world,"worldq10param.nc")
#----------------------------------------------------------------------------

## Q10 (all oscillations)
display("Q10 (all oscillations)")
display(plotMAP(q10_world, dmax=4, param=2))

## Q10 (fast oscillations)
display("Q10 (fast oscillations)")
display(plotMAP(q10_world, dmax=3, param=3))

## delta Q10 (Q10fast/Q10)
display("delta Q10 (Q10fast/Q10)")
display(plotMAP(q10_world, dmax=4, param=4))

#----------------------------------------------------------------------------

## R^2 tau/rho
display("R^2 tau/rho")
display(plotMAP(q10_world, dmax=1, param=5))

## R^2 tau_fast/rho_fast
display("R^2 tau_fast/rho_fast")
display(plotMAP(q10_world, dmax=1, param=6))

## Rb
##display("Rb")
##display(plotMAP(q10_world, dmax=5, param=1))
#----------------------------------------------------------------------------

using PyPlot
ioff() ## Interactive plotting OFF, necessary for inline plotting in IJulia

fig = figure("test_map", figsize=(10,10)) ## Not strictly required
ax = axes() ## Not strictly required

from __future__ import division
import numpy as np
from mpl_toolkits.basemap import Basemap

def mapformat():

  m = Basemap(projection='robin', lon_0=0,resolution='c')
  ## resolution c, l, i, h, f in that order

  m.drawmapboundary(fill_color='white', zorder=-1)
  m.fillcontinents(color='0.8', lake_color='white', zorder=0)

  m.drawcoastlines(color='0.6', linewidth=0.5)
  m.drawcountries(color='0.6', linewidth=0.5)

  m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
  m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

  return m
#----------------------------------------------------------------------------
