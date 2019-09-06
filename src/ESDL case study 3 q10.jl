
# ## Case study 3: Model-parameter estimation in the ESDL
# ### Example of the temperature sensitivity of ecosystem respiration
# 
# #### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)
# 
# * Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).
# 
# * The NB is written based on Julia 1.2
# 
# * Normal text are explanations referring to notation and equations in the paper
# 
# * `# comments in the code are intended explain specific aspects of the coding`
# 
# * ### New steps in workflows are introduced with bold headers
# 
# Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

## for plotting later on (need to be loaded first, to avoid conflicts)
using PyCall, PyPlot, PlotUtils

## for operating the Earth system data lab
using ESDL, ESDLPlots

## other relevant packages
using Statistics
#----------------------------------------------------------------------------

# ### Select and subet an Earth system data cube
# 
# We need to choose a cube and here select a 8-dayily, 0.25° resolution global cube. The cube name suggests it is chunked such that we have one time chunk and 720x1440 spatial chunks

cube_handle = Cube("../data/subcube")
#----------------------------------------------------------------------------

# In this case it is better to have one cube for the Tair and one for terrestrial ecosystem respiration R$_{eco}$

world_tair = subsetcube(cube_handle, variable = "air_temperature_2m", time = 2001:2012)
world_resp = subsetcube(cube_handle, variable = "terrestrial_ecosystem_respiration", time = 2001:2012)
#----------------------------------------------------------------------------

# The objective is to estimate  $Q_{10}$ from the decomposed time series. For details we refere the reader to 
# Mahecha, M.D. et al. (2010) *Global convergence in the temperature sensitivity of respiration at ecosystem level.* Science, 329, 838-840.
#         
# The first step is to transformation of both variables, so that the $Q_{10}$ model becomes linear and Tair the exponent:

## Element-wise transformations using `map` are done in a lazy manner, so the
## transformation will be applied only when the data is read or further processed
world_τ = map(tair -> (tair - (273.15+15))/10, world_tair)
world_ρ = map(log, world_resp)

## ... and we combine them into a Data Cube again
world_new = concatenateCubes(τ=world_τ, ρ=world_ρ)
#----------------------------------------------------------------------------

# ### Function we need: moving average for decomposition
# 
# First we need a function for time-series filtering. Using a moving average filter is the simplest way to decomposes a singal into fast and slow oscillations by caluclating a moving average over a window of points. This creates a smoothed curve (slow osc.) which can be subtracted from the original singlal to obtain fast oscillations separately. We could halve likewise used FFTs, SSA, EMD, or any other method for discrete time-series decomposition.

## Moving Average decomposes a singal into fast and slow oscillations
## by caluclating a moving average over a window of points.
## This creates a smoothed curve (slow osc.) which can be subtracted from the original singla,
## to obtain fast oscillations separately.
function movingAverage(xout, xin; windowsize = 4)
    Z     = length(xin)
    ## calculate moving average over window
    ## truncating windows for data points at beginning and end
    movAv = map(1:Z) do i
        r = max(1,i-windowsize):min(i+windowsize,Z)
        mean(view(xin,r))
    end
    ## return slow oscillations in col 1 and fast oscillations in col 2
    xout[:,1] .= movAv
    xout[:,2] .= xin .- movAv
    return xout
end
#----------------------------------------------------------------------------

##We define the input and output dimensions for the decomposition
indims  = InDims("Time")
outdims = OutDims("Time", CategoricalAxis("Scale",["Slow","Fast"]))
cube_decomp = mapCube(movingAverage, world_new, indims=indims, outdims=outdims)
#----------------------------------------------------------------------------

# ### For estimating the temperature sensitivities
# 
# The classical $Q_{10}$ estimation could be realized with the following function

function Q10direct(xout_Q10, xout_rb, xin)
    τ, ρ = eachcol(xin)
    ## solve the regression
    b    = cor(τ, ρ)*std(ρ)/std(τ)
    a    = mean(ρ) - b*mean(τ)

    Q10  = exp(b)
    Rb   = exp(a)
    ##The returned Rb is a constant time series
    xout_rb .= Rb
    xout_Q10 .= Q10
end
#----------------------------------------------------------------------------

# For the scale dependent parameter estimation, the function is a bit more complex. And the numbers in the code comment refer to the  supporting online materials in Mahecha et al. (2010)

function Q10SCAPE(xout_Q10, xout_rb, xin)
    ## xin is now a 3D array with dimensions Time x Scale x Variable
    τ_slow = xin[:, 1, 1]
    τ_fast = xin[:, 2, 1]
    ρ_slow = xin[:, 1, 2]
    ρ_fast = xin[:, 2, 2]
    τ      = τ_slow + τ_fast
    ρ      = ρ_slow + ρ_fast

    ## EQ S5
    ## Q10 calculated on fast oscillations only
    d    = cor(τ_fast, ρ_fast)*std(ρ_fast)/std(τ_fast)
    c    = mean(ρ_fast) - d*mean(τ_fast)
    Q10  = exp(d)

    ## EQ S6: Influence of low frequency temperature on Rb
    ρ_sc = (τ_slow .+ mean(τ)) .* d

    ## EQ S7: Time varying estimate for Rb
    ρ_b  = ρ_slow .+ mean(ρ) .- ρ_sc
    Rb_b  = exp.(ρ_b)

    xout_Q10 .= Q10
    xout_rb  .= Rb_b
end
#----------------------------------------------------------------------------

# ### Application of these functions on the prepared cubes
# 
# 

indims_q10 = InDims("Time","Var")
outdims_q10 = OutDims() ## Just a single number, the first output cube
outdims_rb = OutDims("Time") ## The Rb time series, the second output cube
q10_direct, rb_direct = mapCube(Q10direct, world_new, indims=indims_q10, outdims=(outdims_q10, outdims_rb))
#----------------------------------------------------------------------------

# For the SCAPE approach, the parameter estimation on the decomposed appraoch is then

indims_scape = InDims("Time","Scale","Var",filter=ESDL.DAT.AnyMissing())
q10_scape, rb_scape = mapCube(Q10SCAPE,cube_decomp, indims=indims_scape, outdims=(outdims_q10, outdims_rb))
#----------------------------------------------------------------------------

# ### The rest is plotting
# Including some additional analyses, not included in the paper

function plot_robin(titulo, DAT, clbtitle;
        cm = ColorMap(get_cmap("GnBu", 100)),
        dmin=minimum(skipmissing(DAT)),
        dmax=maximum(skipmissing(DAT)),
        spl = (1,1,1),
    )

    ccrs = pyimport_conda("cartopy.crs","cartopy")
    feat = pyimport_conda("cartopy.feature","cartopy")

    DAT = replace(DAT,missing=>NaN)

    #### make new figure
    fig = plt.figure(figsize=[10, 10])

    #### set the projection
    ax = plt.subplot(spl..., projection=ccrs.Robinson())

    #### add title
    plt.title(titulo, fontsize=18)

    #### land and ocean backgrounds
    ax.add_feature(feat.LAND,  color = [0.9, 0.9, 0.9])
    ax.add_feature(feat.OCEAN, color = [0.85, 0.85, 0.85])
    ax.coastlines(resolution = "50m", color = [0, 0, 0], lw = 0.5)

    #### show data
    im = ax.imshow(reverse(DAT', dims = 1), transform = ccrs.PlateCarree(), cmap = cm, vmin = dmin, vmax = dmax)

    #### add colobar
    clb = plt.colorbar(im,
        pad = 0.05,
        shrink = 0.7,
        aspect = 30,
        orientation = "horizontal",
        extend = "max")

    clb.ax.set_title(clbtitle)
    ####plt.show()

    return fig

end
#----------------------------------------------------------------------------

p1 = plot_robin("a) Confounded Parameter Estimation", q10_direct[:,:], "Q10",dmin=1, dmax=2.5);
savefig("../figures/q10_confounded.pdf",orientation="landscape",bbox_inches="tight")
p2 = plot_robin("b) Scale Dependent Parameter Estimation", q10_scape[:,:], "Q10",dmin=1,dmax=2.5);
savefig("../figures/q10_scape.pdf",orientation="landscape",bbox_inches="tight")
#----------------------------------------------------------------------------

##Construct a new cube
ds = concatenateCubes(tair=world_tair, rb=rb_scape)
## And compute the correlation between Air temperature and Base respiration
using Statistics
cor_tair_rb = mapslices(i->cor(eachcol(i)...),ds, dims=("Time","Variable"))

q10_diff = q10_direct - q10_scape
#----------------------------------------------------------------------------

plot_robin("Correlation Tair and Rb",cor_tair_rb[:,:],"Coefficient",dmin=-1.0, dmax=1.0)
plot_robin("Ratio of Q10conv and Q10Scape",q10_diff[:,:],"Ratio",dmin=-1.0,dmax=1.0);
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
