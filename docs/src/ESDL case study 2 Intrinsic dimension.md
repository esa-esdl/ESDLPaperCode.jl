```@meta
EditURL = "@__REPO_ROOT_URL__/"
```

## Case study 2: Intrinsic dimensions of ecosystem dynamics
### As estimate based on PCAs

#### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)

* Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).

* The NB is written based on Julia 1.1

* Normal text are explanations referring to notation and equations in the paper

* `# comments in the code are itended explain specific aspects of the coding`

* ### New steps in workflows are introduced with bold headers

Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

### Load required packages

```@example ESDL case study 2 Intrinsic dimension
# for plotting later on (need to be loaded first, to avoid conflicts)
using PyCall, PyPlot

# for operating the Earth system data lab
using ESDL

# for parallel computing
using Distributed
```

In this study we investigate the redundancy the different variables in each pixel. Therefore we calculate a linear dimensionality reduction (PCA) and check how many dimensions are needed to explain 90% of the variance of a cube that contained originally 6 variables.  First we check out the variables from the cube and add some processors, because we want to do a global study

### Select and prepare (subset/gapfill) an Earth system data cube

We need to choose a cube and here select a 8-dayly, 0.25° resolution global cube. The cube name suggests it is chunked such that we have one time chunk and 720x1440 spatial chunks

```@example ESDL case study 2 Intrinsic dimension
cube_handle = Cube("../data/subcube")
```

As we see here, we have a data cube of the form (compare Table 1 in the paper):

```math
    \mathcal{C}(\{lat, lon, time, var\})
```

There is a command that returns the metadata fro the variable axis for better orientation:

```@example ESDL case study 2 Intrinsic dimension
cubeinfo(cube_handle)
```

```@example ESDL case study 2 Intrinsic dimension
# if we want the names of the variables:
println(getAxis("Var", cube_handle).values)
```

Having the variable names allows us to make a selection, such that we can subset the global cube. We should also take care that the variables are as complete as possible in the time window we analyze. This has been explored a priori.

```@example ESDL case study 2 Intrinsic dimension
# vector of variables we will work with
vars = ["evaporative_stress",
    "latent_energy",
    "black_sky_albedo_avhrr",
    "fapar_tip",
    "root_moisture",
    "transpiration",
    "white_sky_albedo_avhrr",
    "sensible_heat",
    "bare_soil_evaporation",
    "net_radiation",
    "net_ecosystem_exchange",
    "evaporation",
    "terrestrial_ecosystem_respiration",
    "land_surface_temperature",
    "leaf_area_index",
    "white_sky_albedo",
    "gross_primary_productivity",
    "black_sky_albedo"];

# time window where most of them are complete
timespan = Date("2003-01-01")..Date("2011-12-31")

# subset the grand cube and get the cube we will analyse here
cube_subset = subsetcube(cube_handle, time = timespan, variable = vars)
```

An important preprocessing step is gapfilling. We do not want to enter the debate on the optimal gapfilling method. What we do here is gapfilling first with the mean seasonal cycle (where it can be estimated), and interpolating long-recurrent gaps (typically in winter seasons).

```@example ESDL case study 2 Intrinsic dimension
# use the ESDL buit-in function
cube_fill = gapFillMSC(cube_subset)
```

The interpolation of wintergpas needs a function that we code here an call `LinInterp`.

```@example ESDL case study 2 Intrinsic dimension
using Interpolations

function LinInterp(y)
  try
    # find the values we need to input
    idx_nan = findall(ismissing, y)
    idx_ok  = findall(!ismissing, y)

    # make sure to have a homogenous input array
    y2 = Float32[y[i] for i in idx_ok]

    # generate an interpolation object based on the good data
    itp = interpolate((idx_ok,), y2, Gridded(Linear()))

    # fill the missing values based on a linter interpolation
    y[idx_nan] = itp(idx_nan)
    return y
  catch
    return y
  end
end


# short test
x = [2.5,missing,3.8,missing,8.9]
LinInterp(x)
```

The function `LiInterp` can now be applied on each time series, so we would have a rather trival mapping of the form:

\begin{equation}
  f_{\{time\}}^{\{time}\} : \mathcal{C}(\{lat, lon, time, var\}) \rightarrow \mathcal{C}(\{lat, lon, time, var\}).
\end{equation}

For operations of this kind, the best is to use the `mapslices` function. In the ESDL package, this function needs the input function, the cube handle, and an indication on which dimension we would apply it. The function can then infer that the output dimension here is also an axis of type `Time`:

```@example ESDL case study 2 Intrinsic dimension
cube_fill_itp = mapslices(LinInterp, cube_fill, dims = "Time")
```

As we describe in the paper, we estimate the intrinsic dimensions from the raw, yet gapfilled, data cube (`cube_fill_itp`), but also based on spectrally decomposed data. The decomposition via discrete FFTs is an atomic operation of the following form (Eq. 12),

\begin{equation}
  f_{\{time\}}^{\{time, freq\}} : \mathcal{C}(\{lat, lon, time, var\}) \rightarrow \mathcal{C}(\{lat, lon, time, var, freq\}).
\end{equation}

which can be done using a pre-implemented ESDL function:

```@example ESDL case study 2 Intrinsic dimension
import Zarr
ESDL.ESDLDefaults.compressor[] = Zarr.BloscCompressor(clevel=1)
cube_decomp = filterTSFFT(cube_fill_itp)
```

### Estimate intrinic dimension via PCA

For estimating the intrinsic estimation via PCA from a multivariate time series we need essenntially two atomic functions. First, dimensionality reduction,

\begin{equation}
     f_{\{time, var\}}^{\{time, princomp \}} : \mathcal{C}(\{time, var\}) \rightarrow \mathcal{C}(\{time, princomp\})
\end{equation}

And second estimating from the reduced space the number of dimensions that represent more variance than the threshold (for details see paper):
\begin{equation}
     f_{\{time, princomp\}}^{\{ \}} : \mathcal{C}(\{time, var\}) \rightarrow \mathcal{C}(\{int dim\})
\end{equation}

However, we as both steps emerge from the same analysis it is more efficient to wrap these two steps in a single atomic functions which has the structure:
\begin{equation}
     f_{\{time, var\}}^{\{ \}} : \mathcal{C}(\{time, var\}) \rightarrow \mathcal{C}(\{\})
\end{equation}

We can now apply this to the cube: The latter was the operation described in the paper (Eq. 11) as

\begin{equation}
     f_{\{time, var\}}^{\{ \}} : \mathcal{C}(\{lat, lon, time, var\}) \rightarrow \mathcal{C}(\{lat, lon\})
\end{equation}

```@example ESDL case study 2 Intrinsic dimension
# packages needed on each core
using MultivariateStats, Statistics

function sufficient_dimensions(xin::AbstractArray, expl_var::Float64 = 0.95)

  any(ismissing,xin) && return NaN
  npoint, nvar = size(xin)
  means = mean(xin, dims = 1)
  stds  = std(xin,  dims = 1)
  xin   = broadcast((y,m,s) -> s>0.0 ? (y-m)/s : one(y), xin, means, stds)
  pca = fit(PCA, xin', pratio = 0.999, method = :svd)
  return findfirst(cumsum(principalvars(pca)) / tprincipalvar(pca) .> expl_var)
end
```

We first apply the function `cube_decomp`to the standard data cube with the threshold of 95% of retained variance. as we see from the description of the atomic function above, we need as minimum input dimension `Time` and `Variable`. We call the output cube `cube_int_dim`, which efficiently is a map.

```@example ESDL case study 2 Intrinsic dimension
cube_int_dim = mapslices(sufficient_dimensions, cube_fill_itp, 0.95, dims = ("Time","Variable"))
```

Saving intermediate results can save CPU later, not needed to guarantee reproducability tough

```@example ESDL case study 2 Intrinsic dimension
saveCube(cube_int_dim, "../data/IntDim", overwrite=true)
```

Now we apply the same function

\begin{equation}
    f_{\{time, var\}}^{\{ \}} : \mathcal{C}(\{time, var\}) \rightarrow \mathcal{C}(\{\})
\end{equation}

to the spectrally decomposed cube (Eq. 13):

\begin{equation}
       f_{\{time, var\}}^{\{\}} : \mathcal{C}(\{lat, lon, time, var, freq\})\rightarrow \mathcal{C}(\{lat, lon, freq\})
\end{equation}

```@example ESDL case study 2 Intrinsic dimension
cube_int_dim_dec = mapslices(sufficient_dimensions, cube_decomp, 0.95, dims = ("Time","Variable"))
```

```@example ESDL case study 2 Intrinsic dimension
saveCube(cube_int_dim_dec, "../data/IntDimDec", overwrite=true)
```

### Visualizing results is not part of the ESDL package.
Here we rely on PyPlot to use the neat `cartopy`pacakge

```@example ESDL case study 2 Intrinsic dimension
ccrs = pyimport_conda("cartopy.crs","cartopy")
feat = pyimport_conda("cartopy.feature","cartopy")

# function to plot global map
# not generic - for this application only
function plot_robin(titulo, DAT, clbtitle)

    # make new figure
    fig = plt.figure(figsize=[10, 10])

    # set the projection
    ax = plt.subplot(1, 1, 1, projection=ccrs.Robinson())

    # add title
    plt.title(titulo, fontsize=18)

    # land and ocean backgrounds
    ax.add_feature(feat.LAND,  color = [0.9, 0.9, 0.9])
    ax.add_feature(feat.OCEAN, color = [0.85, 0.85, 0.85])
    ax.coastlines(resolution = "50m", color = [0, 0, 0], lw = 0.5)

    # show data
    im = ax.imshow(reverse(DAT', dims = 1), transform = ccrs.PlateCarree(), cmap = cm, vmin = 1.5, vmax = 12.5)

    # add colobar
    clb = plt.colorbar(im,
        pad = 0.05,
        shrink = 0.7,
        aspect = 30,
        orientation = "horizontal",
        drawedges = "false",
        extend = "both",
        ticks = 2:12)

    clb.ax.set_title(clbtitle)
    ##plt.show()

    return fig
end
```

We write a little loop to plot and save all subplots of Fig. 4. One key information we need is, in which order the frequencies are saved so we do

```@example ESDL case study 2 Intrinsic dimension
cube_int_dim_dec
```

```@example ESDL case study 2 Intrinsic dimension
println(getAxis("Scale", cube_int_dim_dec).values)
```

```@example ESDL case study 2 Intrinsic dimension
scale_name = getAxis("Scale", cube_int_dim_dec).values
```

```@example ESDL case study 2 Intrinsic dimension
# corresponding colormap
cm = ColorMap(get_cmap("magma_r", 11))

prelab = ["a) ", "b) ", "c) ", "d )"]

# overwrite the first name
scale_name = ["Original Data", "Long-term variability", "Seasonal variability", "Short-term variability"]

save_name = ["Original", "Long", "Annual", "Fast"]

# go through the four subplots
for i in 1:4

    # array access to get the results out of the cube
    if i == 1
        # cube_int_dim has only dimensions lat, lon so
        VAL = cube_int_dim[:, :]
    else
        # cube_int_dim_dec has dimensions lat, lon, freq so
        VAL = cube_int_dim_dec[i, :, :]
    end

    # missings -> NaN for PyPlot
    DAT = zeros(size(VAL))./0.0 ## matrix of NaN
    idx = findall(!ismissing, VAL) ## index if real vals
    DAT[idx] = VAL[idx] ## only insert these

    name = prelab[i]*scale_name[i]
    fig = plot_robin(prelab[i]*scale_name[i], DAT, "Intrinsic DImensions")
    mkpath("../figures")
    savefig("../figures/IntDim_" * save_name[i] * ".pdf",
        orientation = "landscape",
        bbox_inches = "tight")
end
```

```@example ESDL case study 2 Intrinsic dimension
fig  = figure("histplot", figsize = (10, 10))

lableg = ["a) Original Data", "b) Long-term variability", "c) Seasonal variability", "d) Short-term variability"]

lat_ax_vals = getAxis(LatAxis, cube_int_dim.axes).values
lon_ax_vals = getAxis(LonAxis, cube_int_dim.axes).values
weights = cosd.(repeat(lat_ax_vals, inner = length(lon_ax_vals)))

for i = 1:4
    ax   = subplot(410 + i)

    # lon x lat
    data = i == 1 ? cube_int_dim[:, :] : cube_int_dim_dec[i, :, :]

    data = vec(data)
    idx_use =  findall(x -> !ismissing(x) && !isnan(x), data)

    data = data[idx_use]
    latw  = weights[idx_use]

    N, bins, patches = ax.hist(
        data,
        weights = latw,
        density = true,
        alpha = 1,
        bins = 0.5:1:13,
        color = "r",
        linewidth = 0,
        label = lableg[i]
    )

    for j in eachindex(patches)
        patches[j].set_facecolor(cm(max(j-1), 1))
    end

    hist(
        data,
        weights = latw,
        density = true,
        alpha = 1,
        color = [0,0,0,1],
        bins = 0.5:1:13,
        linewidth = 1,
        histtype="step"
    )


    ##legend(frameon= "False",fontsize = 12)
    ax.set_ylim(0, 0.6)
    ax.set_xlim(0.5, 12.5)
    xticks(1:12)
    text(0.02, 0.85, lableg[i],fontsize = 14, transform = ax.transAxes)

    if i == 1
        ylabel("Weighted Frequency", fontsize = 14)
    end
end

xlabel("Intrinsic dimension", fontsize = 14)

savefig("../figures/IntDim_Hist.pdf", bbox_inches = "tight")
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

