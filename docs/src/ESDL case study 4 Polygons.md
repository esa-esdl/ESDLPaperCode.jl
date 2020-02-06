```@meta
EditURL = "<unknown>/src/ESDL case study 4 Polygons.jl"
```

## Case study 4: Extracting a polygon and aggregate data split by a polygon mask
###

#### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)

* Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).

* The NB is written based on Julia 1.3

* Normal text are explanations referring to notation and equations in the paper

* `# comments in the code are intended explain specific aspects of the coding`

* ### New steps in workflows are introduced with bold headers

Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

```@example ESDL case study 4 Polygons
# for operating the Earth system data lab
using ESDL, WeightedOnlineStats
```

### Load the pre-downloaded Earth-System Datacube

```@example ESDL case study 4 Polygons
cube_handle = Cube("../data/subcube")
```

Here we define two subcubes for Gross primary productivity and for surface moisture

```@example ESDL case study 4 Polygons
gpp = subsetcube(cube_handle, variable = "gross_primary_productivity", time = 2003:2012)
moisture = subsetcube(cube_handle, variable = "surface_moisture", time = 2003:2012)
```

The objective is to estimate histograms of gross_primary_productivity and surface moisture and split them by AR5 region. We first download a shapefile defining these regions.

```@example ESDL case study 4 Polygons
cd("../data") do
    if !isfile("referenceRegions.shp")
        p = download("https://www.ipcc-data.org/documents/ar5/regions/referenceRegions.zip")
        run(`unzip $p`)
    end
end
```

After this we can use the shapefile and apply a rasterization method to convert it to a cube. The `labelsym` argument specifies which field to transfer to the cubes metadata.

```@example ESDL case study 4 Polygons
srex = cubefromshape("../data/referenceRegions.shp",gpp,labelsym=:LAB)
```

In order to compute some aggregate statistics over our datasets we join the 3 data cubes into a single iterable table. The data is not loaded but can be iterated over in an efficient manner which is chunk-aware. Additionally we need the latitude values of the Table to compute the weights of our aggregation which represent the grid cell size.

```@example ESDL case study 4 Polygons
t = CubeTable(gpp = gpp, moisture=moisture, region=srex, include_axes=("lat",))
```

If the concept of this table is still a bit opaque, we demonstrate this by converting a small part of the table to a DataFrame. We just apply a filter to sort out missings and then take 10 values of the Table.

```@example ESDL case study 4 Polygons
using DataFrames, Base.Iterators
DataFrame(take(Iterators.filter(r->!any(ismissing,(r.gpp,r.moisture,r.region)),t),10))
```

Now comes the actual aggregation. First we generate an empty `WeightedHist` for every SREX region. Then we loop through all the entries in our table and fit the gpp/moisture pair into the respective histogram. Never will the whole cube be loaded into memory, but only one chunk is read at a time. In the end we create a new (in-memory) data cube from the resulting histograms.

```@example ESDL case study 4 Polygons
using ProgressMeter
function aggregate_by_mask(t,labels)
    n_classes = length(labels)
    # Here we create an empty 2d histogram for every SREX region

    ####hists = [WeightedHist((0.0:1:12,0:0.1:1)) for i=1:n_labels]
    hists = [WeightedHist((0.0:0.1:12,0:0.01:1)) for i=1:n_classes]

    # Now loop through every data point (in space and time)
    @showprogress for row in t
        # If all data are there
        if !any(ismissing,(row.gpp, row.moisture, row.region))
            ####We select the appropriate histogram according to the region the data point belongs to
            h = hists[row.region[]]
            ####And we fit the two data points to the histogram, weight by cos of lat
            fit!(h,(row.gpp,row.moisture),cosd(row.lat))
        end
    end
    ########We create the axes for the new output data cube
    midpointsgpp   = 0.05:0.1:11.95
    midpointsmoist = 0.005:0.01:0.995
    newaxes = CubeAxis[
        CategoricalAxis("SREX",[labels[i] for i in 1:33]),
        RangeAxis("GPP",midpointsgpp),
        RangeAxis("Moisture",midpointsmoist),
    ]
    # And create the new cube object
    data = [WeightedOnlineStats.pdf(hists[reg],(g,m)) for reg in 1:33, g in midpointsgpp, m in midpointsmoist]
    CubeMem(newaxes,data)
end
r = aggregate_by_mask(t,srex.properties["labels"])
```

To illustrate the output we plot the density for the region "Eastern Africa":

```@example ESDL case study 4 Polygons
import Plots
Plots.heatmap(0.005:0.01:0.995,0.05:0.1:11.95,r[srex="EAF"][:,:], clim=(0,5e-7), xlabel="Moisture", ylabel="GPP")
```

```@example ESDL case study 4 Polygons
saveCube(r, "../data/srex_aggregate.zarr", overwrite = true)
```

## Python plotting

To generate the publication-quality plots we use python plotting tools with the following code, which does not demonstrate any ESDL capabilities but is included here for reproducbility:

```@example ESDL case study 4 Polygons
# for plotting
using PyCall, PyPlot, PlotUtils
```

```@example ESDL case study 4 Polygons
lat = getAxis("lat", srex).values
lon = getAxis("lon", srex).values

cm = ColorMap(get_cmap("gray", 33))
ccrs = pyimport_conda("cartopy.crs","cartopy")
feat = pyimport_conda("cartopy.feature","cartopy")

######## make new figure
fig = plt.figure(figsize=[10, 10])

ax = subplot(313, )


######## set the projection
ax = plt.subplot(projection=ccrs.Robinson())

######## add title
plt.title("IPCC AR5 regions", fontsize=20)

######## land and ocean backgrounds
ax.add_feature(feat.LAND,  color = [0.9, 0.9, 0.9])
ax.add_feature(feat.OCEAN, color = [0.85, 0.85, 0.85])
ax.coastlines(resolution = "50m", color = [0, 0, 0], lw = 0.5)

######## show data
DAT = srex[:, :]
DAT = replace(DAT, missing=>NaN)
DAT2 = reverse(DAT', dims = 1)

im = ax.imshow(DAT2, transform = ccrs.PlateCarree(), cmap = cm, vmin = 0, vmax = 30)

for i in 1:33
    if  !(srex.properties["labels"][i] in ["NTP*", "ETP*", "STP*"])
        idx_lat = mean(lat[j.I[2]] for j in CartesianIndices(DAT) if DAT[j]==i)
        idx_lon = mean(lon[j.I[1]] for j in CartesianIndices(DAT) if DAT[j]==i)
        ax.text(idx_lon, idx_lat, weight="bold", color = "red", srex.properties["labels"][i], horizontalalignment = "center", verticalalignment = "center", transform=ccrs.Geodetic())
    end
end

####savefig("SREX.pdf", bbox_inches = "tight")
savefig("../figures/SREX.png", bbox_inches = "tight")
```

And we create some heatmap plots:

```@example ESDL case study 4 Polygons
gpp = collect(getAxis("GPP", r).values)
moi = collect(getAxis("Moisture", r).values);
nothing #hide
```

```@example ESDL case study 4 Polygons
# Draw a heatmap with the numeric values in each cell
####f, ax = plt.subplots(figsize=(9, 6))
sns.set_style("white")
sns.axes_style("darkgrid")

plt.ioff()

for i in 1:33
    reg = srex.properties["labels"][i]
    figure(figsize = (5, 5))
    data = r[srex=reg][:,:]
    data = sqrt.(replace((data')./maximum(data),0.0=>NaN))

    data = transpose(data)

    ax = sns.heatmap(data,  cmap = sns.cm.rocket_r)#, vmin = 0.1, vmax=1)
    ax.set_yticks(1:10:120)
    ax.set_yticklabels(round.(gpp[1:10:120]; digits = 0))

    ax.set_xticks(1:10:100)
    ax.set_xticklabels(round.(moi[1:10:100]; digits = 1))
    #ax.set_xlabel("Surface Moisture []")
    #ax.set_ylabel("Gross Primary Production [g C m⁻² d⁻¹]")
    ax.invert_yaxis()
    ax.set_title(reg, fontsize = 20)
    savefig("normal" * reg * ".png", bbox_inches = "tight")
end
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

