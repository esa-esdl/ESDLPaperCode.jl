# ## Case study 4: Extracting a polygon and aggregate data split by a polygon mask
# ###
# 
# #### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)
# 
# * Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).
# 
# * The NB is written based on Julia 1.3
# 
# * Normal text are explanations referring to notation and equations in the paper
# 
# * `# comments in the code are intended explain specific aspects of the coding`
# 
# * ### New steps in workflows are introduced with bold headers
# 
# Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

## for operating the Earth system data lab
using ESDL, WeightedOnlineStats
#----------------------------------------------------------------------------

# ### Select and subet an Earth system data cube
# 
# We need to choose a cube and here select a 8-dayily, 0.25° resolution global cube. The cube name suggests it is chunked such that we have one time chunk and 720x1440 spatial chunks

cube_handle = Cube("../data/subcube")
#----------------------------------------------------------------------------

# Here we define two subcubes for Gross primary productivity and for surface moisture

gpp = subsetcube(cube_handle, variable = "gross_primary_productivity", time = 2003:2012)
moisture = subsetcube(cube_handle, variable = "surface_moisture", time = 2003:2012)
#----------------------------------------------------------------------------

# The objective is to estimate histograms of gross_primary_productivity and surface moisture and split them by AR5 region. We first download a shapefile defining these regions.

cd("../data") do 
    if !isfile("referenceRegions.shp")
        p = download("https://www.ipcc-data.org/documents/ar5/regions/referenceRegions.zip")
        run(`unzip $p`)
    end
end
#----------------------------------------------------------------------------

# After this we can use the shapefile and apply a rasterization method to convert it to a cube. The `labelsym` argument specifies which field to transfer to the cubes metadata. 

srex = cubefromshape("../data/referenceRegions.shp",gpp,labelsym=:LAB)
using ESDLPlots
plotMAP(srex)
#----------------------------------------------------------------------------

# In order to compute some aggregate statistics over our datasets we join the 3 data cubes into a single iterable table. The data is not loaded but can be iterated over in an efficient manner which is chunk-aware. Additionally we need the latitude values of the Table to compute the weights of our aggregation which represent the grid cell size. 

t = CubeTable(gpp = gpp, moisture=moisture, region=srex, include_axes=("lat",))
#----------------------------------------------------------------------------

# If the concept of this table is still a bit opaque, we demonstrate this by converting a small part of the table to a DataFrame. We just apply a filter to sort out missings and then take 10 values of the Table. 

using DataFrames, Base.Iterators
DataFrame(take(Iterators.filter(r->!any(ismissing,(r.gpp,r.moisture,r.region)),t),10))
#----------------------------------------------------------------------------

# Now comes the actual aggregation. First we generate an empty `WeightedHist` for every SREX region. Then we loop through all the entries in our table and fit the gpp/moisture pair into the respective histogram. In the end we create a new (in-memory) data cube from the resulting histograms. 

using ProgressMeter
function aggregate_by_mask(t,labels)
    hists = [WeightedHist((0.0:1:12,0:0.1:1)) for i=1:33]
    @showprogress for row in t
        if !any(ismissing,(row.gpp, row.moisture, row.region))
            h = hists[row.region[]]
            fit!(h,(row.gpp,row.moisture),cosd(row.lat))
        end
    end
    ##We create the axes for the new output data cube
    midpointsgpp   = 0.5:1.0:11.5
    midpointsmoist = 0.05:0.1:0.95
    newaxes = CubeAxis[
        CategoricalAxis("SREX",[labels[i] for i in 1:33]),
        RangeAxis("GPP",midpointsgpp),
        RangeAxis("Moisture",midpointsmoist),
    ]
    ## And create the new cube object
    data = [WeightedOnlineStats.pdf(hists[reg],(g,m)) for reg in 1:33, g in midpointsgpp, m in midpointsmoist]
    CubeMem(newaxes,data)
end
r = aggregate_by_mask(t,srex.properties["labels"])
#----------------------------------------------------------------------------

# And we create some heatmap plots:

for reg in ("AMZ","CEU","EAF")
    data = r[srex=reg][:,:]
    p = heatmap(0.5:1.0:11.5,0.05:0.1:0.95,(data')./maximum(data), lw = 0, xlab = "gross primary productivity", ylab="surface moisture", title=reg) 
    Plots.savefig(p,"../figures/heatmap_$reg.png")
end
#----------------------------------------------------------------------------

# *This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
