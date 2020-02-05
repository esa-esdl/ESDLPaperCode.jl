```@meta
EditURL = "<unknown>/../src/ESDL case study 1 seasonality.jl"
```

## Case study 1: Seasonal dynamics on the land surface
### Based on simple seasonal statistics

#### Miguel D. Mahecha, Fabian Gans et al. (correspondence to: mmahecha@bgc-jena.mpg.de and fgans@bgc-jena.mpg.de)

* Notebook to reproduce and understand examples in the paper *Earth system data cubes unravel global multivariate dynamics* (sub.).

* The NB is written based on Julia 1.1

* Normal text are explanations referring to notation and equations in the paper

* `# comments in the code are itended explain specific aspects of the coding`

* ### New steps in workflows are introduced with bold headers

Sept 2019, Max Planck Institute for Biogeochemistry, Jena, Germany

### Load required packages

```@example ESDL case study 1 seasonality
# for plotting later on (need to be loaded first, to avoid conflicts)
using PyCall, PyPlot, PlotUtils

# for operating the Earth system data lab
using ESDL
```

Next we get a handle to the Earth System Data Cube we want to use, which provides a description of the cube:

```@example ESDL case study 1 seasonality
cube_handle = Cube("../data/subcube")
```

Check which variables are avaiable in the data cube:

```@example ESDL case study 1 seasonality
cubeinfo(cube_handle)
```

```@example ESDL case study 1 seasonality
# if we want the names of the variables:
println(getAxis("Var", cube_handle).values)
```

We decide which variables to plot:

```@example ESDL case study 1 seasonality
vars = ["gross_primary_productivity", "air_temperature_2m", "surface_moisture"]
```

So we "virtually get" the cube data virually:

```@example ESDL case study 1 seasonality
cube_subset = subsetcube(cube_handle, variable = vars)
```

The next function estimates the median seasonal cycle. This changes the dimension of the cube, as the time domain is replaced by day of year (doy); Eq. 9 in the manuscript:
```math
    f_{\{time\}}^{\{doy\}} : \mathcal{C}(\{lat, lon, time, var\}) \rightarrow \mathcal{C}(\{lat, lon, doy, var\})
```

```@example ESDL case study 1 seasonality
# median seasonal cycle built-in function
cube_msc = getMedSC(cube_subset)
```

The resulting cube `cube_msc` has is of the form $\mathcal{C}(\{lat, lon, doy, var\})$. On this cube we want to apply function `nan_med` (see below) to estimate latitudinal averages for all variables. The atomic function (Eq. 10) needs to have the form, i.e. expecting a longitude and returning a scalar:
```math
    f_{\{lon\}}^{\{\}} : \mathcal{C}(\{lat, lon, doy, var\}) \rightarrow \mathcal{C}(\{lat, doy, var\})
```

```@example ESDL case study 1 seasonality
import Statistics.median

# Applied to the dimension "Lon"
cube_msc_lat = mapslices(median ∘ skipmissing, cube_msc, dims = "Lon")
```

The result of each operation on a data cube is a data cube. Here the resulting cube has the form $\mathcal{C}(\{doy, lat, var\})$
as expected but in different order, which is, irrelevant as axes have no natural order.

### Visualization

At this point we leave the `ESDL` and and go for visualizations. Using PyPlot we can generate fig. 3 of the paper, the data can be exatracted from the resutls cube via array-access `A[:, :]`.

```@example ESDL case study 1 seasonality
# Create a plot and make it a polar plot
function zonal_polar_plot(d_msc_lat, sbp, it, vari, lab)

    ax = subplot(sbp, polar = "true")
    ourcmap = ColorMap(get_cmap("Spectral_r", 720))

    # set polar ticks and labels
    month_ang = 0:11
    month_ang = month_ang .* (360/12)
    month_lab = ["Jan"; "Feb"; "Mar"; "Apr"; "May"; "Jun"; "Jul"; "Aug"; "Sep"; "Oct"; "Nov"; "Dec"]

    # tuning
    ax.set_yticklabels([])
    ax.set_thetagrids(angles = month_ang,
                      labels = month_lab,
                      rotation = month_ang)

    # set Jan to the top of the plot
    ax.set_theta_zero_location("N")

    # switch to clockwise
    ax.set_theta_direction(-1)

    ######## color setup
    if isequal(vari, "gross_primary_productivity")
        N_min = 0
        N_max = 8
        N_var = 8
    elseif isequal(vari, "air_temperature_2m")
        N_min = -35
        N_max = 35
        N_var = 15
    elseif isequal(vari, "surface_moisture")
        N_min = 0
        N_max = 50
        N_var = 10
    end

    ######## background setup

    # define time (all could be done more elegantly in julia I guess... )
    time_ang = range(0, stop = 12, step = 0.25)
    N_time   = length(time_ang)
    time_ang = time_ang .* (360/12)
    time_rad = time_ang .* (pi / 180)

    # create a continous var for the background
    N_div  = N_time*10
    y      = range(N_min, stop = N_max, length = N_div)

    # grid to fill the polar plot
    xgrid  = repeat(time_rad', N_div, 1)
    ygrid  = repeat(y, 1, N_time)

    # a grid of NaNs to make an extra colorbar later
    nangrid  = zeros(size(ygrid)).*NaN
    levels   = range(-90, stop = 90, step = 10)
    ticks    = range(-80, stop = 80, step = 20)

    axsurf2 = ax.contourf(xgrid, ygrid, ygrid.*NaN, N_max,
                          cmap = ourcmap,
                          levels = levels,
                          ticks = ticks)

    # background to the range of values
    if isequal(vari, "gross_primary_productivity")
        axsurf = ax.contourf(xgrid, ygrid, ygrid, N_max,
                         cmap = ColorMap("gray_r"),
                         extend = "max")
    elseif isequal(vari, "air_temperature_2m")
        axsurf = ax.contourf(xgrid, ygrid, ygrid, N_max,
                         cmap = ColorMap("gray_r"),
                         extend = "both")
    elseif isequal(vari, "surface_moisture")
        axsurf = ax.contourf(xgrid, ygrid, ygrid, N_max,
                         cmap = ColorMap("gray_r"),
                         extend = "max")
    end

    #### colorbar setup
    if isodd(sbp)
        # add forground colorbar
        cban = colorbar(axsurf2, fraction = 0.05, shrink = 0.5, pad = 0.18)
        cban.ax.set_title(label = "Latitude")
        cban.set_ticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])
        cban.set_ticklabels(["80°N", "60°N", "40°N", "20°N", "0°", "20°S", "40°S", "60°S", "80°S"])
        cban.ax.invert_yaxis()
    else
        # add background colorbar
        cbsurf = colorbar(axsurf, fraction = 0.05, shrink = 0.5, pad = 0.18)
        if  vari == "gross_primary_productivity"
            cbsurf.ax.set_title(label = "GPP")
            cbsurf.set_label(label = "g C / (m2 d)", rotation=270, labelpad=+20)
        elseif vari == "air_temperature_2m"
            cbsurf.ax.set_title(label = "Tair")
            cbsurf.set_label(label = "°C", rotation = 270, labelpad=+20)
        elseif vari == "surface_moisture"
            cbsurf.ax.set_title(label = "Surf. Moist.")
            cbsurf.set_label(label = "[]", rotation = 270, labelpad=+20)
            cbsurf.set_ticks([0, 10, 20, 30, 40, 50])
            cbsurf.set_ticklabels(["0", "0.1", "0.2", "0.3", "0.4", "0.5"])
        end
    end

    ######## forground setup

    # plot the real data
    N_msc = size(d_msc_lat)[1]
    time_ang_dat = range(1/N_msc, step = 1/N_msc, length = N_msc)
    time_ang_dat = time_ang_dat .* (360)
    time_rad_dat = time_ang_dat .* (pi / 180)
    time_rad_dat = [time_rad_dat; time_rad_dat[1]]

    # add your data
    for j = it
      jj = convert(Int, j)
      try
      var_idx = findall(vari .== getAxis(VariableAxis, d_msc_lat).values)[1]
      ts = d_msc_lat[:, jj, var_idx]
      va = [ts; ts[1]]
      # correction for temperature
      if isequal(vari, "air_temperature_2m")
        va = va.-273.15
      elseif isequal(vari, "surface_moisture")
        va = va.*100
      end
      p  = ax.plot(time_rad_dat, va,
      color = ourcmap(jj),
      linewidth = 0.8)
      catch
      end
    end

    ##ax.annotate(lab,
     #       xy = (thistheta, thisr),  ## theta, radius
      #      xytext=(0.05, 0.05),    ## fraction, fraction
       #     textcoords='figure fraction',
        #    arrowprops=dict(facecolor='black', shrink=0.05),
         #   horizontalalignment='left',
          #  verticalalignment='bottom',
           # )
end
```

```@example ESDL case study 1 seasonality
# create a new figure
figure("polar_lineplot", figsize = (10, 15))

# get the latitude values for which we have data
L = collect(getAxis(LatAxis, caxes(cube_msc_lat)).values)

sbps = 321:2:332
labtoshow = ["a)", "b)", "c)", "d)", "e)", "f)"]

for (sbp, lab, vari) in zip(sbps,labtoshow,getAxis(VariableAxis, caxes(cube_msc_lat)).values)

  it1 = range(720/2, stop = 1, step = -2)
  it2 = range(720/2+1, stop = 720, step = 2)

  zonal_polar_plot(cube_msc_lat, sbp, it1, vari, lab)
  zonal_polar_plot(cube_msc_lat, sbp+1, it2, vari, lab)

end

mkpath("../figures")
savefig("../figures/zonalmeans.png",
        bbox_inches = "tight");
nothing #hide
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

