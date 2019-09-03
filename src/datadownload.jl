using ESDL
import Zarr

#function get_data_nb1()
function download_data()
  c = S3Cube()

  vars = ["air_temperature_2m",
  "surface_moisture",
  "evaporative_stress",
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


  csub = subsetcube(c,time=2000:2017, var=vars)

  destp = joinpath(@__DIR__,"..","data","subcube")

  cs = saveCube(permutedims(csub,(3,1,2,4)),destp,compressor=Zarr.BloscCompressor(clevel=1), overwrite=true, max_cache=1e9)
end
