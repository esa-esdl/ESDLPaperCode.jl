module ESDLPaperCode
export run_case_study, run_all

include("datadownload.jl")
include("literate.jl")

const case_study_files = ["ESDL case study 1 seasonality.jl","ESDL case study 2 Intrinsic dimension.jl","ESDL case study 3 q10.jl"]

function run_case_study(i;redownload=false)
  isdir(joinpath(@__DIR__,"..","data","subcube")) || redownload || download_data()
  cd(@__DIR__) do
    include(case_study_files[i])
  end
end

function run_all(;redownload=false)
  foreach(1:3) do i
    run_case_study(i,redownload=redownload)
  end
end


end
