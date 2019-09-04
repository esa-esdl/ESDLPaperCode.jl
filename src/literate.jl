import Literate

function generate_notebooks()
  sources = case_study_files
  foreach(sources) do s
    Literate.notebook(joinpath(@__DIR__,s),joinpath(@__DIR__,"..","notebooks"), execute=false)
  end
end

function generate_markdown()
  sources = ["ESDL case study 1 seasonality.jl"]#,"ESDL case study 2 Intrinsic dimension.jl","ESDL case study 3 q10.jl"]
  foreach(sources) do s
    Literate.markdown(joinpath(@__DIR__,s),joinpath(@__DIR__,"..","docs","src"))
  end
end
