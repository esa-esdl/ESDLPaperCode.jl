import Literate

sources = filter(i->startswith(i,"ESDL case"),readdir(@__DIR__))



function generate_notebooks()
  sources = ["ESDL case study 1 seasonality.jl"] #,"ESDL case study 2 Intrinsic dimension.jl","ESDL case study 3 q10.jl"]
  foreach(sources) do s
    Literate.notebook(joinpath(@__DIR__,s),joinpath(@__DIR__,"..","notebooks"))
  end
end
