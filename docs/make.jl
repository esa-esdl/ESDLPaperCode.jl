using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using Documenter

makedocs(
  format = Documenter.HTML(),
  sitename = "ESDLExampleCode.jl",
  authors = "Miguel Mahecha, Fabian Gans",
  pages = [
    "Home" => "index.md",
    "Case studies" => [
      "Seasonality" => "ESDL case study 1 seasonality.md",
      "Intrinsic Dimensions" => "ESDL case study 2 intrinsic dimension.md",
      "Q10" => "ESDL case study 3 q10.md",
      ]
  ],
  workdir = @__DIR__
)
