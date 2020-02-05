using Documenter

makedocs(
  format = Documenter.HTML(),
  sitename = "ESDLPaperCode.jl",
  authors = "Miguel Mahecha, Fabian Gans",
  pages = [
    #"Home" => "index.md",
    "Case studies" => [
      "Seasonality" => "ESDL case study 1 seasonality.md",
      "Intrinsic Dimensions" => "ESDL case study 2 Intrinsic dimension.md",
      "Q10" => "ESDL case study 3 q10.md",
      "Polygons" => "ESDL case study 4 Polygons.md",
      ]
  ],
  workdir = @__DIR__
)
