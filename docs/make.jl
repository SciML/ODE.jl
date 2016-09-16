using Documenter, ODE

makedocs(
         format = Documenter.Formats.HTML,
         modules = [ODE],
         clean = false,
         sitename = "ODE.jl",
         pages = [
                  "Home" => "index.md",

                  "Manual" => [
                               "Basics" => "man/basics.md",
                               "Base" => "man/base.md"
                               ]
                  ]

         )

deploydocs(
           repo = "github.com/JuliaODE/ODE.jl.git",
           target = "build",
           julia = "0.5",
           deps = nothing,
           make = nothing
           )
