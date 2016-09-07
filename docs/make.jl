using Documenter, ODE

makedocs(
         format = Documenter.Formats.HTML,
         sitename = "ODE",
         pages = [
                  "Home" => "index.md",

                  "Manual" => [ "Basics" => "man/basics.md",
                                "Base" => "man/base.md" ]
                  ]

         )

deploydocs(
           repo = "github.com/JuliaODE/ODE.jl.git",
           deps = Deps.pip("pygments", "python-markdown-math")
           )
