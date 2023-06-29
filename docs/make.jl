using Documenter
using IRMO

makedocs(
    sitename = "IRMO.jl",
    modules = [IRMO],
    pages = [
            "Home" => "index.md"
            ],
    format = Documenter.HTML()
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/michel-mata/IRMO.jl.git",
    devbranch = "main",
    target = "build",
    versions = ["stable" => "v^"]
)