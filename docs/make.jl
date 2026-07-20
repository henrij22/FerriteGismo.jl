using Documenter, FerriteGismo

const liveserver = "liveserver" in ARGS

if liveserver
    using Revise
    Revise.revise()
end

DocMeta.setdocmeta!(FerriteGismo, :DocTestSetup, :(using FerriteGismo); recursive = true)

makedocs(;
    format = Documenter.HTML(;
        canonical = "https://henrij22.github.io/FerriteGismo.jl/stable",
        collapselevel = 1,
        assets = String[],
    ),
    repo = Documenter.Remotes.GitHub("henrij22", "FerriteGismo.jl"),
    modules = [FerriteGismo],
    sitename = "FerriteGismo.jl",
    warnonly = true, checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Tutorials overview" => "tutorials/index.md",
            "tutorials/heat_equation.md",
            "tutorials/linear_elasticity.md",
        ],
        "API Reference" => "api_reference.md",
        "Developer documentation" => "developer.md",
    ],
)

if !liveserver
    deploydocs(;
        repo = "github.com/henrij22/FerriteGismo.jl.git",
        push_preview = true,
        versions = [
            "stable" => "v^",
            "dev" => "dev",
        ],
    )
end
