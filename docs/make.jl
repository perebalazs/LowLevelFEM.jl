using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter
using LowLevelFEM

const DOC_PAGES = [
    "Introduction" => "index.md",
    "API Reference" => [
        "General" => "General.md",
        "Linear" => "Linear.md",
        "Heat" => "Heat.md",
        "Nonlinear" => "Nonlinear.md",
        "Operators" => "Operators.md",
    ],
    "Examples" => "Examples.md",
    "News & Updates" => "news/index.md",
]

# HTML dokumentáció
makedocs(
    sitename = "LowLevelFEM",
    authors = "Balázs Pere",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
        sidebar_sitename = true,
        # thresholds are in BYTES (integers)
        size_threshold_warn = 300_000,  # ~293 KiB
        size_threshold = 600_000,       # ~586 KiB
        assets = ["assets/analytics.html"],
    ),
    pages = DOC_PAGES,
    doctest = false,
)

# LaTeX (PDF) dokumentáció
#@info "STARTING LATEX BUILD"
#makedocs(
#    sitename = "LowLevelFEM",
#    authors = "Balázs Pere",
#    format = Documenter.LaTeX(),
#    pages = DOC_PAGES,
#    doctest = false,
#)

try
    include(joinpath(@__DIR__, "src/news/generate_feed.jl"))
catch err
    @warn "RSS feed generation failed" err
end

# --- Copy RSS feed into the build directory ---
try
    src_feed = joinpath(@__DIR__, "src/news/feed.xml")
    dest_dir = joinpath(@__DIR__, "build/news")
    mkpath(dest_dir)
    cp(src_feed, joinpath(dest_dir, "feed.xml"); force=true)
    @info "Copied feed.xml to build/news/"
catch err
    @warn "Could not copy feed.xml to build directory" err
end
# ----------------------------------------------

deploydocs(
    repo = "github.com/perebalazs/LowLevelFEM.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    push_preview = true
)
