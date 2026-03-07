using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter
using LowLevelFEM

const DOC_PAGES = [
    "Home" => "index.md",
    "Getting Started" => [
        "Overview" => "getting-started/index.md",
        "Installation" => "getting-started/installation.md",
        "First Problem" => "getting-started/first-problem.md",
        "Mesh and Physical Groups" => "getting-started/mesh-and-physical-groups.md",
        "Fields and Results" => "getting-started/fields-and-results.md",
    ],
    "Tutorials" => [
        "Overview" => "tutorials/index.md",
        "Linear Elasticity (2D)" => "tutorials/linear-elasticity-2d.md",
        "Heat Conduction (Transient)" => "tutorials/heat-conduction-transient.md",
        "Modal and Buckling" => "tutorials/modal-and-buckling.md",
        "Nonlinear Large Deformation" => "tutorials/nonlinear-large-deformation.md",
        "Multifield Weak-Form DSL" => "tutorials/multifield-weak-form-dsl.md",
        "Poisson and Custom Operators" => "tutorials/poisson-and-custom-operators.md",
        "Legacy Example Gallery" => "tutorials/legacy-examples.md",
    ],
    "Manual" => [
        "Overview" => "manual/index.md",
        "Core Types" => "manual/core-types.md",
        "Boundary Conditions and Loads" => "manual/boundary-and-loads.md",
        "Assembly and Solvers" => "manual/assembly-and-solvers.md",
        "Postprocessing and Visualization" => "manual/postprocessing-and-visualization.md",
        "Operators and Fields" => "manual/operators-and-fields.md",
        "Coordinate Systems" => "manual/coordinate-systems.md",
        "Performance and Parallel" => "manual/performance-and-parallel.md",
        "Troubleshooting" => "manual/troubleshooting.md",
    ],
    "API Reference" => [
        "Overview" => "reference/index.md",
        "Preprocessing" => "reference/preprocessing.md",
        "Fields" => "reference/fields.md",
        "Operators" => "reference/operators.md",
        "Multifield" => "reference/multifield.md",
        "Linear" => "reference/linear.md",
        "Nonlinear" => "reference/nonlinear.md",
        "Heat" => "reference/heat.md",
        "Poisson (Legacy Single-Field)" => "reference/poisson.md",
        "Postprocessing" => "reference/postprocessing.md",
        "Extra" => "reference/extra.md",
    ],
    "Explanations" => [
        "Overview" => "explanations/index.md",
        "Weak-Form DSL Design" => "explanations/weak-form-dsl-design.md",
        "Nodal vs Element Fields" => "explanations/nodal-vs-element-fields.md",
        "Matrix-Level Workflow" => "explanations/matrix-level-workflow.md",
    ],
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
        assets = [
         "assets/plausible.js",
        ]

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
