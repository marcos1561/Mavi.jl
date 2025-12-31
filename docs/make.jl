using Documenter, Mavi

makedocs(
    sitename="Mavi Documentation", 
    repo=Remotes.GitHub("marcos1561", "Mavi.jl"),
    pages = [
        "index.md",
        "Manual" => [
            "manual/installation.md",
            "manual/quick_start.md",
            "manual/philosophy.md",
            "manual/expanding_mavi.md",
            "manual/physical_quantities.md",
            "manual/experiments.md",
            "manual/visual_interface.md",
        ],
        "Blog" => [
            "Posts" => [
                "posts/2025-11-11_rings.md",
                "posts/2025-09-14-space_system.md",
            ],
        ],
    ]
    
)

deploydocs(
    repo = "github.com/marcos1561/Mavi.jl.git",
)