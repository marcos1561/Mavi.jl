using Documenter, Mavi

makedocs(
    sitename="Mavi Documentation", 
    repo=Remotes.GitHub("marcos1561", "Mavi.jl"),
    pages = [
        "index.md",
        "manual.md",
        "Blog" => [
            "_posts/2025-11-11_rings.md",
            # "_posts/2025-10-05-variable_particle_number.md",
            "_posts/2025-09-14-space_system.md",
        ]
    ]
    
)

deploydocs(
    repo = "github.com/marcos1561/Mavi.jl.git",
)