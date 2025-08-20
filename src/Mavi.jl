module Mavi

module Utils
include("utils.jl")
end

# Core Stuff
include("configs.jl")
include("states.jl")
include("init_states.jl")
include("space_checks.jl")
include("chunks.jl")
include("systems.jl")
include("integration.jl")
include("quantities.jl")
include("serder.jl")

using .Systems
export System 

using .MaviSerder
export save, load_system

using .States

# Special Systems 
include("rings/rings.jl")

include("visualization.jl")

end
