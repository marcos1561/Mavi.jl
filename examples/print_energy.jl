"""
Example: Extending the Information UI

Calculates kinetic, potential, and total energies during integration, displaying
the results in real time on the UI.
"""
module Example

using Printf

using Mavi
using Mavi.Integration
using Mavi.Configs
using Mavi.Quantities
using Mavi.Visualization

function main()
    state = Mavi.SecondLawState{Float64}(
        pos = [1 -2.5 3.3 -4 5; -1.7 2.1 -3.8 4.4 -5.4],
        vel = [0.3 2 5.7 9.8 3.; 1 0 7.8 .12 2.2],
    )

    system = Mavi.System(
        state=state, 
        space_cfg=SpaceCfg(
            wall_type=RigidWalls(),
            geometry_cfg=CircleCfg(radius=10),
        ), 
        dynamic_cfg=LenJonesCfg(sigma=2,epsilon=4),
        int_cfg=IntCfg(dt=0.01),
    )

    "Calculates and returns energies formatted for printing"
    function get_energy(system, exec_info)
        energys = []

        pe = potential_energy(system, system.dynamic_cfg)
        ke = kinetic_energy(state)
        
        push!(energys, ("Potencial", @sprintf("%7.3f", pe)))
        push!(energys, ("Kinetic", @sprintf("%7.3f", ke)))
        push!(energys, ("Total", @sprintf("%7.3f", pe + ke)))

        return energys
    end

    anim_cfg = AnimationCfg(
        info_cfg=DefaultInfoUICfg(custom_items=get_energy),
    )

    animate(system, Integration.newton_step!, anim_cfg)
end

end

import .Example
Example.main()