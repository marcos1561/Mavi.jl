"Calculation of thermodynamic quantities."
module Quantities

using Mavi.States
using Mavi.Systems
using Mavi.Integration
using Mavi.Configs

export kinetic_energy, potential_energy

"Return system's kinetic energy."
function kinetic_energy(state::State)
    ke = 0
    for vel in state.vel
        ke += sum(abs2, vel)
    end
    return ke / 2
end

"Return system's potential energy."
function potential_energy(system::System, dynamic_cfg::HarmTruncCfg)
    # Aliases
    pos = system.state.pos
    space_cfg = system.space_cfg
    N = get_num_total_particles(system)
    ko = dynamic_cfg.ko
    ro = dynamic_cfg.ro
    ra = dynamic_cfg.ra

    # Calculation
    pot = 0.0
    for i in 1:N
        for j in i+1:N
            dr = calc_diff(pos[i], pos[j], space_cfg)
            dist = sqrt(sum(dr.^2))

            if dist <= ra
                pot += (dist-ro)^2
            end
        end
    end
    pot *= ko/2
    return pot
end

function potential_energy(system::System, dynamic_cfg::LenJonesCfg)
    # Aliases
    pos = system.state.pos
    space_cfg = system.space_cfg
    N = get_num_total_particles(system)
    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Calculation
    pot = 0.0
    for i in 1:N
        for j in i+1:N
            dr = calc_diff(pos[i], pos[j], space_cfg)
            dist = sqrt(sum(dr.^2))
            pot += ((sigma/dist)^12 - (sigma/dist)^6)
        end
    end
    pot *= 4*epsilon

    return pot
end

end