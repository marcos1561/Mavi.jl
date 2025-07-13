"Calculation of thermodynamic quantities."
module Quantities

using Mavi.States: State
using Mavi.Systems: System
using Mavi.Configs

export kinetic_energy, potential_energy

"Return system's kinetic energy."
function kinetic_energy(state::State)
    return sum(state.vel .^2)/2
end

"Return system's potential energy."
function potential_energy(system::System, dynamic_cfg::HarmTruncCfg)
    # Aliases
    N = system.num_p
    r = system.dists
    ko = dynamic_cfg.ko
    ro = dynamic_cfg.ro
    ra = dynamic_cfg.ra

    # Calculation
    pot = 0.0
    for i in 1:N
        for j in i+1:N
            if r[i,j] <= ra
                pot += (r[i,j]-ro)^2
            end
        end
    end
    pot *= ko/2
    return pot
end

function potential_energy(system::System, dynamic_cfg::LenJonesCfg)
    # Aliases
    N = system.num_p
    r = system.dists
    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Calculation
    pot = 0.0
    for i in 1:N
        for j in i+1:N
            pot += ((sigma/r[i,j])^12 - (sigma/r[i,j])^6)
        end
    end
    pot *= 4*epsilon

    return pot
end

end