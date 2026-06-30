"Calculation of thermodynamic quantities."
module Quantities

using Mavi.States
using Mavi.Systems
using Mavi.Integration
using Mavi.Configs
using Mavi.ChunksMod

export kinetic_energy, potential_energy

kinetic_energy(system) = kinetic_energy(system.state)
potential_energy(system) = potential_energy(system, system.dynamic_cfg)

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

function potential(dynamic_cfg::SzaboCfg, dr_norm)
    r_eq = dynamic_cfg.r_eq
    r_max = dynamic_cfg.r_max
    k_rep = dynamic_cfg.k_rep
    k_adh = dynamic_cfg.k_adh

    u = 0
    if dr_norm < r_eq
        u = r_eq * k_rep / 2 * (dr_norm / r_eq - 1)^2 
    elseif dr_norm < r_max
        u = r_eq * k_adh / 2 * (dr_norm / r_eq - 1)^2 
    end

    return u
end

potential_energy(system, dynamic_cfg::SzaboCfg) = potential_energy(system, dynamic_cfg, system.chunks)

function potential_energy(system, dynamic_cfg::SzaboCfg, chunks::ChunksMod.Chunks)
    space_cfg = system.space_cfg
    pos = system.state.pos
    
    total_u = 0
    function pair_potential(i, j, chunk_id)
        dr = Integration.calc_diff(pos[i], pos[j], space_cfg)
        dr_norm = sqrt(sum(dr.^2))
        total_u += potential(dynamic_cfg, dr_norm)
    end
    iterate_over_pais(chunks, pair_potential)

    return total_u
end

function potential_energy(system, dynamic_cfg::SzaboCfg, chunks::Nothing)
    pos = system.state.pos
    space_cfg = system.space_cfg
    
    N = length(pos)
    pot = 0
    for i in 1:N
        for j in i+1:N
            dr = calc_diff(pos[i], pos[j], space_cfg)
            dr_norm = sqrt(sum(dr.^2))
            pot += potential(dynamic_cfg, dr_norm)
        end
    end

    return pot
end

end