using Mavi: Mavi
using Mavi.Integration
using Mavi.Configs
using Mavi.Info

state = Mavi.State{Float64}(
    pos = [1 -2.5 3.3 -4 5; -1.7 2.1 -3.8 4.4 -5.4],
    vel = [0.3 2 5.7 9.8 3.; 1 0 7.8 .12 2.2],
)

system = Mavi.System(
    state=state, 
    space_cfg=CircleCfg(radius=10), 
    dynamic_cfg=LenJonesCfg(sigma=2,epsilon=4),
    int_cfg=IntegrationCfg(dt=0.01),
)

# Make one step to first calculte distances
Integration.step!(system)
# Time evolution
for i in 1:100
    Integration.step!(system)
    if i % 5 == 0 # print every 5 steps
        pe = potential_energy(system,system.dynamic_cfg)
        ke = kinetic_energy(state)
        println("Potential: $pe, kinetic: $ke, total: $(pe+ke)")
    end
end