# Physical Quantities
The `Quantities` module contains functions that return physical quantities of interest. Currently, two of them are defined:

- **Kinetic energy:**  
    The system's kinetic energy is calculated by the function `kinetic_energy(state)`, which takes the system's state as a parameter. Since it depends only on the velocities, it is independent of the chosen dynamic configuration.

- **Potential energy:**  
    The system's potential energy is calculated by the function `potential_energy(system, dynamic_cfg)`, which takes the system and the dynamic configuration as arguments, since it depends on the chosen interaction between particles.  
    The user should ensure that the correct dispatch for the `potential_energy()` function is defined for any custom dynamic configuration.

The example [print_energy.jl](https://github.com/marcos1561/Mavi.jl/blob/main/examples/print_energy.jl) uses these functions.