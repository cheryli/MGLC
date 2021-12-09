# Parallel LBM codes for complex flow

This project focus on the multi-GPU parallel implementation of LBM(Lattice Boltzmann Method). We use  OpenACC to accelerate codes on single GPU and MPI for inter-GPU communication.

**Blod** for completed, *italic* for coded but not documented. <u> underlined </u> means we are now working on this.
* `MPI` (Pure) MPI implementation for multi processor, Including communication overlap(non-blocked version), and scalable decomposition of sub-domains. (sub-domains arranged in 3 dimensions.)
contains:
    * `Laplace` *Well-know jacobi iteration for testing*
    * `Lid_driven_cavity` <u> Lid-driven cavity flow. </u> 
    * `Thermal_flow`   Thermal flow.
    * `Particle_flow`   Particle flow.
    * `Multi_Phase_flow`    Multiphase flow.

* `OpenACC` OpenACC accelerated codes for single GPU.

* `MPI + OpenACC` Multi-GPU solver.

