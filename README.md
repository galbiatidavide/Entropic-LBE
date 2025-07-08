# C++ Entropic-Lattice-Boltzmann Scheme

I developed a high-performance MPI-parallelized code to efficiently solve classical computational fluid dynamics (CFD) benchmarks, including the 2D lid-driven cavity problem. The implementation includes tailored boundary conditions relevant to microfluidic regimes, where the particulate nature of the flow becomes dominant. In such non-continuum settings, the classical Navier–Stokes equations break down, and conventional no-slip boundary conditions are no longer valid.

To address these challenges, I implemented the Entropic Lattice Boltzmann Method (ELBM) — an unconditionally stable and robust approach particularly well-suited for modeling rarefied and microscale flows. ELBM naturally captures non-equilibrium effects and ensures numerical stability even in regimes where traditional methods fail.

The project includes modular, scalable code designed for extensibility and high-resolution simulations, with current applications targeting flow configurations where the continuum hypothesis breaks down and kinetic effects play a crucial role.

To protect the originality and performance advantages of the full codebase, I have opted to publish only a simplified 2D version on GitHub. This public release is intentionally basic, designed primarily to illustrate the core algorithmic structure and key steps of the Entropic Lattice Boltzmann Method, and relies on traditional dynamic memory allocation rather than the optimized PETSc-based backend used in the complete implementation.

🔬 The project is actively evolving, with ongoing developments focused on extending the method to more rarefied regimes.

This project was developed in the context of the course *Methods and Models for Statistical Mechanics*.

---

✨🚀🧪

**Powered by Davide Galbiati**
