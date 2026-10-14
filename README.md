# Static and Dynamic Aeroelastic Analysis of a Strut-Braced Wing Configuration

Projects developed within the Structural Dynamics and Aeroelasticity course (M.Sc. in Aeronautical Engineering at Politecnico di Milano, A.Y. 2025-2026).

## Project Overview
This repository contains two comprehensive assignments focused on the static and dynamic aeroelastic characterization, numerical modeling, and structural optimization of a cantilevered wing reinforced by an intermediate elastic strut, featuring a concentrated tip fuel tank in the dynamic phase. All numerical simulations, matrix formulations, and post-processing analyses were implemented in MATLAB using the Ritz-Galerkin assumed-modes method.

---

## Assignment 1: Static Aeroelasticity and Structural Optimization
This module addresses the static aeroelastic response and design optimization of the braced wing.

### Key Analyses and Features:
* **Problem Modeling & Ritz-Galerkin Discretization:** Formulation of coupled bending and torsional governing equations using polynomial trial functions satisfying clamped-free root and tip boundary conditions.
* **Principle of Virtual Work (PVW):** Derivation of internal and external virtual work incorporating wing strain energy, discrete strut elastic contributions, and aerodynamic loads.
* **Divergence Dynamic Pressure Analysis:** Evaluation of critical divergence speed ($q_D$) as a function of the strut chordwise attachment position ($x_B$), demonstrating complete immunity to divergence for forward strut locations.
* **Control Effectiveness & Reversal:** Mapping of aileron control effectiveness ($E_c$) across dynamic pressures and calculation of control reversal speed ($q_R$).
* **Deformed Shape Comparison:** Reconstruction and comparison of spanwise bending and torsional displacements between the unbraced (clean) and braced configurations.
* **Convergence Study:** Verification of numerical stability and convergence of the Ritz polynomial basis.
* **Design Optimization:** Parametric sweep and optimization of the strut inclination angle ($\gamma$) to achieve a 40% increase in divergence pressure while minimizing mass penalties.

---

## Assignment 2: Dynamic Aeroelasticity, Flutter, and Gust Response
This module extends the formulation to unsteady aerodynamics, structural dynamics, flutter boundaries, and gust loads in the presence of a tip tank.

### Key Analyses and Features:
* **Generalized Mass, Stiffness, and Damping Matrices:** Rigorous derivation incorporating distributed wing mass, structural stiffness, strut equivalent spring, and a concentrated tip fuel tank modeled as a point mass with polar inertia.
* **Free-Vibration Modal Convergence:** Iterative modal analysis to determine optimal truncation limits ($N_w = 4$ bending modes, $N_{\theta} = 8$ torsion modes) based on natural frequency relative error thresholds ($\epsilon < 10^{-3}$).
* **Frequency Response Function (FRF):** Computation of the aeroelastic Bode plot evaluating harmonic aileron deflection input ($\beta$) to tip leading-edge displacement ($w_P$).
* **Flutter Analysis:** Parametric state-space flutter speed ($V_F$) evaluation across varying strut chordwise attachment positions ($x_B/c$), identifying critical instability boundaries and mode coupling mechanisms.
* **Step Gust Response:** State-space transient integration of root bending moment $M_x(0,t)$ subjected to a vertical step gust at subcritical flight speeds ($75\% V_F$).
* **Stochastic White-Noise Gust Analysis:** Evaluation of steady-state root bending moment variance ($\sigma^2$) under continuous uniform vertical white-noise turbulence via the algebraic Lyapunov equation.

---

## Repository Structure
* `Assignment1graziani.m`: MATLAB script containing all static aeroelasticity, divergence, reversal, and optimization routines.
* `Assignment2graziani.m`: MATLAB script containing dynamic mass/stiffness matrix assembly, modal convergence, FRF Bode plots, flutter sweep, step gust, and Lyapunov stochastic analysis.
* `numReportGraziani.pdf` / `reportGraziani.pdf`: Detailed technical report for Assignment 1.
* `numericalreportGraziani.pdf` / `analiticreportGraziani.pdf`: Detailed technical report for Assignment 2.

## Author
* **Matteo Marwan Graziani** (Matricola: 307177)

*Institution:* Politecnico di Milano — School of Industrial and Information Engineering  
*Course Professor:* Prof. Giuseppe Quaranta
