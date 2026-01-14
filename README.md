# Coordinated Chromosome Oscillation

Simulation and analysis code supporting the study of coordinated chromosome
oscillations in the metaphase spindle.

## Overview
This repository contains MATLAB code used to:

1. Simulate kinetochore-driven oscillations of individual sister chromosomes.
2. Model coordinated chromosome oscillations arising from mechanical coupling
   between chromosomes, implemented as stochastic inter-chromosomal springs.

The code accompanies the manuscript:

> **Coordinated chromosome motion emerges from mechanical coupling mediated by the physical spindle environment**

and was used to generate the simulation results and figures reported therein.

---

## Repository structure

### 1. Single sister-chromosome oscillation model
**`Model1_OneChromosome/`**

- `One_Oscillating_Chromosome_main.m`  
  → Figure S2C

- `One_Oscillating_Chromosome_noise.m`  
  → Figure 2B

- `One_Oscillating_Chromosome_noise_multi_iteration.m`  
  → Figure 2C, Figure 3D

---

### 2. Negative control: Independent chromosome oscillations
**`Model2_IndependentChromosomes/`**

These simulations model multiple chromosomes oscillating independently, serving
as a negative control for coordination analyses.

- `multiple_chromosomes_noise.m`
- `multiple_chromosomes_noise_iteration.m`

---

### 3. Coordinated multi-chromosome oscillation model
**`Model3_InterConnectedChromosomes/`**

These simulations include stochastic mechanical coupling between chromosome
centers-of-mass.

- `Multiple_connected_chromosome_noise.m`  
  → Figure 4C (top, middle), Figure S3C (top, middle)

- `Multiple_connected_chromosome_noise_iteration.m`  
  → Figure 4C (bottom), Figure S3C (bottom), Figure 6B

- `Multiple_connected_chromosome_noise_iteration_l0.m`  
  → Figure 6C

- `Multiple_connected_chromosome_noise_iteration_k0.m`  
  → Figure 6D

---

### 4. Parameter sweeps and coordination analysis

**Single-chromosome coordination sweeps**  
(located in `Model1_OneChromosome/`)
- `Model1_OneChromosome/determinate_Sweep_Activity_bifurcation.m`  
  → Figure S2B, Figure S2D (left)

- `Model1_OneChromosome/noise_Sweep_Activity_bifurcation.m`  
  → Figure 3C, Figure S2D (right)
  
**Multi-chromosome coordination sweeps**  
(located in `Model3_InterConnectedChromosomes/`)
- `Model3_InterConnectedChromosomes/Alpha_Kct_Sweep_MSV.m`  
  → Figure S3B

- `Model3_InterConnectedChromosomes/Alpha_Kct_Sweep_Pearson.m`  
  → Figure S3A

- `Model3_InterConnectedChromosomes/kon_koff_sweep_connectivity.m`  
  → Figure S4A, Figure S4C

- `Model3_InterConnectedChromosomes/Kon_koff_Sweep_1stCrosscorrelation.m`  
  → Figure S4A

- `Model3_InterConnectedChromosomes/l0_K_Sweep_Pearson.m`  
  → Figure 4B

- `Model3_InterConnectedChromosomes/konkoff.m`  
  → Figure S4B

- `Model3_InterConnectedChromosomes/spring_connect.m`  
  → Dependency function for stochastic inter-chromosomal coupling

---

## Requirements
- MATLAB **R2024a or later**
- Signal Processing Toolbox (required for peak/valley detection using `findpeaks`)

---

## Usage
1. Set the repository root directory as your MATLAB working directory.
2. Run the main simulation scripts within each model directory to reproduce the
   corresponding analyses and figures.
3. Analysis functions are called internally by the simulation scripts and do not
   need to be executed separately.

---

## Reproducibility notes
- Each simulation script documents parameter values used for the corresponding
  manuscript figure.
- Stochastic simulations are averaged over multiple iterations as described in
  script headers and in the manuscript.
- Output files (e.g. CSV, figures) are intentionally excluded from version
  control; only code is tracked.

---

## License
This repository is released under the MIT License.
