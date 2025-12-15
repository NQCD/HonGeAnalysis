# HonGeAnalysis

[![Build Status](https://github.com/NQCD/HonGeAnalysis/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/NQCD/HonGeAnalysis/actions/workflows/CI.yml)



The repository is developed by [Xuexun (or Hokseon in Cantonese)](https://louhokseson.github.io) and aims to reproduce the simulations and plots showed in the manuscript:
> [A Haldane-Anderson Model Hamiltonian for Hyperthermal Hydrogen Scattering from a Semiconductor Surface](https://about:blank) by Xuexun Lu, Nils Hertl, Sara Oregioni, Riley J Preston, Samuel Louis Rudge, Michael Thoss, Rocco Martinazzo and Reinhard J. Maurer.

<details>
<summary><strong>Schematic of H-Atom Scattering from Ge(111)</strong></summary>


<p align="center">
  <img src="README_fig/schematics/HonGeschematic.png" width="500" style="margin-right:20px;" />
</p><p align="center">
  <em>Figure 1. 3-D Schematic of an H atom scattering from Ge(111) surface.</em>
</p>

However, this project applies 1-D analytical model to represent the nonadiabatic effects during the scattering process. **Figure 1. is only for illustrative purpose.**
</details>

#### ONLY FOR Hokseon

<details>
<summary><strong>Hokseon needs to read this</strong></summary>

Hope the future me would thank the Hokseon who is wrapping up this repository. This repo should be ready in some machines that Hokseon can use by the time he comes back to it. i.e. Archer2 and Taskfarm.
- If you need to configure a new machine, please follow the step you wrote below. Because with the given UUID of NQCDynamics and NQCModels, julia package git installer would automatically download v0.15.0. Make sure that you do
```julia
(HonGeAnalysis) pkg> add NQCDynamics@0.13.4
```
in your configure script for that machine.Same applies to `NQCModels.jl` which has to be either v0.8.20 or v0.8.19.
- Make sure `data` folder has folder `sims`. Within `sims` should have `Ehrenfest`, `IESH` and `Individual-Large`. The first two are used for storing the raw `.h5` output from NQCD simulation. Each of the `.h5` should be named after the simulation parameters ([parameters_IESH.jl](/scripts/simulation_IESH/parameters_IESH.jl) for reference). The given run scripts (Ehrenfest/IESH) would skip taht simulation if the conjugate output `.h5` exist in folder `Ehrenfest` or `IESH`.
- When you need to do repeating simulations (mostly IESH for energy loss/ sticking probability convergence), make sure you turn the `is_dividual_large_saving = true` in the parameters_IESH.jl for simualtions. For example, 500 trajectories for a `.h5` for 1000 times.
- When you have generated massive amount of .h5 file under folder `Individual-Large/your_simulation_parameter...`, you need to process/extract the useful properties by 
1. Run [traj2kineticloss.jl](/scripts/data_engineering/traj2kineticloss.jl) and [traj2nstick.jl](/scripts/data_engineering/traj2nstick.jl) (order of executing does not matter since they are independent). These two generate folder `scattering_counting` and `scattered_kinetic_loss` containing `.csv` with the desired properities (including confidence errors) from the NQCD simulated `.h5` data.  The `.csv` is easy for storage and rsyncing between machines.
- Rsyncing the whole folder from HPCs to storage machine to process the [traj2kineticloss.jl](/scripts/data_engineering/traj2kineticloss.jl) and [traj2nstick.jl](/scripts/data_engineering/traj2nstick.jl). `rsync -avn` is a dry run. testing whether the stuff can be send to desitnation. The actually syncing need `rsync -av` 

**Whole folder rsyncing**
```bash
rsync -avn your-path-to-the-folder destination
```
**Files inside folder rsyncing**
```bash
rsync -avn your-path-to-the-folder/ destination
```

- Rsyncing the processed `.csv` to a local laptop:
**Exclude .h5 files**
```bash
rsync -avn --exclude='*.h5' source destination
```
</details>




## Step-by-step instructions from the author

### Activate the Julia environment
<details>
<summary><strong>Steps</strong></summary>

After cloning this repository to your machine, you need to the followings
1. Turn on the Julia REPL where you terminal need to be at the root of this repository
   ```julia
   (@v1.10) pkg> activate .
   ```
   Add `HokseonRegistry` to the julia environment so that you could access some personal packages Hokseon has developed for filing and plotting purposes.
   ```julia
   (HonGeAnalysis) pkg> registry add https://github.com/Louhokseson/HokseonRegistry
   ```
   Also, the NQCDRegistry is needed to access the NQCD family packages.
   ```julia
   (HonGeAnalysis) pkg> registry add https://github.com/NQCD/NQCRegistry
   ```
2. If you managed to activate the project environment, you should see the prompt
   ```julia
   (HonGeAnalysis) pkg> 
   ```
3. Resolve the dependencies / equivalent to rewrite the [Manifest.toml](./Manifest.toml) based on the [Project.toml](./Project.toml)
   ```julia
   (HonGeAnalysis) pkg> resolve
   ```
4. Instantate the dependencies / equivalent to download the packages listed on [Project.toml](./Project.toml)
   ```julia
   (HonGeAnalysis) pkg> instantiate
   ```
5. Last step is to double check you have the correct version of the `NQCDynamics.jl` (v0.13.4) and its companion package, `NQCModels.jl` (v0.8.19) in the julia enviroment.
   ```julia
   (HonGeAnalysis) pkg> status NQCDynamics
   ```
   If you see the version is not correct, you can use the following command to update the package.
   ```julia
   (HonGeAnalysis) pkg> add NQCDynamics@0.13.4
   ```
   Exact same procedure applies to `NQCModels.jl` package.
</details>

For those who are confident with Julia and skip the folded instructions, you are highly recommended to check your installed version of `NQCDynamics.jl` and `NQCModels.jl` packages are v0.13.4 and v0.8.19 respectively. If you are not sure, please follow the folded content in **Steps**.

#### Data availibility
While you are configurating your Julia environment with [HonGeAnalysis](https://github.com/NQCD/HonGeAnalysis), you are recommended to download the conjugate data of this repository in the meantime. Please unzip the downloaded compressed file under directory `HonGeAnalysis/`. It contains the following data:

##### Figure data
In the directory [figure_data](./figue_data/), you can find the data used to generate the figures in the manuscript. Those data are stored in `.txt` files, which are easy to be read by users.
##### NQCD Simulated data
You can find the simulated data from the `data/sims` directory.

1. You can find the raw data files with extension of `.h5` which contain the raw ouput from the [NQCDynamics.jl](https://github.com/NQCD/NQCDynamics.jl) simulation from the folder Ehrenfest and IESH. The name of the files indicate the simulation parameters, e.g.,
```
centre=0_couplings_rescale=1.95_discretisation=GapGaussLegendre_dt=0.01_gap=0.49_impuritymodel=Hokseon_incident_energy=0.3_is_Wigner_initial=false_mass=1.01_method=EhrenfestNA_nstates=150_temperature=300.0_tmax=1000_trajectories=1_width=50.h5
```
contains the simulation of 500 trajctories of the Ehrenfest dynamics with desired output variables along the dynamics.
2. The Individual_Large directory contains folders with name illustrating the simulation parameters, e.g.,
```
centre=0_couplings_rescale=1.95_decoherence=EDC_discretisation=GapGaussLegendre_dt=0.05_gap=0.49_impuritymodel=Hokseon_incident_energy=0.25_is_Wigner=false_mass=1.01_method=AdiabaticIESH_nstates=150_temperature=300.0_tmax=1001_trajectories=500_width=50
```
Inside, it contains the processed data from a large set of simulations and stored in the `.csv` files. Those `.csv` files are named with the distinct job id of each julia simulation. 


##### Ab-initio calculated data
The ab-initio data is stored in directory `data/ab-initio_cals`. The specs of the DFT calculations is illustrated in the repo conjugate paper.


##### Experimental data
The H atom scattering on Ge(111) surface data is published and available in 
> 
[Krüger, K., Wang, Y., Tödter, S. et al. Hydrogen atom collisions with a semiconductor efficiently promote electrons to the conduction band. Nat. Chem. 15, 326–331 (2023).](https://doi.org/10.1038/s41557-022-01085-x)




<details>
<summary><strong>Acknowledgments</strong></summary>

**Fundings**

- EPSRC Doctoral Training Partnership
- MSCA Postdoctoral Fellowship 
- Erasmus+ Traineeship Mobility 
- Alexander von Humboldt Research Fellowship
- DFG Grant
- UKRI Future Leaders Fellowship
- UKRI Frontier Research Grant

**Computing Resources**
- Scientific Computing Research Technology Platform (SCRTP) in University of Warwick
- Archer2 UK National Supercomputing Service
- Sulis HPC Midlands+ Computing Centre 
- bwHPC (Baden-Württemberg High Performance Computing)

**Hosting Institutions**

<p align="center">
  <a href="https://warwick.ac.uk" target="_blank">
    <img src="https://warwick.ac.uk/static_war/render/id7/images/wordmark.svg.136055278947" width="200" style="margin:10px;" alt="University of Warwick" />
  </a>
  <a href="https://www.uni-freiburg.de" target="_blank">
    <img src="https://www.bwhpc.de/img/Logos_Unis/UFR-Wortmarke-Grundform_Blau_CMYK.jpg" width="200" style="margin:10px;" alt="University of Freiburg" />
  </a>
  <a href="https://www.univie.ac.at" target="_blank">
    <img src="https://www.univie.ac.at/_assets/12d6b27de3745b4c3b247e6ea4f9bbcc/Images/Logos/univie.svg" width="200" style="margin:10px;" alt="University of Vienna" />
  </a>
  <a href="https://www.unimi.it" target="_blank">
    <img src="https://4euplus.eu/cuni_new_web/dist/images/4eu/logo_detail_milano_2x.png?v=1.1.1" width="200" style="margin:10px;" alt="University of Milan" />
  </a>
</p>



</details>



