# HonGeAnalysis

This code base is using the [Julia Language](https://julialang.org/).

The repository is developed by [Hokseon(Xuexun)](https://louhokseson.github.io) and aims to reproduce the simulations and plots showed in the manuscript:
> [A Haldane-Anderson Model Hamiltonian for Hyperthermal Hydrogen Scattering from a Semiconductor Surface](https://doi.org/10.1016/j.cpc.2023.108885) by Xuexun Lu, Nils Hertl and Reinhard Maurer.

This project is part of Xuexun's PhD proejct funded by the


## Step-by-step instructions from the author
#### Activate the Julia environment
After cloning this repository to your machine, you need to the followings
1. Turn on the Julia REPL where you terminal need to be at the root of this repository
   ```julia
   (@v1.10) pkg> activate .
   ```
   Add `HokseonRegistry` to the julia environment so that you could access some personal packages Hokseon has developed for filing and plotting purposes.
   ```julia
   (HonGeAnalysis) pkg> registry add https://github.com/Louhokseson/HokseonRegistry
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

