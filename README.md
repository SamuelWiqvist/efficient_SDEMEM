# Code for the paper: "Efficient inference for stochastic differential mixed-effects models using correlated particles pseudo-marginal algorithms"

**privet in-progress repo**

### File structure

/src/SDEMEM OU process: source code for the Ornstein-Uhlenbeck SDEMEM

/src/SDEMEM tumor growth: source code for Tumor growth SDEMEM

/analyses: notebooks and scripts used to analyse the results for the Ornstein-Uhlenbeck SDEMEM

/lunarc: run scripts for the LUNARC cluster (http://www.lunarc.lu.se/) (used for the Ornstein-Uhlenbeck SDEMEM)

/data/SDEMEM OU: Simulation results Ornstein-Uhlenbeck SDEMEM (this folder is empty but the simulation results can be generated from the code)

### Software 

The algorithms in the paper are implemented in `Julia 1.0.0` (for the Ornstein-Uhlenbeck SDEMEM), and in `R xyz` (for the Tumor growth SDEMEM).

### Data 

The data used for two experiments can be generated from the code. 


