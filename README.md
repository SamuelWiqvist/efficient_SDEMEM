# Code for the paper: "Efficient inference for stochastic differential mixed-effects models using correlated particles pseudo-marginal algorithms"

Link to paper paper: https://arxiv.org/abs/1907.09851


N.B.: The results in the pre-print at arXiv are computed using the version of the code at tag preprint_v1.

### File structure

/src/SDEMEM OU process: source code for the Ornstein-Uhlenbeck SDEMEM

/src/SDEMEM tumor growth: source code for tumor growth SDEMEM

/analyses: notebooks and scripts used to analyse the results for the Ornstein-Uhlenbeck SDEMEM

/lunarc: run scripts for the LUNARC cluster (http://www.lunarc.lu.se/) (used for the Ornstein-Uhlenbeck SDEMEM)

/data/SDEMEM OU: Simulation results Ornstein-Uhlenbeck SDEMEM (this folder is empty but the simulation results can be generated from the code)

### Software

The algorithms in the paper are implemented in `Julia` (for the Ornstein-Uhlenbeck SDEMEM), and in `R` (for the tumor growth SDEMEM).

### Data

The data used for two experiments can be generated from the code.
