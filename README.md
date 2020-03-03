# Code for the paper: "Efficient inference for stochastic differential mixed-effects models using correlated particles pseudo-marginal algorithms"

Link to paper: https://arxiv.org/abs/1907.09851


N.B.: The results in the pre-print v2 at arXiv are computed using the version of the code at tag preprint_v2.

### File structure

/src/SDEMEM OU process: source code for the Ornstein-Uhlenbeck SDEMEM

/src/SDEMEM tumor growth: source code for tumor growth SDEMEM

/src/SDEMEM OU neuron data: source code for the neuronal data example 

/analyses: notebooks and scripts used to analyse the results for the Ornstein-Uhlenbeck SDEMEM adn the neuronal data example

/lunarc: run scripts for the LUNARC cluster (http://www.lunarc.lu.se/) (used for the Ornstein-Uhlenbeck SDEMEM and the neuronal data example)

/data/SDEMEM OU: Simulation results for the Ornstein-Uhlenbeck SDEMEM (this folder is empty but the simulation results can be generated from the code)

/data/SDEMEM OU neuron data: Simulation results for the neuronal data example (this folder is empty but the simulation results can be generated from the code)

### Software

The algorithms in the paper are implemented in `Julia` (Ornstein-Uhlenbeck SDEMEM, and neuronal data example), and in `R` (for the tumor growth SDEMEM).

### Data

The data used for for the  Ornstein-Uhlenbeck SDEMEM and the tumor growth model can be generated from the code. To acess the neuronal data, please contact the authers.
