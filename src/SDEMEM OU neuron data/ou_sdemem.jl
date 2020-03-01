# load packages
using CSV
using Random

include(pwd()*"/src/SDEMEM OU neuron data/load_neuron_data.jl")

"""
    set_up(M::Int=312, seed::Int=100)

Returns data and prior distributions.
"""
function set_up(M::Int=312, seed::Int=100)

    # define priors
    μ_0_1 = log(1/10); M_0_1 = 1; α_1 = 2; β_1 = 1
    μ_0_2 = log(1.5); M_0_2 = 1; α_2 = 2; β_2 = 1
    μ_0_3 = log(0.5); M_0_3 = 1; α_3 = 2; β_3 = 1

    prior_parameters_η = [μ_0_1 M_0_1 α_1 β_1;μ_0_2 M_0_2 α_2 β_2;μ_0_3 M_0_3 α_3 β_3]
    prior_parameters_σ_ϵ = [-1; 1/1]

    # load data
    ISI_all, nbr_ISE = load_neuron_data()
    y = ISI_all[1:M]

    # return: data and parameteres for priors
    return y,prior_parameters_η,prior_parameters_σ_ϵ

end
