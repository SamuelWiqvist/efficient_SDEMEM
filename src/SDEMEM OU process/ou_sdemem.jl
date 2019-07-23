# load packages
using CSV
using Random

"""
    set_up(;μ_θ_1::Real = -0.7, τ_θ_1::Real = 1.5, μ_θ_2::Real = 1.5, τ_θ_2::Real = 1.2, μ_θ_3::Real = -0.9, τ_θ_3::Real = 1.5, σ_ϵ::Real = 0.3, dt::Real=0.1, M::Int=5, N::Int=100)

Returns ground-truth parameter values and simulated data from the OU SDEMEM model.
"""
function set_up(;μ_θ_1::Real = -0.7, τ_θ_1::Real = 4,
                 μ_θ_2::Real = 2.3, τ_θ_2::Real = 10,
                 μ_θ_3::Real = -0.9, τ_θ_3::Real = 4,
                 σ_ϵ::Real = 0.3, dt::Real=0.05, M::Int=5, N::Int=100,seed::Int)


    # sample random effects
    Random.seed!(seed)

    ϕ = zeros(M,3)

    for i = 1:M
        ϕ[i,1] = μ_θ_1 + sqrt(1/τ_θ_1)*randn()
        ϕ[i,2] = μ_θ_2 + sqrt(1/τ_θ_2)*randn()
        ϕ[i,3] = μ_θ_3 + sqrt(1/τ_θ_3)*randn()
    end

    # model paramters
    η = [μ_θ_1; μ_θ_2; μ_θ_3; τ_θ_1; τ_θ_2; τ_θ_3]

    # define priors
    μ_0_1 = 0
    M_0_1 = 1
    α_1 = 2
    β_1 = 1

    μ_0_2 = 1
    M_0_2 = 1
    α_2 = 2
    β_2 = 1/2

    μ_0_3 = 0
    M_0_3 = 1
    α_3 = 2
    β_3 = 1

    prior_parameters_η = [μ_0_1 M_0_1 α_1 β_1;μ_0_2 M_0_2 α_2 β_2;μ_0_3 M_0_3 α_3 β_3]

    prior_parameters_σ_ϵ = [1; 1/0.4] # Gamma with scale parameter α = 1 and rate parameter β = 1/0.4

    # generate data

    Random.seed!(seed)
    y,x,t_vec = generate_data(N, M, η, ϕ, σ_ϵ,dt)

    # return: data, model parameteres, and parameteres for priors
    return y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ

end

"""
    generate_data(N::Int, M::Int, η::Array, ϕ::Array, σ_ϵ::Real,dt::Real)

Generate data from the OU SDEMEM model.
"""
function generate_data(N::Int, M::Int, η::Array, ϕ::Array, σ_ϵ::Real,dt::Real)

    x_0 = zeros(M) # pre-allocate matriceis
    x = zeros(M,N)
    y = zeros(M,N)
    t_vec = zeros(N)

    for m = 1:M

        x[m,1] = x_0[m] # set start value
        y[m,1] = x[m,1] + σ_ϵ*randn()

        θ_1 = exp(ϕ[m,1]) # set parameters for subject m
        θ_2 = exp(ϕ[m,2])
        θ_3 = exp(ϕ[m,3])

        σ_cond = sqrt((θ_3^2/(2*θ_1))*(1-exp(-2*θ_1*dt))) # set latent process std for subject m

        # simulate process for subject m
        for t = 2:N
            t_vec[t] = t_vec[t-1] + dt
            μ_cond = θ_2 + (x[m,t-1]-θ_2)*exp(-θ_1*dt)
            x[m,t] = μ_cond + σ_cond*randn()
            y[m,t] = x[m,t] + σ_ϵ*randn()
        end
    end

    # return data
    return y,x,t_vec

end
