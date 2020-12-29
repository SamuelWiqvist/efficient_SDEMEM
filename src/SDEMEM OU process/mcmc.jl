using LinearAlgebra
using Statistics
using Printf
using StatsFuns
#using SortingAlgorithms
#using Distributions
using StatsBase


################################################################################
# Gibbs algorithms
################################################################################

"""
Gibbs sampler using Kalman filter for the OU SDEMEM model.
"""
function gibbs_exact(R::Int,
                     y::Array,
                     dt::Real,
                     Σ_i_σ_ϵ::Real,
                     Σ_i_ϕ::Array,
                     γ_ϕ_0::Real,
                     γ_σ_ϵ_0::Real,
                     μ_i_ϕ::Array,
                     μ_i_σ_ϵ::Real,
                     α_star_ϕ::Real,
                     α_star_σ_ϵ::Real,
                     log_λ_i_ϕ::Vector,
                     log_λ_i_σ_ϵ::Real,
                     update_interval::Real,
                     start_update::Real,
                     α_power::Real,
                     startval_ϕ::Array,
                     startval_σ_ϵ::Real,
                     prior_parameters_η::Array,
                     prior_parameters_σ_ϵ::Array,
					 start_from_true_param::Bool=false)

    println("Starting Gibbs sampler using Kalman.")

    # dim for problem
    M,N = size(y)

    # initilization
    chain_ϕ = zeros(M, 3, R)
    chain_σ_ϵ = zeros(1, R)
    chain_η = zeros(6, R)
    ϕ_star = zeros(size(startval_ϕ))
    accept_vec = zeros(3,R)
    accept = false

    # first iteration
    chain_ϕ[:,:,1] = startval_ϕ
    chain_σ_ϵ[1] = startval_σ_ϵ
    print_interval = 1000

	if start_from_true_param
		chain_η[:,1] = η
	else
		# sample star values for η
		for j = 1:3

		    # get paramaters
		    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

		    chain_η[j,1] = μ_0_j # set start value to mean
		    chain_η[j+3,1] = (α_j - 1)/β_j

		end

	end 
    accept_vec[:,1] = [M;1;3] # we accapt the start values

    α_prob_accept_σ_ϵ = 0.
    α_prob_accept_ϕ = zeros(M)

    # compute likelihood value
    loglik_old_y = kalman(y, startval_σ_ϵ, startval_ϕ, dt)

    # main loop
    for r = 2:R

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,:] = chain_ϕ[m,:,r-1] + cholesky(exp(log_λ_i_ϕ[m])*Σ_i_ϕ[m]).U*randn(3)
        end

        loglik_star_y = kalman(y, chain_σ_ϵ[r-1], ϕ_star, dt)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M

            # TODO remove jacobian!
            #jacobian_old =  sum(chain_ϕ[m,:,r-1])
            #jacobian_star = sum(ϕ_star[m,:])

            #log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)
            log_α = (loglik_star_y[m] + loglik_star_ϕ[m]) - (loglik_old_y[m] + loglik_old_ϕ[m])

            α_prob_accept_ϕ[m] = min(1, exp(log_α))


            if any(isnan,ϕ_star[m,:]) == true || loglik_star_y[m] == NaN
                log_α = log(0)
                α_prob_accept_ϕ[m] = min(1, exp(log_α))
            end

            accept = log(rand()) < log_α


            if accept
                chain_ϕ[m,:,r] = ϕ_star[m,:]
                accept_vec[1,r] = accept_vec[1,r] +1
                loglik_old_y[m] = loglik_star_y[m]
            else
                chain_ϕ[m,:,r] = chain_ϕ[m,:,r-1]
            end

        end

        # secound gibbs stage, update σ_ϵ
        σ_ϵ_star = chain_σ_ϵ[r-1] + sqrt(exp(log_λ_i_σ_ϵ)*Σ_i_σ_ϵ)*randn()

        loglik_star_y = kalman(y, σ_ϵ_star, chain_ϕ[:,:,r], dt)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_old_σ_ϵ)
        α_prob_accept_σ_ϵ = min(1, exp(log_α))

        # correct if we have  NaNs
        if any(isnan,loglik_star_y) == true || σ_ϵ_star == NaN
            log_α = log(0)
            α_prob_accept_σ_ϵ = min(1, exp(log_α))
        end


        accept = log(rand()) < log_α

        if accept
            chain_σ_ϵ[r] = σ_ϵ_star
            accept_vec[2,r] = accept_vec[2,r] + 1
            loglik_old_y = loglik_star_y
        else
            chain_σ_ϵ[r] = chain_σ_ϵ[r-1]
        end

        # third stage, update η
        for j = 1:3

            # get paramaters
            μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

            ϕ_j = chain_ϕ[:,j,r-1]
            ϕ_j_mean = 1/M*sum(ϕ_j)

            α = α_j + div(M,2)
            β = β_j + 0.5*sum((ϕ_j .- ϕ_j_mean).^2) + (M*M_0_j)/(2*(M+M_0_j))*(ϕ_j_mean-μ_0_j)^2

            τ_j = rand_gamma(α)
            τ_j = τ_j/β

            μ = (M*τ_j*ϕ_j_mean)/(M*τ_j + M_0_j*τ_j) + (M_0_j*τ_j*μ_0_j)/(M*τ_j + M_0_j*τ_j)
            τ = M*τ_j + M_0_j*τ_j

            μ_j = μ + sqrt(1/τ)*randn()

            chain_η[j,r] = μ_j
            chain_η[j+3,r] = τ_j

            accept_vec[3,r] = accept_vec[3,r] +1

        end

        print_progress(r,print_interval,accept_vec,M)

        log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ = adaptive_update_σ_ϵ(r,update_interval,start_update,log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ,γ_σ_ϵ_0,α_power,chain_σ_ϵ,α_prob_accept_σ_ϵ,α_star_σ_ϵ)
        adaptive_update_ϕ!(r,update_interval,start_update,log_λ_i_ϕ,μ_i_ϕ,Σ_i_ϕ,γ_ϕ_0,α_power,M,chain_ϕ,α_prob_accept_ϕ,α_star_ϕ)

    end

    println("Ending Gibbs sampler using Kalman.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec


end



"""
Gibbs sampler using CPMMH for the OU SDEMEM model.

The naive version of the Gibbs scheme (Gibbs scheme #1) where we do not integrate out the random numbers.
"""
function naive_gibbs_cpmmh(R::Int,
                           y::Array,
                           dt::Real,
                           Σ_i_σ_ϵ::Real,
                           Σ_i_ϕ::Array,
                           γ_ϕ_0::Real,
                           γ_σ_ϵ_0::Real,
                           μ_i_ϕ::Array,
                           μ_i_σ_ϵ::Real,
                           α_star_ϕ::Real,
                           α_star_σ_ϵ::Real,
                           log_λ_i_ϕ::Vector,
                           log_λ_i_σ_ϵ::Real,
                           update_interval::Real,
                           start_update::Real,
                           α_power::Real,
                           startval_ϕ::Array,
                           startval_σ_ϵ::Real,
                           prior_parameters_η::Array,
                           prior_parameters_σ_ϵ::Array,
                           nbr_particles::Int,
                           ρ::Real)

    println("Starting Gibbs sampler using CPMMH.")

    # dim for problem
    M,N = size(y)

    # initilization
    chain_ϕ = zeros(M, 3, R)
    chain_σ_ϵ = zeros(1, R)
    chain_η = zeros(6, R)
    ϕ_star = zeros(size(startval_ϕ))
    accept_vec = zeros(3,R)
    accept = false

    # first iteration
    chain_ϕ[:,:,1] = startval_ϕ
    chain_σ_ϵ[1] = startval_σ_ϵ

    if ρ > 0
        run_sort = true
    else
        run_sort = false
    end


    # sample star values for η
    for j = 1:3

        # get paramaters
        μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

        chain_η[j,1] = μ_0_j # set start value to mean
        chain_η[j+3,1] = (α_j - 1)/β_j # set start value to mode

    end

    accept_vec[:,1] = [M;1;3] # we accapt the start values

    α_prob_accept_σ_ϵ = 0.
    α_prob_accept_ϕ = zeros(M)

    print_interval = 100
    # TODO check if the random numbers are correct, now, we have ONE set of random numbers
    # but they are updated for each block

    # set random numbers
    u_prop_old = randn(nbr_particles,N+1,M)
    u_resample_old = randn(N,2,M)

    # estimate likelihood valie
    loglik_old_y = cpf(y, startval_σ_ϵ, startval_ϕ, dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

    # main loop
    for r = 2:R

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,:] = chain_ϕ[m,:,r-1] + cholesky(exp(log_λ_i_ϕ[m])*Σ_i_ϕ[m]).U*randn(3)
        end


        # update random numbers for block 1
        u_prop_new_block_1 = randn(nbr_particles,N+1,M)
        u_resample_new_block_1 = randn(N,2,M_subjects)

        u_prop_star_block_1 = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new_block_1 # TODO check this!
        u_resample_star_block_1 = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new_block_1

        # calc log-lik
        loglik_star_y = cpf(y, chain_σ_ϵ[r-1], ϕ_star, dt, u_prop_star_block_1, u_resample_star_block_1,  nbr_particles, run_sort)
        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M
            
            # TODO remove jacobian!
            #jacobian_old =  sum(chain_ϕ[m,:,r-1])
            #jacobian_star = sum(ϕ_star[m,:])

            #log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)
            log_α = (loglik_star_y[m] + loglik_star_ϕ[m]) - (loglik_old_y[m] + loglik_old_ϕ[m])
            
            α_prob_accept_ϕ[m] = min(1, exp(log_α))

            if any(isnan,ϕ_star[m,:]) == true || loglik_star_y[m] == NaN
                log_α = log(0)
                α_prob_accept_ϕ[m] = min(1, exp(log_α))
            end

            accept = log(rand()) < log_α


            if accept
                chain_ϕ[m,:,r] = ϕ_star[m,:]
                accept_vec[1,r] = accept_vec[1,r] +1
                loglik_old_y[m] = loglik_star_y[m]
                u_prop_old[:,:,m] = u_prop_star_block_1[:,:,m] # update random numbers for block 1
                u_resample_old[:,:,m] = u_resample_star_block_1[:,:,m]
            else
                chain_ϕ[m,:,r] = chain_ϕ[m,:,r-1]
            end

        end

        # secound gibbs stage, update σ_ϵ
        σ_ϵ_star = chain_σ_ϵ[r-1] + sqrt(exp(log_λ_i_σ_ϵ)*Σ_i_σ_ϵ)*randn()

        # update random numbers for block 2
        u_prop_new_block_2 = randn(nbr_particles,N+1,M)
        u_resample_new_block_2 = randn(N,2,M_subjects)

        u_prop_star_block_2 = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new_block_2
        u_resample_star_block_2 = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new_block_2

        # calc loglik
        loglik_star_y = cpf(y, σ_ϵ_star, chain_ϕ[:,:,r], dt, u_prop_star_block_2, u_resample_star_block_2,  nbr_particles, run_sort)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_old_σ_ϵ)
        α_prob_accept_σ_ϵ = min(1, exp(log_α))

        # correct if we have  NaNs
        if any(isnan,loglik_star_y) == true || σ_ϵ_star == NaN
            log_α = log(0)
            α_prob_accept_σ_ϵ = min(1, exp(log_α))
        end

        accept = log(rand()) < log_α

        if accept
            chain_σ_ϵ[r] = σ_ϵ_star
            accept_vec[2,r] = accept_vec[2,r] + 1
            loglik_old_y = loglik_star_y
            u_prop_old[:,:,:] = u_prop_star_block_2[:,:,:] # update random numbers for block 2
            u_resample_old[:,:,:] = u_resample_star_block_2[:,:,:]
        else
            chain_σ_ϵ[r] = chain_σ_ϵ[r-1]
        end

        # third stage, update η
        for j = 1:3

            # get paramaters
            μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

            ϕ_j = chain_ϕ[:,j,r-1]
            ϕ_j_mean = 1/M*sum(ϕ_j)

            α = α_j + div(M,2)
            β = β_j + 0.5*sum((ϕ_j .- ϕ_j_mean).^2) + (M*M_0_j)/(2*(M+M_0_j))*(ϕ_j_mean-μ_0_j)^2

            τ_j = rand_gamma(α) # TODO check scaling for Gamma
            τ_j = τ_j/β

            μ = (M*τ_j*ϕ_j_mean)/(M*τ_j + M_0_j*τ_j) + (M_0_j*τ_j*μ_0_j)/(M*τ_j + M_0_j*τ_j)
            τ = M*τ_j + M_0_j*τ_j

            μ_j = μ + sqrt(1/τ)*randn()

            chain_η[j,r] = μ_j
            chain_η[j+3,r] = τ_j

            accept_vec[3,r] = accept_vec[3,r] +1

        end

        print_progress(r,print_interval,accept_vec,M)

        log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ = adaptive_update_σ_ϵ(r,update_interval,start_update,log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ,γ_σ_ϵ_0,α_power,chain_σ_ϵ,α_prob_accept_σ_ϵ,α_star_σ_ϵ)
        adaptive_update_ϕ!(r,update_interval,start_update,log_λ_i_ϕ,μ_i_ϕ,Σ_i_ϕ,γ_ϕ_0,α_power,M,chain_ϕ,α_prob_accept_ϕ, α_star_ϕ)



    end

    println("Ending Gibbs sampler using CPMMH.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec

end


"""
Gibbs sampler using CPMMH for the OU SDEMEM model.

Gibbs scheme #2 where we integrate out the random numbers.
"""
function gibbs_cpmmh(R::Int,
                     y::Array,
                     dt::Real,
                     Σ_i_σ_ϵ::Real,
                     Σ_i_ϕ::Array,
                     γ_ϕ_0::Real,
                     γ_σ_ϵ_0::Real,
                     μ_i_ϕ::Array,
                     μ_i_σ_ϵ::Real,
                     α_star_ϕ::Real,
                     α_star_σ_ϵ::Real,
                     log_λ_i_ϕ::Vector,
                     log_λ_i_σ_ϵ::Real,
                     update_interval::Real,
                     start_update::Real,
                     α_power::Real,
                     startval_ϕ::Array,
                     startval_σ_ϵ::Real,
                     startval_η::Array,
                     prior_parameters_η::Array,
                     prior_parameters_σ_ϵ::Array,
                     nbr_particles::Int,
                     ρ::Real)

    println("Starting Gibbs sampler using CPMMH.")

    # dim for problem
    M,N = size(y)

    # initilization
    chain_ϕ = zeros(M, 3, R)
    chain_σ_ϵ = zeros(1, R)
    chain_η = zeros(6, R)
    ϕ_star = zeros(size(startval_ϕ))
    accept_vec = zeros(3,R)
    accept = false

    # first iteration
    chain_ϕ[:,:,1] = startval_ϕ
    chain_σ_ϵ[1] = startval_σ_ϵ

    if ρ > 0
        run_sort = true
    else
        run_sort = false
    end

    chain_η[:,1] = startval_η 

    accept_vec[:,1] = [M;1;3] # we accapt the start values

    α_prob_accept_σ_ϵ = 0.
    α_prob_accept_ϕ = zeros(M)

    print_interval = 100

    # set random numbers
    u_prop_old = randn(nbr_particles,N+1,M)
    u_resample_old = randn(N,2,M)

    # estimate likelihood valie
    loglik_old_y = cpf(y, startval_σ_ϵ, startval_ϕ, dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

    # main loop
    for r = 2:R

        for m = 1:M
            ϕ_star[m,:] = chain_ϕ[m,:,r-1] + cholesky(exp(log_λ_i_ϕ[m])*Σ_i_ϕ[m]).U*randn(3)
        end

        u_prop_new = randn(nbr_particles,N+1,M)
        u_resample_new = randn(N,2,M_subjects)

        u_prop_star = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new # TODO check this!
        u_resample_star = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new

        loglik_star_y = cpf(y, chain_σ_ϵ[r-1], ϕ_star, dt, u_prop_star, u_resample_star,  nbr_particles, run_sort)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M


            # TODO remove jacobian!
            # jacobian correct 
            #jacobian_old =  sum(chain_ϕ[m,:,r-1])
            #jacobian_star = sum(ϕ_star[m,:])

            # no prior here since we have a normal prior on the log scale 
            #log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)
            log_α = (loglik_star_y[m] + loglik_star_ϕ[m]) - (loglik_old_y[m] + loglik_old_ϕ[m])

            α_prob_accept_ϕ[m] = min(1, exp(log_α))

            if any(isnan,ϕ_star[m,:]) == true || loglik_star_y[m] == NaN
                log_α = log(0)
                α_prob_accept_ϕ[m] = min(1, exp(log_α))
            end

            accept = log(rand()) < log_α

            if accept
                chain_ϕ[m,:,r] = ϕ_star[m,:]
                accept_vec[1,r] = accept_vec[1,r] +1
                loglik_old_y[m] = loglik_star_y[m]
                u_prop_old[:,:,m] = u_prop_star[:,:,m]
                u_resample_old[:,:,m] = u_resample_star[:,:,m]
            else
                chain_ϕ[m,:,r] = chain_ϕ[m,:,r-1]
            end

        end

        # secound gibbs stage, update σ_ϵ
        σ_ϵ_star = chain_σ_ϵ[r-1] + sqrt(exp(log_λ_i_σ_ϵ)*Σ_i_σ_ϵ)*randn()

        loglik_star_y = cpf(y, σ_ϵ_star, chain_ϕ[:,:,r], dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        # TODO fix???
        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_old_σ_ϵ)
        α_prob_accept_σ_ϵ = min(1, exp(log_α))

        # correct if we have  NaNs
        if any(isnan,loglik_star_y) == true || σ_ϵ_star == NaN
            log_α = log(0)
            α_prob_accept_σ_ϵ = min(1, exp(log_α))
        end


        accept = log(rand()) < log_α

        if accept
            chain_σ_ϵ[r] = σ_ϵ_star
            accept_vec[2,r] = accept_vec[2,r] + 1
            loglik_old_y = loglik_star_y
        else
            chain_σ_ϵ[r] = chain_σ_ϵ[r-1]
        end

        # third stage, update η
        for j = 1:3

            # get paramaters
            μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

            ϕ_j = chain_ϕ[:,j,r-1]
            ϕ_j_mean = 1/M*sum(ϕ_j)

            α = α_j + div(M,2)
            β = β_j + 0.5*sum((ϕ_j .- ϕ_j_mean).^2) + (M*M_0_j)/(2*(M+M_0_j))*(ϕ_j_mean-μ_0_j)^2

            τ_j = rand_gamma(α) # TODO check scaling for Gamma
            τ_j = τ_j/β

            μ = (M*τ_j*ϕ_j_mean)/(M*τ_j + M_0_j*τ_j) + (M_0_j*τ_j*μ_0_j)/(M*τ_j + M_0_j*τ_j)
            τ = M*τ_j + M_0_j*τ_j

            μ_j = μ + sqrt(1/τ)*randn()

            chain_η[j,r] = μ_j
            chain_η[j+3,r] = τ_j

            accept_vec[3,r] = accept_vec[3,r] +1

        end


        print_progress(r,print_interval,accept_vec,M)

        log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ = adaptive_update_σ_ϵ(r,update_interval,start_update,log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ,γ_σ_ϵ_0,α_power,chain_σ_ϵ,α_prob_accept_σ_ϵ,α_star_σ_ϵ)
        adaptive_update_ϕ!(r,update_interval,start_update,log_λ_i_ϕ,μ_i_ϕ,Σ_i_ϕ,γ_ϕ_0,α_power,M,chain_ϕ,α_prob_accept_ϕ,α_star_ϕ)

    end

    println("Ending Gibbs sampler using CPMMH.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec

end


################################################################################
# help functions for Gibbs algorithms
################################################################################

# print progress to consol
function print_progress(r::Int,print_interval::Int,accept_vec::Array,M::Int)

    if mod(r-1,print_interval) == 0
        println("-------------")
        @printf "Percentage done:  %.2f %%\n" r/R*100
        @printf "Accaptance rate on %.0f to %.0f for ϕ: %.2f %%:\n" r-print_interval r sum(accept_vec[1,r-print_interval:r])/(M*print_interval)*100
        @printf "Accaptance rate on %.0f to %.0f for σ_ϵ: %.2f %%:\n" r-print_interval r sum(accept_vec[2,r-print_interval:r])/print_interval*100
    end

end

# Generlized AM adaptaton for the sigma_epsilon parameter,
# algorithm 4 in https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf)
function adaptive_update_σ_ϵ(r::Int,update_interval::Int,start_update::Int,
    log_λ_i_σ_ϵ::Real,μ_i_σ_ϵ::Real,Σ_i_σ_ϵ::Real,γ_σ_ϵ_0::Real,α_power::Real,
    chain_σ_ϵ::Array,α_prob_accept_σ_ϵ::Real,α_star_σ_ϵ::Real)

    if mod(r-1,update_interval) == 0 && r-1 >= start_update && r < 60000

        #γ = γ_σ_ϵ_0/(r-start_update)^α_power
        γ = γ_σ_ϵ_0/(r)^α_power

        log_λ_i_σ_ϵ = log_λ_i_σ_ϵ + γ*(α_prob_accept_σ_ϵ - α_star_σ_ϵ)
        μ_i_σ_ϵ = μ_i_σ_ϵ + γ*(chain_σ_ϵ[r] - μ_i_σ_ϵ)
        Σ_i_σ_ϵ = Σ_i_σ_ϵ + γ*((chain_σ_ϵ[r]-μ_i_σ_ϵ)^2 - Σ_i_σ_ϵ)

    end

    return log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ
end


# Generlized AM adaptaton for the random effects ϕ,
# algorithm 4 in https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf)
function adaptive_update_ϕ!(r::Int,update_interval::Int,start_update::Int,
    log_λ_i_ϕ::Vector,μ_i_ϕ::Array,Σ_i_ϕ::Array,γ_ϕ_0::Real,α_power::Real,
    M::Int,chain_ϕ::Array,α_prob_accept_ϕ::Vector,α_star_ϕ::Real)

    if mod(r-1,update_interval) == 0 && r-1 >= start_update && r < 60000

        #γ = γ_ϕ_0/(r-start_update)^α_power
        γ = γ_ϕ_0/(r)^α_power

        for m in 1:M
            Σ_i_ϕ_tmp = Σ_i_ϕ[m]
            log_λ_i_ϕ[m] =  log_λ_i_ϕ[m]  + γ*(α_prob_accept_ϕ[m] - α_star_ϕ)
            μ_i_ϕ[m,:] = μ_i_ϕ[m,:] + γ*(chain_ϕ[m,:,r] - μ_i_ϕ[m,:])
            Σ_i_ϕ_tmp = Σ_i_ϕ_tmp + γ*((chain_ϕ[m,:,r]-μ_i_ϕ[m,:])*(chain_ϕ[m,:,r]-μ_i_ϕ[m,:])' - Σ_i_ϕ_tmp)
            Σ_i_ϕ[m] = (Σ_i_ϕ_tmp + Σ_i_ϕ_tmp')/2 # add the transpose and diviade by 2 for numerical stability (to make sure that Σ_i_ϕ[m] is symmertic)
        end

    end

end

# calc loglik for ϕ
function calc_loglik_ϕ(ϕ, η)

    loglik = zeros(size(ϕ,1))

    for m = 1:size(ϕ,1)
        for i in 1:size(ϕ,2)
            loglik[m] = loglik[m] + normlogpdf(η[i], sqrt(1/η[i+3]), ϕ[m,i])
        end
    end

    return loglik

end

# Sample from Gamma(α,1), code adapted from https://www.cs.toronto.edu/~radford/csc2541.F04/gamma.html
function rand_gamma(α::Real)

    x = 0

    # TODO: We should have some check here such that alpha is a natrual number
    for i = 1:Int64(α)
        x = x + log(rand())
    end

    return -x

end


################################################################################
# Kalman/pf
################################################################################

# TODO rewrite the kalman and cpf in the same "style" as for the neuronal data model, to see if parallelization works


"""
Parallel Kalman filter
"""
function kalman(y_full::Array, σ_ϵ::Real, ϕ::Array, τ_s::Real, diagnostics::Bool=false)

    M = size(y_full,1)

    loglik_est = zeros(M)

    if diagnostics; path = Vector{Float64}[]; end #[zeros(length(100))]; end


    for m in 1:M

        if diagnostics
            loglik_est[m],x = kalman_filter(y_full[m,:], σ_ϵ, ϕ[m,:], τ_s,diagnostics)
            append!(path,[x])
        else
            loglik_est[m] = kalman_filter(y_full[m,:], σ_ϵ, ϕ[m,:], τ_s,diagnostics)
        end

    end

    diagnostics == true ? (return loglik_est, path) : (return loglik_est)


end

# Kalman filter. Estimate the likelihood for the OU process (LTI) system. The code follows the notation in "Grey-box pharmacokinetic/pharmacodynamic modelling of euglycaemic clamp study".
function kalman_filter(y::Vector, σ_ϵ::Real, ϕ::Vector, τ_s::Real, diagnostics::Bool=false)

    T = length(y)

    if diagnostics; x_hat_kalman = zeros(T); end

    #println(loglik_est[m])
    θ_1 = exp(ϕ[1])
    θ_2 = exp(ϕ[2])
    θ_3 = exp(ϕ[3])

    # set model
    B = θ_1*θ_2
    A = -θ_1
    σ = θ_3
    C = 1
    S = σ_ϵ^2


    # start values for the kalman filter
    P_start = var(y)
    x_hat_start = 0 # TODO: This is not correct we can fix this by sampling the start values form N(0, σ_ϵ)

    P_k = P_start
    x_k = x_hat_start

    loglik_est = 0.


    # main loop
    for k = 1:T

        x_k = exp(A*τ_s)*x_k + (1/A)*(exp(τ_s*A)-1)*B
        P_k = exp(A*τ_s)*P_k*exp(A*τ_s) + σ^2*(1/(2*A))*(exp(2*A*τ_s)-1)

        ϵ_k = y[k]-C*x_k
        R_k = C*P_k*C + S
        K = P_k*C*inv(R_k)
        x_k = x_k + K*ϵ_k
        P_k = P_k - K*R_k*K

        loglik_est = loglik_est -0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)

        if diagnostics; x_hat_kalman[k] = x_k; end


    end

    diagnostics == true ? (return loglik_est, x_hat_kalman) : (return loglik_est)


end



"""
Parallel cpf
"""
function cpf(y_full::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Int, run_sort::Bool,diagnostics::Bool=false)

    M = size(y_full,1)

    loglik_est = zeros(M)

    if diagnostics; path = Matrix{Float64}[]; end

    for m in 1:M

        if diagnostics
            loglik_est[m],x = particlefilter(y_full[m,:], σ_ϵ, ϕ[m,:], dt, u_prop[:, :, m], u_resample[:, :, m], N, run_sort, diagnostics)
            push!(path, x)
        else
            loglik_est[m] = particlefilter(y_full[m,:], σ_ϵ, ϕ[m,:], dt, u_prop[:, :, m], u_resample[:, :, m], N, run_sort, diagnostics)
        end

    end

    diagnostics == true ? (return loglik_est, path) : (return loglik_est)

end


# Estimates the loglikelihood using the correlated pf.
function particlefilter(y::Vector, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample_temp::Array, N::Real, run_sort::Bool, diagnostics::Bool=false)

    # pre-allocation

    loglik = 0.

    T = length(y)
    u_resample_calc = zeros(T)

    θ_1 = exp(ϕ[1]) # set model parameters
    θ_2 = exp(ϕ[2])
    θ_3 = exp(ϕ[3])

    # convert standard normal to standard uniforms
    #u_resample_temp = u_resample
    u_prop_temp = u_prop[:,2:end]
    #u_prop_temp_init = u_prop[:,1]

    for i in 1:T
      u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
    end

    # pre-allocation
    x = zeros(N,T) # particels
    w = zeros(N,T) # weigts

    # set start values
    xinit = zeros(N) + std(y)*u_prop[:,1] #(std(y[m,:]) + σ_ϵ)*randn(N)

    #==    
    if N > 20
        xinit = zeros(N) + std(y)*u_prop[:,1] #(std(y[m,:]) + σ_ϵ)*randn(N)
    elseif N <= 20 && N > 5
        xinit = zeros(N) + std(y)/10*u_prop[:,1] #(std(y[m,:]) + σ_ϵ)*randn(N)
    else
        xinit = zeros(N) + std(y)/20*u_prop[:,1] #(std(y[m,:]) + σ_ϵ)*randn(N)
    end
    ==#

    for t in 1:T

        if t == 1 # first iteration

          #x[:,1] = xinit

          # propagate particelsend
          #x[:,1] = stateprop(xinit, θ_1, θ_2, θ_3, dt)
          #for i = 1:N; w[i,1] = 1/N; end;
          #x[:,1] = stateprop(xinit, θ_1, θ_2, θ_3, dt)

          # propagate particels
          x[:,t] = stateprop(xinit, θ_1, θ_2, θ_3, dt, u_prop_temp[:,t]) #r*x_resample.*exp.(-x_resample .+ e[:,t])

        else

            # sort particles and weigths
            # alg=RadixSort
            if run_sort == true
                sorted_idx = sortperm(x[:,t-1]; alg = QuickSort) # ensure that we use QuickSort
                w[:,t-1] = w[sorted_idx,t-1]
                x[:,t-1] = x[sorted_idx,t-1]
            end

            # resample particels using systematic resampling
            ind = sysresample2(w[:,t-1],N,u_resample_calc[t])
            x_resample = x[ind,t-1]

            # propagate particels
            x[:,t] = stateprop(x_resample, θ_1, θ_2, θ_3, dt, u_prop_temp[:,t]) #r*x_resample.*exp.(-x_resample .+ e[:,t])

        end

        # calc weigths and update loglik
        loglik = loglik + calc_weigths(w,x[:,t],y[t],σ_ϵ,t)

    end

    diagnostics == true ? (return loglik, x) : (return loglik)

end


################################################################################
# help functions for kalman/pf
################################################################################

# sample exactly from the OU process
function stateprop(x_old::Vector, θ_1::Real, θ_2::Real, θ_3::Real, dt::Real, u::Array)

    x_prop = zeros(length(x_old))

    σ_cond = sqrt((θ_3^2/(2*θ_1))*(1-exp(-2*θ_1*dt)))

    for i in 1:length(x_old)
        μ_cond = θ_2 + (x_old[i] - θ_2)*exp(-θ_1*dt)
        x_prop[i] = μ_cond + σ_cond*u[i]
    end

    return x_prop

end



# calc weigths and update loglik
function calc_weigths(w::Array, x::Array, y::Real, σ_ϵ::Real, t::Int)


    N = length(x) # nbr particels
    logw = zeros(N)
    w_temp = zeros(N)

    # calc w
    for i in 1:N; logw[i] = log_normalpdf(x[i], σ_ϵ, y); end

    # find largets wegith
    constant = maximum(logw)

    # subtract largets weigth
    for i in 1:N; w_temp[i] = exp(logw[i] - constant); end

    # calc sum of weigths
    w_sum = sum(w_temp)

    # normalize weigths
    for i in 1:N; w_temp[i] = w_temp[i]/w_sum; end

    w[:,t] = w_temp

    # return loglik
    return constant + log(w_sum) - log(N)


end


# calc weigths and update loglik
function calc_weigths!(w::Array, loglik::Array, x::Array, y::Real, σ_ϵ::Real, m::Int, t::Int)


    N = length(x) # nbr particels
    logw = zeros(N)
    w_temp = zeros(N)

    # calc w
    @inbounds for i in 1:N; logw[i] = log_normalpdf(x[i], σ_ϵ, y); end

    # find largets wegith
    constant = maximum(logw)

    # subtract largets weigth
    @inbounds for i in 1:N; w_temp[i] = exp(logw[i] - constant); end

    # calc sum of weigths
    w_sum = sum(w_temp)

    # update loglik
    loglik[m] =  loglik[m] + constant + log(w_sum) - log(N)

    # normalize weigths
    @inbounds for i in 1:N; w_temp[i] = w_temp[i]/w_sum; end

    w[:,t] = w_temp


end

# Log-pdf for the normal distribution
function log_normalpdf(μ::Real, σ::Real, x::Real)

  R = σ^2
  ϵ = x-μ
  return -0.5*(log(det(R)) + ϵ*inv(R)*ϵ)
end


# Systematic resampling. Code adapted from Andrew Golightly.
function sysresample2(wts::Array,N::Int64,uni::Real)

  vec = zeros(Int64,N)
  wsum = sum(wts)
  k = 1
  u = uni/N
  wsumtarg = u
  wsumcurr=wts[k]/wsum
  delta = 1/N

  for i = 1:N
    while wsumcurr < wsumtarg
      k = k+1
      wsumcurr=wsumcurr+wts[k]/wsum
    end
    vec[i]=k
    wsumtarg=wsumtarg+delta
  end

  return vec

end
