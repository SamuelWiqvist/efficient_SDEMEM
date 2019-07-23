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
    gibbs_exact(R::Int,
                y::Array,
                dt::Real,
                cov_ϕ::Array,
                cov_σ_ϵ::Real,
                startval_ϕ::Array,
                startval_σ_ϵ::Real,
                prior_parameters_η::Array,
                prior_parameters_σ_ϵ::Array)

Gibbs sampler using Kalman filter for the OU SDEMEM model.
"""
function gibbs_exact(R::Int,
                     y::Array,
                     dt::Real,
                     cov_ϕ::Array,
                     cov_σ_ϵ::Real,
                     startval_ϕ::Array,
                     startval_σ_ϵ::Real,
                     prior_parameters_η::Array,
                     prior_parameters_σ_ϵ::Array)

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

    # sample star values for η
    for j = 1:3

        # get paramaters
        μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

        chain_η[j,1] = μ_0_j # set start value to mean
        chain_η[j+3,1] = (α_j - 1)/β_j

    end

    accept_vec[:,1] = [M;1;3] # we accapt the start values

    # compute likelihood value
    loglik_old_y = kalman(y, startval_σ_ϵ, startval_ϕ, dt)

    # main loop
    for r = 2:R

        print_progress(r, R, 1000)

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,1] = chain_ϕ[m,1,r-1] + cov_ϕ[1]*randn()
            ϕ_star[m,2] = chain_ϕ[m,2,r-1] + cov_ϕ[2]*randn()
            ϕ_star[m,3] = chain_ϕ[m,3,r-1] + cov_ϕ[3]*randn()
        end

        loglik_star_y = kalman(y, chain_σ_ϵ[r-1], ϕ_star, dt)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M

            jacobian_old =  sum(chain_ϕ[m,:,r-1])
            jacobian_star = sum(ϕ_star[m,:])

            log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)

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
        σ_ϵ_star = chain_σ_ϵ[r-1] + cov_σ_ϵ*randn()

        loglik_star_y = kalman(y, σ_ϵ_star, chain_ϕ[:,:,r], dt)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_star_σ_ϵ)

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

    end

    println("Ending Gibbs sampler using Kalman.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec

end



"""
    naive_gibbs_cpmmh(R::Int,
                y::Array,
                dt::Real,
                cov_ϕ::Array,
                cov_σ_ϵ::Real,
                startval_ϕ::Array,
                startval_σ_ϵ::Real,
                prior_parameters_η::Array,
                prior_parameters_σ_ϵ::Array,
                nbr_particles::Int,
                ρ::Real)

Gibbs sampler using CPMMH for the OU SDEMEM model.

The naive version of the Gibbs scheme (Gibbs scheme #1) where we do not integrate out the random numbers.
"""
function naive_gibbs_cpmmh(R::Int,
                           y::Array,
                           dt::Real,
                           cov_ϕ::Array,
                           cov_σ_ϵ::Real,
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

    # TODO check if the random numbers are correct, now, we have ONE set of random numbers
    # but they are updated for each block

    # set random numbers
    u_prop_old = randn(nbr_particles,N+1,M)
    u_resample_old = randn(N,2,M)

    # estimate likelihood valie
    loglik_old_y = cpf(y, startval_σ_ϵ, startval_ϕ, dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

    # main loop
    for r = 2:R

        print_progress(r, R, 1000)

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,1] = chain_ϕ[m,1,r-1] + cov_ϕ[1]*randn()
            ϕ_star[m,2] = chain_ϕ[m,2,r-1] + cov_ϕ[2]*randn()
            ϕ_star[m,3] = chain_ϕ[m,3,r-1] + cov_ϕ[3]*randn()
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

            jacobian_old =  sum(chain_ϕ[m,:,r-1])
            jacobian_star = sum(ϕ_star[m,:])

            log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)

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
        σ_ϵ_star = chain_σ_ϵ[r-1] + cov_σ_ϵ*randn()

        # update random numbers for block 2
        u_prop_new_block_2 = randn(nbr_particles,N+1,M)
        u_resample_new_block_2 = randn(N,2,M_subjects)

        u_prop_star_block_2 = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new_block_2
        u_resample_star_block_2 = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new_block_2

        # calc loglik
        loglik_star_y = cpf(y, σ_ϵ_star, chain_ϕ[:,:,r], dt, u_prop_star_block_2, u_resample_star_block_2,  nbr_particles, run_sort)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_star_σ_ϵ)

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


    end

    println("Ending Gibbs sampler using CPMMH.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec

end


"""
    gibbs_cpmmh(R::Int,
                y::Array,
                dt::Real,
                cov_ϕ::Array,
                cov_σ_ϵ::Real,
                startval_ϕ::Array,
                startval_σ_ϵ::Real,
                prior_parameters_η::Array,
                prior_parameters_σ_ϵ::Array,
                nbr_particles::Int,
                ρ::Real)

Gibbs sampler using CPMMH for the OU SDEMEM model.

Gibbs scheme #2 where we integrate out the random numbers.
"""
function gibbs_cpmmh(R::Int,
                     y::Array,
                     dt::Real,
                     cov_ϕ::Array,
                     cov_σ_ϵ::Real,
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

    # set random numbers
    u_prop_old = randn(nbr_particles,N+1,M)
    u_resample_old = randn(N,2,M)

    # estimate likelihood valie
    loglik_old_y = cpf(y, startval_σ_ϵ, startval_ϕ, dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

    # main loop
    for r = 2:R

        print_progress(r, R, 1000)

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,1] = chain_ϕ[m,1,r-1] + cov_ϕ[1]*randn()
            ϕ_star[m,2] = chain_ϕ[m,2,r-1] + cov_ϕ[2]*randn()
            ϕ_star[m,3] = chain_ϕ[m,3,r-1] + cov_ϕ[3]*randn()
        end

        u_prop_new = randn(nbr_particles,N+1,M)
        u_resample_new = randn(N,2,M_subjects)

        u_prop_star = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new # TODO check this!
        u_resample_star = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new

        loglik_star_y = cpf(y, chain_σ_ϵ[r-1], ϕ_star, dt, u_prop_star, u_resample_star,  nbr_particles, run_sort)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M

            jacobian_old =  sum(chain_ϕ[m,:,r-1])
            jacobian_star = sum(ϕ_star[m,:])

            log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] - jacobian_old)

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
        σ_ϵ_star = chain_σ_ϵ[r-1] + cov_σ_ϵ*randn()

        loglik_star_y = cpf(y, σ_ϵ_star, chain_ϕ[:,:,r], dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

        prior_old_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_σ_ϵ[r-1])
        prior_star_σ_ϵ = gammalogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ) - (sum(loglik_old_y) + prior_star_σ_ϵ)

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


    end

    println("Ending Gibbs sampler using CPMMH.")

    return chain_ϕ, chain_σ_ϵ, chain_η, accept_vec

end


################################################################################
# help functions for Gibbs algorithms
################################################################################


function print_progress(i::Int, N::Int, print_interval::Int=10000)
    if mod(i,print_interval) == 0
        @printf("Percentage done:  %.2f %%\n", i/N*100)
    end
end


function calc_loglik_ϕ(ϕ, η)

    loglik = zeros(size(ϕ,1))

    for m = 1:size(ϕ,1)
        for i in 1:size(ϕ,2)
            loglik[m] = loglik[m] + normlogpdf(η[i], sqrt(1/η[i+3]), ϕ[m,i])
        end
    end

    return loglik

end


# sample from Gamma(α,1)
# code adapted from https://www.cs.toronto.edu/~radford/csc2541.F04/gamma.html
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

"""
    kalman(y::Array, σ_ϵ::Real, ϕ::Array, τ_s::Real, return_path_est::Bool=false, return_partial_loglik::Bool=false)


Estimate the likelihood for the OU process (LTI) system. The code follows the notation in "Grey-box pharmacokinetic/pharmacodynamic modelling of euglycaemic clamp study".
"""
function kalman(y::Array, σ_ϵ::Real, ϕ::Array, τ_s::Real, return_path_est::Bool=false, return_partial_loglik::Bool=false)

    M, N = size(y)

    x_hat_kalman = zeros(size(y))
    y_hat_kalman = zeros(size(y))
    residuals = zeros(size(y))

    loglik_est = zeros(M)

    if return_partial_loglik
        partial_loglik = zeros(M,N+1)
    end

    for m = 1:M

        #println(loglik_est[m])
        θ_1 = exp(ϕ[m,1])
        θ_2 = exp(ϕ[m,2])
        θ_3 = exp(ϕ[m,3])

        # set model
        B = θ_1*θ_2
        A = -θ_1
        σ = θ_3
        C = 1
        S = σ_ϵ^2

        # start values for the kalman filter
        P_start = var(y[m,:])
        x_hat_start = 0 # TODO: This is not correct we can fix this by sampling the start values form N(0, σ_ϵ)

        P_k = P_start
        x_k = x_hat_start

        # main loop
        for k = 1:N

            if k == 1
                # dont update
            else
                x_k = exp(A*τ_s)*x_k + (1/A)*(exp(τ_s*A)-1)*B
                P_k = exp(A*τ_s)*P_k*exp(A*τ_s) + σ^2*(1/(2*A))*(exp(2*A*τ_s)-1)
            end

            ϵ_k = y[m,k]-C*x_k
            R_k = C*P_k*C + S
            K = P_k*C*inv(R_k)
            x_k = x_k + K*ϵ_k
            P_k = P_k - K*R_k*K

            y_hat_kalman[m,k] = C*x_k
            x_hat_kalman[m,k] = x_k

            # calc loglik
            #loglik_est[m] = loglik_est[m] - log(det(R_k)) - ϵ_k*inv(R_k)*ϵ_k
            loglik_est[m] = loglik_est[m]-0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)

            if return_partial_loglik
                partial_loglik[m,k+1] = loglik_est[m] #-0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)
            end

            #println(loglik_est[m])

            residuals[m,k] = ϵ_k
        end


    end

    if return_path_est
        return loglik_est, y_hat_kalman, x_hat_kalman, residuals
    elseif return_partial_loglik
        return loglik_est, partial_loglik
    else
        return loglik_est
    end

end

"""
    cpf(y::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Int = 1000, return_path_est::Bool=false, return_partial_loglik::Bool=false)


Estimates the loglikelihood using the correlated pf.
"""
function cpf(y::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Int, run_sort::Bool, return_path_est::Bool=false, return_partial_loglik::Bool=false)

    # pre-allocation
    M, T = size(y)
    loglik = zeros(M)
    u_resample_calc = zeros(T)

    if return_path_est
        x_return = zeros(M,N,T)
        w_return = zeros(M,N,T)
        residuals = zeros(size(y))
    end

    if return_partial_loglik
        partial_loglik = zeros(M,T+1)
    end

    for m in 1:M

        θ_1 = exp(ϕ[m,1]) # set model parameters
        θ_2 = exp(ϕ[m,2])
        θ_3 = exp(ϕ[m,3])

        # convert standard normal to standard uniforms
        u_resample_temp = u_resample[:,:,m]
        u_prop_temp = u_prop[:,2:end,m]
        u_prop_temp_init = u_prop[:,1,m]

        for i in 1:T
          u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
        end

        # pre-allocation
        x = zeros(N,T) # particels
        w = zeros(N,T) # weigts

        # set start values
        xinit = zeros(N) + std(y[m,:])*u_prop_temp_init #(std(y[m,:]) + σ_ϵ)*randn(N)

        for t in 1:T

            if t == 1 # first iteration

              x[:,1] = xinit

              # propagate particelsend
              #x[:,1] = stateprop(xinit, θ_1, θ_2, θ_3, dt)
              #for i = 1:N; w[i,1] = 1/N; end;
              #x[:,1] = stateprop(xinit, θ_1, θ_2, θ_3, dt)

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
            if return_path_est
                residuals[m,t] = y[m,t] - mean(x[:,t])
            end

            calc_weigths!(w,loglik,x[:,t],y[m,t],σ_ϵ,m,t)

            if return_partial_loglik
                partial_loglik[m,t+1] = loglik[m]
            end


        end

        if return_path_est
            x_return[m,:,:] = x
            w_return[m,:,:] = w
        end
    end

    if return_partial_loglik && return_path_est
        println("test")
        return loglik, w_return, x_return, partial_loglik
    elseif return_path_est
        return loglik, w_return, x_return, residuals
    elseif return_partial_loglik
        return loglik, partial_loglik
    else
        return loglik
    end


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
