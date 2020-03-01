import LinearAlgebra.det
import Statistics.var
import Statistics.std
using Printf
using StatsFuns
#using SortingAlgorithms
#using Distributions
using StatsBase
#using SharedArrays

################################################################################
# Gibbs algorithms
################################################################################

"""
Gibbs sampler using Kalman filter for the OU SDEMEM model.
"""
function gibbs_exact(R::Int,
                     y::Vector,
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
                     startval_log_σ_ϵ::Real,
                     startval_η::Vector,
                     prior_parameters_η::Array,
                     prior_parameters_σ_ϵ::Array)

    println("Starting Gibbs sampler using Kalman.")

    # dim for problem
    M = length(y)

    # initilization
    chain_ϕ = zeros(M, 3, R)
    chain_log_σ_ϵ = zeros(1, R)
    chain_η = zeros(6, R)
    ϕ_star = zeros(size(startval_ϕ))
    accept_vec = zeros(3,R)
    accept = false

    loglik_vec = zeros(R)

    α_prob_accept_σ_ϵ = 0.
    α_prob_accept_ϕ = zeros(M)

    # first iteration
    chain_ϕ[:,:,1] = startval_ϕ
    chain_log_σ_ϵ[1] = startval_log_σ_ϵ
    print_interval = 1000

    # set star values for η
    chain_η[:,1] = startval_η

    #=
    for j = 1:3

        # get paramaters
        μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

        chain_η[j,1] = μ_0_j # set start value to mean
        chain_η[j+3,1] = (α_j - 1)/β_j

    end
    =#

    accept_vec[:,1] = [M;1;3] # we accapt the start values

    # compute likelihood value
    loglik_old_y = kalman(y, exp(startval_log_σ_ϵ), startval_ϕ, dt)

    loglik_vec[1] = sum(loglik_old_y)

    # main loop
    for r = 2:R

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,:] = chain_ϕ[m,:,r-1] + cholesky(exp(log_λ_i_ϕ[m])*Σ_i_ϕ[m]).U*randn(3)
        end

        loglik_star_y = kalman(y, exp(chain_log_σ_ϵ[r-1]), ϕ_star, dt)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M

            jacobian_old =  sum(chain_ϕ[m,:,r-1])
            jacobian_star = sum(ϕ_star[m,:])

            log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] + jacobian_old)
            α_prob_accept_ϕ[m] = min(1, exp(log_α))

            # correct if we have  NaNs
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
        log_σ_ϵ_star = chain_log_σ_ϵ[r-1] + sqrt(exp(log_λ_i_σ_ϵ)*Σ_i_σ_ϵ)*randn()

        jacobian_old =  chain_log_σ_ϵ[r-1]
        jacobian_star = log_σ_ϵ_star

        loglik_star_y = kalman(y, exp(log_σ_ϵ_star), chain_ϕ[:,:,r], dt)

        prior_old_σ_ϵ = normlogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_log_σ_ϵ[r-1])
        prior_star_σ_ϵ = normlogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],log_σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_σ_ϵ + jacobian_star) - (sum(loglik_old_y) + prior_star_σ_ϵ + jacobian_old)
        α_prob_accept_σ_ϵ = min(1, exp(log_α))

        # correct if we have  NaNs
        if any(isnan,loglik_star_y) == true || log_σ_ϵ_star == NaN
            log_α = log(0)
            α_prob_accept_σ_ϵ = min(1, exp(log_α))
        end

        accept = log(rand()) < log_α

        if accept
            chain_log_σ_ϵ[r] = log_σ_ϵ_star
            accept_vec[2,r] = accept_vec[2,r] + 1
            loglik_old_y = loglik_star_y
        else
            chain_log_σ_ϵ[r] = chain_log_σ_ϵ[r-1]
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

        loglik_vec[r] = sum(loglik_old_y)

        print_progress(r,print_interval,accept_vec,M)

        log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ = adaptive_update_σ_ϵ(r,update_interval,start_update,log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ,γ_σ_ϵ_0,α_power,chain_log_σ_ϵ,α_prob_accept_σ_ϵ)
        adaptive_update_ϕ!(r,update_interval,start_update,log_λ_i_ϕ,μ_i_ϕ,Σ_i_ϕ,γ_ϕ_0,α_power,M,chain_ϕ,α_prob_accept_ϕ)



    end

    println("Ending Gibbs sampler using Kalman.")

    return chain_ϕ, chain_log_σ_ϵ, chain_η, accept_vec, loglik_vec

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
                     startval_log_σ_ϵ::Real,
                     startval_η::Vector,
                     prior_parameters_η::Array,
                     prior_parameters_σ_ϵ::Array,
                     nbr_particles::Vector,
                     T_vec::Vector,
                     ρ::Real,
                     particel_filter::Function)

    println("Starting Gibbs sampler using CPMMH.")

    # dim for problem
    M = length(y)

    # initilization
    chain_ϕ = zeros(M, 3, R)
    chain_log_σ_ϵ = zeros(1, R)
    chain_η = zeros(6, R)
    ϕ_star = zeros(size(startval_ϕ))
    accept_vec = zeros(3,R)
    accept = false
    loglik_vec = zeros(R)

    α_prob_accept_σ_ϵ = 0.
    α_prob_accept_ϕ = zeros(M)

    # first iteration
    chain_ϕ[:,:,1] = startval_ϕ
    chain_log_σ_ϵ[1] = startval_log_σ_ϵ
    print_interval = 1000

    if ρ > 0
        run_sort = true
    else
        run_sort = false
    end


    # set star values for η
    chain_η[:,1] = startval_η

    #=
    # sample star values for η
    for j = 1:3

        # get paramaters
        μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

        chain_η[j,1] = μ_0_j # set start value to mean
        chain_η[j+3,1] = (α_j - 1)/β_j # set start value to mode

    end
    =#

    accept_vec[:,1] = [M;1;3] # we accapt the start values

    # set random numbers
    #u_prop_old = randn(nbr_particles,N+1,M)
    #u_resample_old = randn(N,2,M)

    u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
    u_resample_old = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
        append!(u_resample_old, [randn(T_vec[i],2)])
    end

    # estimate likelihood valie
    loglik_old_y = particel_filter(y, exp(startval_log_σ_ϵ), startval_ϕ, dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)

    loglik_vec[1] = sum(loglik_old_y)

    #nbr_nan  = 0
    # main loop
    for r = 2:R

        # first gibbs stage, update \phi
        for m = 1:M
            ϕ_star[m,:] = chain_ϕ[m,:,r-1] + cholesky(exp(log_λ_i_ϕ[m])*Σ_i_ϕ[m]).U*randn(3)
        end

        u_prop_new = [randn(nbr_particles[1], T_vec[1]+1)]
        u_resample_new = [randn(T_vec[1], 2)]

        for i in 2:M_subjects;
            append!(u_prop_new, [randn(nbr_particles[i],T_vec[i]+1)])
            append!(u_resample_new, [randn(T_vec[i],2)])
        end


        u_prop_star = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new # TODO check this!
        u_resample_star = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new

        loglik_star_y = particel_filter(y, exp(chain_log_σ_ϵ[r-1]), ϕ_star, dt, u_prop_star, u_resample_star,  nbr_particles, run_sort)

        loglik_old_ϕ = calc_loglik_ϕ(chain_ϕ[:,:,r-1], chain_η[:,r-1])
        loglik_star_ϕ = calc_loglik_ϕ(ϕ_star, chain_η[:,r-1])

        for m = 1:M


            jacobian_old =  sum(chain_ϕ[m,:,r-1])
            jacobian_star = sum(ϕ_star[m,:])

            log_α = (loglik_star_y[m] + loglik_star_ϕ[m] + jacobian_star) - (loglik_old_y[m] + loglik_old_ϕ[m] + jacobian_old)
            α_prob_accept_ϕ[m] = min(1, exp(log_α))

            # correct if we have  NaNs
            if loglik_star_y[m] == NaN
                log_α = log(0)
                α_prob_accept_ϕ[m] = min(1, exp(log_α))
                nbr_nan = nbr_nan +1
            end


            accept = log(rand()) < log_α

            if accept
                chain_ϕ[m,:,r] = ϕ_star[m,:]
                accept_vec[1,r] = accept_vec[1,r] +1
                loglik_old_y[m] = loglik_star_y[m]
                u_prop_old[m] = u_prop_star[m]
                u_resample_old[m] = u_resample_star[m]
            else
                chain_ϕ[m,:,r] = chain_ϕ[m,:,r-1]
            end

        end

        # secound gibbs stage, update log_σ_ϵ
        log_σ_ϵ_star = chain_log_σ_ϵ[r-1] + sqrt(exp(log_λ_i_σ_ϵ)*Σ_i_σ_ϵ)*randn()

        jacobian_old =  chain_log_σ_ϵ[r-1]
        jacobian_star = log_σ_ϵ_star

        loglik_star_y = particel_filter(y, exp(log_σ_ϵ_star), chain_ϕ[:,:,r], dt, u_prop_old, u_resample_old,  nbr_particles, run_sort)


        prior_old_log_σ_ϵ = normlogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],chain_log_σ_ϵ[r-1])
        prior_star_log_σ_ϵ = normlogpdf(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2],log_σ_ϵ_star)

        log_α = (sum(loglik_star_y) + prior_star_log_σ_ϵ + jacobian_star) - (sum(loglik_old_y) + prior_star_log_σ_ϵ + jacobian_old)
        α_prob_accept_σ_ϵ = min(1, exp(log_α))

        # correct if we have  NaNs
        if any(isnan,loglik_star_y) == true
            #nbr_nan = nbr_nan +1
            #println(r)
            #println(log_σ_ϵ_star)
            #println(log_λ_i_σ_ϵ)
            #println(Σ_i_σ_ϵ)
            #println(chain_ϕ[:,:,r])
            #println(loglik_star_y)
            #println(particel_filter(y, exp(log_σ_ϵ_star), chain_ϕ[:,:,r], dt, u_prop_old, u_resample_old,  nbr_particles, run_sort))
            log_α = log(0)
            α_prob_accept_σ_ϵ = min(1, exp(log_α))
        end

        accept = log(rand()) < log_α


        if accept
            chain_log_σ_ϵ[r] = log_σ_ϵ_star
            accept_vec[2,r] = accept_vec[2,r] + 1
            loglik_old_y = loglik_star_y
        else
            chain_log_σ_ϵ[r] = chain_log_σ_ϵ[r-1]
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

        loglik_vec[r] = sum(loglik_old_y)

        print_progress(r,print_interval,accept_vec,M)

        log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ = adaptive_update_σ_ϵ(r,update_interval,start_update,log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ,γ_σ_ϵ_0,α_power,chain_log_σ_ϵ,α_prob_accept_σ_ϵ)
        adaptive_update_ϕ!(r,update_interval,start_update,log_λ_i_ϕ,μ_i_ϕ,Σ_i_ϕ,γ_ϕ_0,α_power,M,chain_ϕ,α_prob_accept_ϕ)

    end

    #println(nbr_nan)
    println("Ending Gibbs sampler using CPMMH.")

    return chain_ϕ, chain_log_σ_ϵ, chain_η, accept_vec, loglik_vec

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
    chain_σ_ϵ::Array,α_prob_accept_σ_ϵ::Real)

    if mod(r-1,update_interval) == 0 && r-1 >= start_update

        γ = γ_σ_ϵ_0/r^α_power

        log_λ_i_σ_ϵ = log_λ_i_σ_ϵ + γ*(α_prob_accept_σ_ϵ - α_star_σ_ϵ)
        μ_i_σ_ϵ = μ_i_σ_ϵ + γ*(chain_σ_ϵ[r] - μ_i_σ_ϵ)
        Σ_i_σ_ϵ = Σ_i_σ_ϵ + γ*((chain_σ_ϵ[r]-μ_i_σ_ϵ)^2 - Σ_i_σ_ϵ)
        # perhaps we need the scalar version of Σ_i_ϕ[m] = (Σ_i_ϕ_tmp + Σ_i_ϕ_tmp')/2 here # add the transpose and diviade by 2 for numerical stability (to make sure that Σ_i_ϕ[m] is symmertic)


    end

    return log_λ_i_σ_ϵ,μ_i_σ_ϵ,Σ_i_σ_ϵ
end


# Generlized AM adaptaton for the random effects ϕ,
# algorithm 4 in https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf)
function adaptive_update_ϕ!(r::Int,update_interval::Int,start_update::Int,
    log_λ_i_ϕ::Vector,μ_i_ϕ::Array,Σ_i_ϕ::Array,γ_ϕ_0::Real,α_power::Real,
    M::Int,chain_ϕ::Array,α_prob_accept_ϕ::Vector)

    if mod(r-1,update_interval) == 0 && r-1 >= start_update

        γ = γ_ϕ_0/r^α_power

        for m in 1:M
            Σ_i_ϕ_tmp = Σ_i_ϕ[m]
            log_λ_i_ϕ[m] =  log_λ_i_ϕ[m]  + γ*(α_prob_accept_ϕ[m] - α_star_ϕ)
            μ_i_ϕ[m,:] = μ_i_ϕ[m,:] + γ*(chain_ϕ[m,:,r] - μ_i_ϕ[m,:])
            Σ_i_ϕ_tmp = Σ_i_ϕ_tmp + γ*((chain_ϕ[m,:,r]-μ_i_ϕ[m,:])*(chain_ϕ[m,:,r]-μ_i_ϕ[m,:])' - Σ_i_ϕ_tmp)
            Σ_i_ϕ[m] = (Σ_i_ϕ_tmp + Σ_i_ϕ_tmp')/2 # add the transpose and diviade by 2 for numerical stability (to make sure that Σ_i_ϕ[m] is symmertic)
        end

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




# parallel kalman/pf filter code (without diagnostics)

################################################################################
# Kalman/pf
################################################################################

"""
Parallel Kalman filter
"""
function kalman(y_full::Array, σ_ϵ::Real, ϕ::Array, τ_s::Real, diagnostics::Bool=false)

    M = length(y_full)

    #loglik_est = zeros(M)

    if diagnostics; path = Vector{Float64}[]; end #[zeros(length(100))]; end

    #loglik_est = SharedArray{Float64}(M)
    loglik_est = zeros(M)

    #@sync @distributed for m in 1:M

    Threads.@threads for m in 1:M

        if diagnostics
            loglik_est[m],x = kalman_filter(y_full[m].mV, σ_ϵ, ϕ[m,:], τ_s, diagnostics)
            append!(path,[x])
        else
            loglik_est[m] = kalman_filter(y_full[m].mV, σ_ϵ, ϕ[m,:], τ_s, diagnostics)
        end

    end

    #diagnostics == true ? (return Array(loglik_est), path) : (return Array(loglik_est))
    diagnostics == true ? (return loglik_est, path) : (return loglik_est)


end

# Kalman filter. Estimate the likelihood for the OU process (LTI) system. The code follows the notation in "Grey-box pharmacokinetic/pharmacodynamic modelling of euglycaemic clamp study".
function kalman_filter(y::Vector, σ_ϵ::Real, ϕ::Vector, τ_s::Real, diagnostics::Bool=false)

    T = length(y)

    if diagnostics; x_hat_kalman = zeros(T); end

    #println(loglik_est[m])
    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3]) # this is actually the correct parameterization for the Kalman filter

    B = ν
    A = -λ
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

        #loglik_est = loglik_est - 0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)
        loglik_est = loglik_est + log_normalpdf(C*x_k, sqrt(R_k), y[k])

        x_k = x_k + K*ϵ_k
        P_k = P_k - K*R_k*K

        if diagnostics; x_hat_kalman[k] = x_k; end


    end

    diagnostics == true ? (return loglik_est, x_hat_kalman) : (return loglik_est)


end




"""
Parallel cpf
"""
function cpf(y_full::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Vector, u_resample::Vector, N_vec::Vector, run_sort::Bool)

    M = length(y_full)


    loglik_est = zeros(M)

    Threads.@threads for m in 1:M

        loglik_est[m] = particlefilter(y_full[m].mV, σ_ϵ, ϕ[m,:], dt, u_prop[m], u_resample[m], N_vec[m], run_sort, false)

    end

    return loglik_est

end


function cpf(y_full::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Vector, u_resample::Vector, N_vec::Vector, run_sort::Bool,diagnostics::Bool)

    M = length(y_full)


    if diagnostics;
        path = Matrix{Float64}[];
        nbr_resamlping = zeros(M)
    end

    loglik_est = zeros(M)

    for m in 1:M

        loglik_est[m],x,nbr_resamlping[m] = particlefilter(y_full[m].mV, σ_ϵ, ϕ[m,:], dt, u_prop[m], u_resample[m], N_vec[m], run_sort, diagnostics)
        push!(path, x)

    end

    diagnostics == true ? (return loglik_est, path, nbr_resamlping) : (return loglik_est)

end



# Estimates the loglikelihood using the correlated pf.
function particlefilter(y::Vector, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Real, run_sort::Bool, diagnostics::Bool=false)

    # pre-allocation

    loglik = 0.

    T = length(y)
    u_resample_calc = zeros(T)

    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3]) # this is actually the correct parameterization for the Kalman filter

    θ_1 = λ
    θ_2 = 1/λ*ν
    θ_3 = σ

    # convert standard normal to standard uniforms
    u_resample_temp = u_resample
    u_prop_temp = u_prop[:,2:end]
    u_prop_temp_init = u_prop[:,1]

    for i in 1:T
      u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
    end

    # pre-allocation
    x = zeros(N,T) # particels
    w = zeros(N,T) # weigts

    nbr_resamlping = T

    # set start values
    xinit = zeros(N) + std(y)/10*u_prop_temp_init #(std(y[m,:]) + σ_ϵ)*randn(N)

    for t in 1:T

        if t == 1 # first iteration

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
            w[:,t-1] = 1/N*ones(N)

            # propagate particels
            x[:,t] = stateprop(x_resample, θ_1, θ_2, θ_3, dt, u_prop_temp[:,t]) #r*x_resample.*exp.(-x_resample .+ e[:,t])

        end

        # calc weigths and update loglik
        loglik = loglik + calc_weigths(w,x[:,t],y[t],σ_ϵ,t)

    end

    diagnostics == true ? (return loglik, x,nbr_resamlping) : (return loglik)

end




function bridge_pf(y_full::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Vector, u_resample::Vector, N_vec::Vector, run_sort::Bool)

    M = length(y_full)

    loglik_est = zeros(M)


    Threads.@threads for m in 1:M
        loglik_est[m] = diffusionbridgeparticlefilter(y_full[m].mV, σ_ϵ, ϕ[m,:], dt, u_prop[m], u_resample[m], N_vec[m], run_sort,false)
    end

    return loglik_est

end



function bridge_pf(y_full::Array, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Vector, u_resample::Vector, N_vec::Vector, run_sort::Bool, diagnostics::Bool)

    M = length(y_full)

    if diagnostics;
        path = Matrix{Float64}[];
        nbr_resamlping = zeros(M)
    end

    loglik_est = zeros(M)

    for m in 1:M

        if diagnostics
            loglik_est[m],x,nbr_resamlping[m] = diffusionbridgeparticlefilter(y_full[m].mV, σ_ϵ, ϕ[m,:], dt, u_prop[m], u_resample[m], N_vec[m], run_sort,diagnostics)
            push!(path, x)
        else
            loglik_est[m] = diffusionbridgeparticlefilter(y_full[m].mV, σ_ϵ, ϕ[m,:], dt, u_prop[m], u_resample[m], N_vec[m], run_sort,diagnostics)
        end

    end

    diagnostics == true ? (return loglik_est, path, nbr_resamlping) : (return loglik_est)

end

# Estimates the loglikelihood using the correlated pf.
function diffusionbridgeparticlefilter(y::Vector, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Real, run_sort::Bool,diagnostics::Bool=false)

    # pre-allocation

    loglik = 0.

    T = length(y)
    u_resample_calc = zeros(T)

    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3]) # this is actually the correct parameterization for the Kalman filter


    # convert standard normal to standard uniforms
    u_resample_temp = u_resample
    u_prop_temp = u_prop[:,2:end]
    u_prop_temp_init = u_prop[:,1]

    for i in 1:T
      u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
    end

    # pre-allocation
    x = zeros(N,T) # particels
    w = zeros(N,T) # weigts

    nbr_resamlping = T

    # set start values
    xinit = zeros(N) + std(y)/10*u_prop_temp_init #(std(y[m,:]) + σ_ϵ)*randn(N)

    for t in 1:T

        if t == 1 # first iteration

          #x[:,1] = xinit

          #loglik = loglik + calc_weigths(w,x[:,t],y[t],σ_ϵ,t)

          # propagate particels
          x_resample = xinit
          x[:,t], p_hat_mean, p_hat_std = sample_partilces(N, y[t], x_resample, u_prop_temp[:,t], dt, λ, ν, σ, σ_ϵ)

        else

            if run_sort == true
                sorted_idx = sortperm(x[:,t-1]; alg = QuickSort) # ensure that we use QuickSort
                w[:,t-1] = w[sorted_idx,t-1]
                x[:,t-1] = x[sorted_idx,t-1]
            end

            # resample particels using systematic resampling
            ind = sysresample2(w[:,t-1],N,u_resample_calc[t])
            x_resample = x[ind,t-1]
            w[:,t-1] = 1/N*ones(N)

            # propagate particels
            x[:,t], p_hat_mean, p_hat_std = sample_partilces(N, y[t], x_resample, u_prop_temp[:,t], dt, λ, ν, σ, σ_ϵ)

        end

        loglik = loglik + calc_apf_weigths(N, w, x[:,t], x_resample, p_hat_mean, p_hat_std, y[t],σ_ϵ,t, dt, λ, ν, σ)

    end

    diagnostics == true ? (return loglik, x,nbr_resamlping) : (return loglik)

end

################################################################################
## Testing of particel filters
################################################################################

# Kalman filter. Estimate the likelihood for the OU process (LTI) system. The code follows the notation in "Grey-box pharmacokinetic/pharmacodynamic modelling of euglycaemic clamp study".
function kalman_filter_test(y::Vector, σ_ϵ::Real, ϕ::Vector, τ_s::Real, diagnostics::Bool=false)

    T = length(y)

    if diagnostics; x_hat_kalman = zeros(T); end

    #println(loglik_est[m])
    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3]) # this is actually the correct parameterization for the Kalman filter

    B = ν
    A = -λ
    C = 1
    S = σ_ϵ^2

    # start values for the kalman filter
    P_start = 0 #var(y)
    x_hat_start = 5 # TODO: This is not correct we can fix this by sampling the start values form N(0, σ_ϵ)

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

        #loglik_est = loglik_est - 0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)
        loglik_est = loglik_est + log_normalpdf(C*x_k, sqrt(R_k), y[k])

        x_k = x_k + K*ϵ_k
        P_k = P_k - K*R_k*K

        if diagnostics; x_hat_kalman[k] = x_k; end


    end

    diagnostics == true ? (return loglik_est, x_hat_kalman) : (return loglik_est)


end


# Estimates the loglikelihood using the correlated pf.
function particlefilter_test(y::Vector, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Real, run_sort::Bool, diagnostics::Bool=false)

    # pre-allocation

    loglik = 0.

    T = length(y)
    u_resample_calc = zeros(T)

    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3])

    θ_1 = λ
    θ_2 = 1/λ*ν
    θ_3 = σ

    # convert standard normal to standard uniforms
    u_resample_temp = u_resample
    u_prop_temp = u_prop[:,2:end]
    u_prop_temp_init = u_prop[:,1]

    for i in 1:T
      u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
    end

    # pre-allocation
    x = zeros(N,T) # particels
    w = zeros(N,T) # weigts

    nbr_resamlping = T

    # set start values
    xinit = 5*ones(N)

    for t in 1:T

        if t == 1 # first iteration

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
            w[:,t-1] = 1/N*ones(N)

            # propagate particels
            x[:,t] = stateprop(x_resample, θ_1, θ_2, θ_3, dt, u_prop_temp[:,t]) #r*x_resample.*exp.(-x_resample .+ e[:,t])

        end

        # calc weigths and update loglik
        loglik = loglik + calc_weigths(w,x[:,t],y[t],σ_ϵ,t)

    end

    diagnostics == true ? (return loglik, x,nbr_resamlping) : (return loglik)

end

# The bridge filter
function diffusionbridgeparticlefilter_test(y::Vector, σ_ϵ::Real, ϕ::Array, dt::Real, u_prop::Array, u_resample::Array, N::Real, run_sort::Bool,diagnostics::Bool=false)

    # pre-allocation

    loglik = 0.

    T = length(y)
    u_resample_calc = zeros(T)

    λ = exp(ϕ[1])
    ν = exp(ϕ[2])
    σ = exp(ϕ[3]) # this is actually the correct parameterization for the Kalman filter


    # convert standard normal to standard uniforms
    u_resample_temp = u_resample
    u_prop_temp = u_prop[:,2:end]
    u_prop_temp_init = u_prop[:,1]

    for i in 1:T
      u_resample_calc[i] = exp(-(u_resample_temp[i,1]^2+u_resample_temp[i,2]^2)/2)
    end

    # pre-allocation
    x = zeros(N,T) # particels
    w = zeros(N,T) # weigts

    nbr_resamlping = T

    # set start values
    xinit = 5*ones(N)

    for t in 1:T

        if t == 1 # first iteration

          #x[:,1] = xinit

          #loglik = loglik + calc_weigths(w,x[:,t],y[t],σ_ϵ,t)

          # propagate particels
          x_resample = xinit
          x[:,t], p_hat_mean, p_hat_std = sample_partilces(N, y[t], x_resample, u_prop_temp[:,t], dt, λ, ν, σ, σ_ϵ)


        else

            if run_sort == true
                sorted_idx = sortperm(x[:,t-1]; alg = QuickSort) # ensure that we use QuickSort
                w[:,t-1] = w[sorted_idx,t-1]
                x[:,t-1] = x[sorted_idx,t-1]
            end

            # resample particels using systematic resampling
            ind = sysresample2(w[:,t-1],N,u_resample_calc[t])
            x_resample = x[ind,t-1]
            w[:,t-1] = 1/N*ones(N)

            # propagate particels
            x[:,t], p_hat_mean, p_hat_std = sample_partilces(N, y[t], x_resample, u_prop_temp[:,t], dt, λ, ν, σ, σ_ϵ)

        end

        loglik = loglik + calc_apf_weigths(N, w, x[:,t], x_resample, p_hat_mean, p_hat_std, y[t],σ_ϵ,t, dt, λ, ν, σ)

    end

    diagnostics == true ? (return loglik, x,nbr_resamlping) : (return loglik)

end

function calc_apf_weigths(N::Int, w::Array, x_new::Vector,x_old::Vector, p_hat_mean::Vector, p_hat_std::Real, y::Real, σ_ϵ::Real, t::Int, dt::Real, λ::Real, ν::Real, σ::Real)


    # compute weigths
    x_transtion_std = sqrt((σ^2/(2*λ))*(1-exp(-2*λ*dt)))

    logw = zeros(N)
    w_temp = zeros(N)

    # here we do not need the t > 1 condition since we t > 1 now
    @inbounds for i in 1:N
        x_transtion_mean = 1/λ*ν + (x_old[i] - 1/λ*ν)*exp(-λ*dt)
        logw[i] = log_normalpdf(x_new[i], σ_ϵ, y) + log_normalpdf(x_transtion_mean, x_transtion_std, x_new[i]) - log_normalpdf(p_hat_mean[i], p_hat_std, x_new[i])
    end

    # find largets wegith
    constant = maximum(logw)

    # subtract largets weigth
    @inbounds for i in 1:N; w_temp[i] = exp(logw[i] - constant); end

    # calc sum of weigths
    w_sum = sum(w_temp)

    # normalize weigths
    for i in 1:N; w_temp[i] = w_temp[i]/w_sum; end
    w[:,t] = w_temp

    return constant + log(w_sum) - log(N)

end

function sample_partilces(N::Int, y::Real, x::Vector, u_prop_temp::Vector, dt::Real, λ::Real, ν::Real, σ::Real, σ_ϵ::Real)

    x_new = zeros(N)
    p_hat_mean = zeros(N)
    α_0 = zeros(N)

    @inbounds for i in 1:N; α_0[i] = x[i]*exp(-λ*dt)+ν/λ*(1-exp(-λ*dt)) ; end
    β_0 = σ^2/(2*λ)*(1-exp(-2*λ*dt))

    @inbounds for i in 1:N; p_hat_mean[i] = α_0[i] + β_0*(β_0 + σ_ϵ^2)^(-1)*(y-α_0[i]); end
    p_hat_std = sqrt(β_0*(1-(β_0 + σ_ϵ^2)^(-1)*β_0))

    @inbounds for i in 1:N; x_new[i] =  p_hat_mean[i] + p_hat_std*u_prop_temp[i]; end

    return x_new, p_hat_mean, p_hat_std

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

    # update loglik
    #loglik =  loglik + constant + log(w_sum) - log(N)

    # normalize weigths
    for i in 1:N; w_temp[i] = w_temp[i]/w_sum; end

    w[:,t] = w_temp

    # return loglik
    return constant + log(w_sum) - log(N)


end



function log_normalpdf(μ::Real, σ::Real, x::Real)

  R = σ^2
  ϵ = x-μ

  return -0.5*(log(2*pi) + log(det(R)) + ϵ*inv(R)*ϵ)

  #return -0.5*(log(det(R)) + ϵ*inv(R)*ϵ)

end

#normlogpdf(1, 1, 1)
#log_normalpdf(1,1,1)

#=
function log_normalpdf(μ::Real, σ::Real, x::Real)

    # normlogpdf(μ::Real, σ::Real, x::Number)
    return normlogpdf(μ, σ, x)

end
=#

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
