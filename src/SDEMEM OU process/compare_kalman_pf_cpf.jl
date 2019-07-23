# script to run inference for the OU SDEMEM model

#using Pkg
using PyPlot
using LinearAlgebra
using KernelDensity
using DataFrames
using KernelDensity
using Distributions

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40

nbr_particles_pf = 1000
nbr_particles_cpf = 1000

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=100)

PyPlot.figure()
PyPlot.plot(Array(x'))

PyPlot.figure()
PyPlot.plot(Array(y'))

@time loglik_kalman, partial_loglik_sum_kalman = kalman(y, σ_ϵ, ϕ, dt, false, true)

@time loglik_pf, w, x, partial_loglik_sum_pf = cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles_cpf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_pf, false, true, true)

@time loglik_corr_pf, w_corr_pf, x_corr_pf, partial_loglik_sum_corr_pf = cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles_cpf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_cpf, true, true, true)


loglik_kalman

loglik_pf

loglik_corr_pf


partial_loglik_sum_kalman

partial_loglik_sum_pf

partial_loglik_sum_corr_pf


idx = 25

PyPlot.figure()
PyPlot.plot(partial_loglik_sum_pf[idx,:], "r")
PyPlot.plot(partial_loglik_sum_kalman[idx,:], "b")
PyPlot.plot(partial_loglik_sum_corr_pf[idx,:], "--r")

partial_loglik_kalman = diff(partial_loglik_sum_kalman,dims=2)
partial_loglik_pf = diff(partial_loglik_sum_pf,dims=2)
partial_loglik_cpf = diff(partial_loglik_sum_corr_pf,dims=2)

PyPlot.figure()
PyPlot.plot(partial_loglik_pf[idx,:], "r")
PyPlot.plot(partial_loglik_kalman[idx,:], "b")
PyPlot.plot(partial_loglik_cpf[idx,:], "--r")


loglik_pf_m = zeros(M_subjects,100)
loglik_cpf_m = zeros(M_subjects,100)

for i in 1:100
    loglik_pf_m[:,i] = cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles_pf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_pf, false)
    loglik_cpf_m[:,i] = cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles_cpf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_cpf, true)
end

loglik_pf_m

loglik_cpf_m

loglik_pf_m = sum(loglik_pf_m, dims = 1)

loglik_cpf_m = sum(loglik_cpf_m, dims = 1)

std(loglik_pf_m)

std(loglik_cpf_m)

mean(loglik_pf_m,dims=2)

mean(loglik_cpf_m,dims=2)



loglik_pg_dist = kde(loglik_pf_m[:])
loglik_pg_dist_cpf = kde(loglik_cpf_m[:])

PyPlot.figure()
PyPlot.plot(loglik_pg_dist.x,loglik_pg_dist.density, "b")
PyPlot.plot(loglik_pg_dist_cpf.x,loglik_pg_dist_cpf.density, "--b")
PyPlot.plot((sum(loglik_kalman), sum(loglik_kalman)), (0, maximum(loglik_pg_dist.density)), "k")


# test for other parameter values
prior_dist_σ_ϵ = Gamma(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2])

prior_dist_τ_1 = Gamma(prior_parameters_η[1,3],1/prior_parameters_η[1,4])
prior_dist_τ_2 = Gamma(prior_parameters_η[2,3],1/prior_parameters_η[2,4])
prior_dist_τ_3 = Gamma(prior_parameters_η[3,3],1/prior_parameters_η[3,4])


# check loglik est val
diff_kalman_pf = zeros(5)

PyPlot.figure()

for i in 1:5
    Random.seed!(i)

    σ_ϵ_new = rand(prior_dist_σ_ϵ)

    τ_1_new  = rand(prior_dist_τ_1)
    τ_2_new  = rand(prior_dist_τ_2)
    τ_3_new  = rand(prior_dist_τ_3)

    μ_1_new  = rand(Normal(prior_parameters_η[1,1], sqrt(1/(prior_parameters_η[1,2]*τ_1_new))))
    μ_2_new  = rand(Normal(prior_parameters_η[2,1], sqrt(1/(prior_parameters_η[2,2]*τ_2_new))))
    μ_3_new  = rand(Normal(prior_parameters_η[3,1], sqrt(1/(prior_parameters_η[3,2]*τ_3_new))))

    log_θ_1_new = μ_1_new .+ sqrt(1/τ_1_new)*randn(M_subjects)
    log_θ_2_new = μ_1_new .+ sqrt(1/τ_2_new)*randn(M_subjects)
    log_θ_3_new = μ_1_new .+ sqrt(1/τ_3_new)*randn(M_subjects)

    η_new = [μ_1_new; μ_2_new; μ_3_new; τ_1_new; τ_2_new; τ_3_new]
    ϕ_new = [log_θ_1_new log_θ_2_new log_θ_3_new]


    @time loglik_kalman, y_hat, x_hat, res_kalman = kalman(y, σ_ϵ_new, ϕ_new, dt, true)

    loglik_kalman = sum(loglik_kalman)

    loglik_pf_m = zeros(M_subjects,100)

    for i in 1:100
        loglik_pf_m[:,i] = cpf(y, σ_ϵ_new, ϕ_new, dt, randn(nbr_particles_pf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_pf, false)
    end

    loglik_pf_m = sum(loglik_pf_m, dims = 1)

    loglik_pg_dist = kde(loglik_pf_m[:])

    loglik_cpf_m = zeros(M_subjects,100)


    for i in 1:100
        loglik_cpf_m[:,i] = cpf(y, σ_ϵ_new, ϕ_new, dt, randn(nbr_particles_cpf,N_time+1,M_subjects), randn(N_time,2,M_subjects), nbr_particles_cpf, true)
    end

    loglik_cpf_m = sum(loglik_cpf_m, dims = 1)

    loglik_pg_dist_cpf = kde(loglik_cpf_m[:])

    PyPlot.subplot(2,2,i)
    PyPlot.plot(loglik_pg_dist.x,loglik_pg_dist.density, "b")
    PyPlot.plot(loglik_pg_dist_cpf.x,loglik_pg_dist_cpf.density, "--b")
    PyPlot.plot((loglik_kalman, loglik_kalman), (0, maximum(loglik_pg_dist.density)), "k")

    diff_kalman_pf[i] = abs(loglik_kalman - mean(loglik_pf_m))
end


diff_kalman_pf


PyPlot.figure()
PyPlot.plot(diff_kalman_pf, "*")

# check loglik tracking

for i in 1:10
    PyPlot.figure()
    plot_counter = 0
    for idx in 1:M_subjects

        Random.seed!(i)
        σ_ϵ_new = rand(prior_dist_σ_ϵ)
        θ_1_new = rand(prior_dist_θ_1)
        θ_3_new = rand(prior_dist_θ_3)
        ϕ_new = rand(dist_θ_2, M_subjects)

        @time loglik_kalman, partial_loglik_sum_kalman = kalman(y, [θ_1_new; θ_3_new], θ_3_new, ϕ_new, dt, false, true)

        @time loglik_pf, w, x, partial_loglik_sum_pf = pf(y, [θ_1_new; θ_3_new], θ_3_new, ϕ_new, dt, 1000, true, true)

        #idx = 10
        partial_loglik_kalman = diff(partial_loglik_sum_kalman,dims=2)
        partial_loglik_pf = diff(partial_loglik_sum_pf,dims=2)

        plot_counter = plot_counter+1
        PyPlot.subplot(10,4,plot_counter)
        PyPlot.plot(partial_loglik_pf[idx,:], "r")
        PyPlot.plot(partial_loglik_kalman[idx,:], "b")

    end

end

# check path tracking

for i in 1:10
    PyPlot.figure()
    plot_counter = 0
    for idx in 1:M_subjects

        Random.seed!(i)
        σ_ϵ_new = rand(prior_dist_σ_ϵ)
        θ_1_new = rand(prior_dist_θ_1)
        θ_3_new = rand(prior_dist_θ_3)
        ϕ_new = rand(dist_θ_2, M_subjects)

        @time loglik_kalman, partial_loglik_sum_kalman = kalman(y, [θ_1_new; θ_3_new], θ_3_new, ϕ_new, dt, false, true)

        @time loglik_pf, w, x, partial_loglik_sum_pf = pf(y, [θ_1_new; θ_3_new], θ_3_new, ϕ_new, dt, 1000, true, true)

        plot_counter = plot_counter+1
        PyPlot.subplot(10,4,plot_counter)
        PyPlot.plot(Array(x[idx,1:10,:]'), "--")
        PyPlot.plot(Array(y[idx,:]), "b")

    end
end
