library(mcmcse)
library(readr)

setwd("~/Documents/projects/cpmmh for sdemems/code")

# function for loading the full results from one experiment 
load_results <- function(seed, nbr_particles = "5", rho = "0.999") {

  R = 60000
  burn_in = 10000
  nbr_na = 0
  
  sim_res <- matrix(nrow = 1, ncol = 5)
  
  N_time = "200"
  M_subjects = "40"
  
  # load restuls 
  
  job = paste(seed, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")
  
  simdata_cpmmh <- read_csv(paste("data/SDEMEM OU/cpmmh for plot mess vs N/sim_data_", job, ".csv", sep = ""))
  run_time_cpmmh = as.numeric((simdata_cpmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/cpmmh for plot mess vs N/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/cpmmh for plot mess vs N/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/cpmmh for plot mess vs N/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)

  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_cpmmh = min(na.omit(ess))/run_time_cpmmh
  ess_avg_per_min_cpmmh = mean(na.omit(ess))/run_time_cpmmh
  
  sim_res = c(run_time_cpmmh, 
              ess_min_per_min_cpmmh*run_time_cpmmh,
              ess_min_per_min_cpmmh,
              nbr_na)
  
  return(list(sim_res))
  
}



output_test = load_results(200, 5, 0.99)

seeds = c(100,200,300,400,500,600,700,800,900,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)

N = c(5,10,20,50,100)


seeds = c(100,200,300,400,500)

N = c(100)



data_temp = matrix(nrow = length(seeds), ncol = 4)
data_temp[2,] = unlist(load_results(200, 5, 0.99), use.names=FALSE) 




for (n in N){
    
  data_temp = matrix(nrow = length(seeds), ncol = 4)
  
  rhp_level = 0.999
  N_val = n
  
  for (i in 1:length(seeds)){ 
    print(i)    
    data_temp[i,] = unlist(load_results(seeds[i], N_val, rhp_level), use.names=FALSE)
  }
  
  
  save_path = paste('./data/SDEMEM OU/cpmmh for plot mess vs N/', 'run_time_ess_mess', '_', 
                    toString(rhp_level), '_', toString(N_val), '.csv', sep = "")
  
  #write.csv(data.frame(data_temp), save_path)

}


# todo loop over N 
# save one file for each rho_level and N combination w all runs for that point

