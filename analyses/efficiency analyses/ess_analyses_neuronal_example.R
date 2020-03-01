library(mcmcse)
library(readr)

setwd("~/Documents/projects/cpmmh for sdemems/code")

load_results <- function(seed) {
  
  R = 100000
  burn_in = 20000
  nbr_na = 0
  
  sim_res <- matrix(nrow = 4, ncol = 3)
  
  M_subjects = "100"
  
  # kalman
  job = paste(seed, "_", M_subjects, sep = "")
  
  simdata_kalman <- read_csv(paste("data/SDEMEM OU neuron data/kalman/sim_data_", job, ".csv", sep = ""))
  run_time_kalman = as.numeric((simdata_kalman)[1,1])/60
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU neuron data/kalman/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU neuron data/kalman/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU neuron data/kalman/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:300]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_kalman = min(na.omit(ess))/run_time_kalman
  
  
  # pmmh
  nbr_particles = "1"
  rho = "0.0"
  
  job = paste(seed, "_", M_subjects, "_", rho,"_", nbr_particles, "_", "bridge", sep = "")
  
  
  simdata_pmmh <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_pmmh = as.numeric((simdata_pmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:300]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_pmmh = min(na.omit(ess))/run_time_pmmh
  
  # CPMMH bridge 0999
  nbr_particles = "1"
  rho = "0.999"
  
  job = paste(seed, "_", M_subjects, "_", rho,"_", nbr_particles, "_", "bridge", sep = "")
  
  simdata_cpmmh <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_cpmmh = as.numeric((simdata_cpmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:300]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_cpmmh = min(na.omit(ess))/run_time_cpmmh

  
  
  # CPMMH bridge 09
  nbr_particles = "1"
  rho = "0.9"
  
  job = paste(seed, "_", M_subjects, "_", rho,"_", nbr_particles, "_", "bridge", sep = "")
  
  simdata_cpmmh09 <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_cpmmh09 = as.numeric((simdata_cpmmh09)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU neuron data/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:300]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_cpmmh09 = min(na.omit(ess))/run_time_cpmmh09
  
  
  sim_res[,1] <- c(run_time_kalman,
                   run_time_pmmh,
                   run_time_cpmmh09,
                   run_time_cpmmh)
  
  sim_res[,2] <- c(ess_min_per_min_kalman*run_time_kalman, 
                   ess_min_per_min_pmmh*run_time_pmmh,
                   ess_min_per_min_cpmmh09*run_time_cpmmh09,
                   ess_min_per_min_cpmmh*run_time_cpmmh)
  
  
  sim_res[,3] <- c(ess_min_per_min_kalman,
                   ess_min_per_min_pmmh,
                   ess_min_per_min_cpmmh09,
                   ess_min_per_min_cpmmh)
  
  
  #sim_res[,4] <- c(ess_min_per_min_kalman/ess_min_per_min_pmmh_naive, 
  #                 ess_min_per_min_cpmmh099/ess_min_per_min_pmmh_naive,
  #                 ess_min_per_min_cpmmh0999/ess_min_per_min_pmmh_naive)
  
  
  
  
  # change to per min and runtime in mie 
  print(paste0("Seed: ", seed))
  
  # print results
  print("Runtime (min) ")
  print(paste0("Kalman: ", run_time_kalman))
  print(paste0("PMMH: ", run_time_pmmh))
  print(paste0("CPMMH-09: ", run_time_cpmmh09))
  print(paste0("CPMMH-0999: ", run_time_cpmmh))
  
  
  print("min ESS (for all chains and blocks)")
  print(paste0("Kalman: ", ess_min_per_min_kalman*run_time_kalman))
  print(paste0("PMMH: ", ess_min_per_min_pmmh*run_time_pmmh))
  print(paste0("CPMMH-09: ", ess_min_per_min_cpmmh09*run_time_cpmmh09))
  print(paste0("CPMMH-099: ", ess_min_per_min_cpmmh*run_time_cpmmh))
  
  print("min ESS/min (for all chains and blocks)")
  print(paste0("Kalman: ", ess_min_per_min_kalman))
  print(paste0("PMMH: ", ess_min_per_min_pmmh))
  print(paste0("CPMMH-09: ", ess_min_per_min_cpmmh09))
  print(paste0("CPMMH-0999: ", ess_min_per_min_cpmmh))
  
  print("Rel")
  print(paste0("Kalman: ", ess_min_per_min_kalman/ess_min_per_min_pmmh))
  print(paste0("PMMH: ", ess_min_per_min_pmmh/ess_min_per_min_pmmh))
  print(paste0("CPMMH-09: ", ess_min_per_min_cpmmh09/ess_min_per_min_pmmh))
  print(paste0("CPMMH-0999: ", ess_min_per_min_cpmmh/ess_min_per_min_pmmh))
  
  print("Nbr NA")
  print(paste0("Nbr NA in ess:s: ", nbr_na))
  
  return(sim_res)
  
}


res = load_results(1*100)

print(round(res,2))


res = array(dim=c(5,4,3))

for (dataset in 1:5){ 
  res[dataset,,] <-load_results(dataset*100)
}

# print mean over the five data sets
print(round(apply(res,c(2,3),mean),2))
