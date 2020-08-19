library(mcmcse)
library(readr)

#setwd("~/Documents/projects/cpmmh for sdemems/code")
setwd("~/Documents/projects/cpmmh for sdemems/public github/efficient_SDEMEM")

load_results <- function(seed) {

  R = 60000
  burn_in = 10000
  nbr_na = 0
  
  sim_res <- matrix(nrow = 5, ncol = 4)
  
  N_time = "200"
  M_subjects = "40"
  
  # kalman
  job = paste(seed, "_", M_subjects, "_", N_time, sep = "")
  
  simdata_kalman <- read_csv(paste("data/SDEMEM OU/kalman/sim_data_", job, ".csv", sep = ""))
  run_time_kalman = as.numeric((simdata_kalman)[1,1])/60
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/kalman/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/kalman/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/kalman/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_kalman = min(na.omit(ess))/run_time_kalman
  
  
  # pmmh - naive 
  
  nbr_particles = "3000"
  rho = "0.0"
  
  job = paste(seed, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")
  
  simdata_pmmh <- read_csv(paste("data/SDEMEM OU/naive_cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_pmmh_naive = as.numeric((simdata_pmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/naive_cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/naive_cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/naive_cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_pmmh_naive = min(na.omit(ess))/run_time_pmmh_naive
  
  
  # pmmh
  nbr_particles = "3000"
  rho = "0.0"
  
  job = paste(seed, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")
  
  simdata_pmmh <- read_csv(paste("data/SDEMEM OU/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_pmmh = as.numeric((simdata_pmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_pmmh = min(na.omit(ess))/run_time_pmmh
  
  
  # cpmmh-099
  nbr_particles = "100"
  rho = "0.99"
  
  job = paste(seed, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")
  
  simdata_pmmh <- read_csv(paste("data/SDEMEM OU/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_cpmmh099 = as.numeric((simdata_pmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_cpmmh099 = min(na.omit(ess))/run_time_cpmmh099
  
  
  # cpmmh-0999
  nbr_particles = "50"
  rho = "0.999"
  
  job = paste(seed, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")
  
  simdata_pmmh <- read_csv(paste("data/SDEMEM OU/cpmmh/sim_data_", job, ".csv", sep = ""))
  run_time_cpmmh0999 = as.numeric((simdata_pmmh)[1,1])/60
  
  
  chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
  chain = chain_sigma_epsilon[(burn_in+1):R,1]
  ess_epsilon = ess(chain)
  
  chain_eta <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_eta_", job, ".csv", sep = ""))
  chain = chain_eta[(burn_in+1):R,1:6]
  ess_eta = ess(chain)
  
  chain_phi <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_phi_", job, ".csv", sep = ""))
  chain = chain_phi[(burn_in+1):R,1:120]
  ess_phi = ess(chain)
  
  ess = c(ess_epsilon, ess_eta, ess_phi)
  nbr_na = nbr_na + sum(is.na(ess))
  ess_min_per_min_cpmmh0999 = min(na.omit(ess))/run_time_cpmmh0999

  sim_res[,1] <- c(run_time_kalman, 
                   run_time_pmmh_naive, 
                   run_time_pmmh, 
                   run_time_cpmmh099, 
                   run_time_cpmmh0999)
  
  sim_res[,2] <- c(ess_min_per_min_kalman*run_time_kalman, 
                   ess_min_per_min_pmmh_naive*run_time_pmmh_naive, 
                   ess_min_per_min_pmmh*run_time_pmmh, 
                   ess_min_per_min_cpmmh099*run_time_cpmmh099, 
                   ess_min_per_min_cpmmh0999*run_time_cpmmh0999)
  
  
  sim_res[,3] <- c(ess_min_per_min_kalman, 
                   ess_min_per_min_pmmh_naive, 
                   ess_min_per_min_pmmh, 
                   ess_min_per_min_cpmmh099, 
                   ess_min_per_min_cpmmh0999)
  
  sim_res[,4] <- c(ess_min_per_min_kalman/ess_min_per_min_pmmh_naive, 
                   ess_min_per_min_pmmh_naive/ess_min_per_min_pmmh_naive, 
                   ess_min_per_min_pmmh/ess_min_per_min_pmmh_naive, 
                   ess_min_per_min_cpmmh099/ess_min_per_min_pmmh_naive,
                   ess_min_per_min_cpmmh0999/ess_min_per_min_pmmh_naive)
  
  
  
  
  # change to per min and runtime in mie 
  print(paste0("Seed: ", seed))
  
  # print results
  print("Runtime (min) ")
  print(paste0("Kalman: ", run_time_kalman))
  print(paste0("PMMH - naive: ", run_time_pmmh_naive))
  print(paste0("PMMH: ", run_time_pmmh))
  print(paste0("CPMMH-099: ", run_time_cpmmh099))
  print(paste0("CPMMH-0999: ", run_time_cpmmh0999))
  
  print("min ESS (for all chains and blocks)")
  print(paste0("Kalman: ", ess_min_per_min_kalman*run_time_kalman))
  print(paste0("PMMH - naive: ", ess_min_per_min_pmmh_naive*run_time_pmmh_naive))
  print(paste0("PMMH: ", ess_min_per_min_pmmh*run_time_pmmh))
  print(paste0("CPMMH-099: ", ess_min_per_min_cpmmh099*run_time_cpmmh099))
  print(paste0("CPMMH-0999: ", ess_min_per_min_cpmmh0999*run_time_cpmmh0999))
  
  
  print("min ESS/min (for all chains and blocks)")
  print(paste0("Kalman: ", ess_min_per_min_kalman))
  print(paste0("PMMH - naive: ", ess_min_per_min_pmmh_naive))
  print(paste0("PMMH: ", ess_min_per_min_pmmh))
  print(paste0("CPMMH-099: ", ess_min_per_min_cpmmh099))
  print(paste0("CPMMH-0999: ", ess_min_per_min_cpmmh0999))
  
  print("Rel")
  print(paste0("Kalman: ", ess_min_per_min_kalman/ess_min_per_min_pmmh_naive))
  print(paste0("PMMH: ", ess_min_per_min_pmmh_naive/ess_min_per_min_pmmh_naive))
  print(paste0("PMMH: ", ess_min_per_min_pmmh/ess_min_per_min_pmmh_naive))
  print(paste0("CPMMH-099: ", ess_min_per_min_cpmmh099/ess_min_per_min_pmmh_naive))
  print(paste0("CPMMH-0999: ", ess_min_per_min_cpmmh0999/ess_min_per_min_pmmh_naive))
  
  print("Nbr NA")
  print(paste0("Nbr NA in ess:s: ", nbr_na))
  
  
  
  
  return(sim_res)
  
}


res = array(dim=c(5,5,4))

for (dataset in 1:5){ 
  res[dataset,,] <-load_results(dataset*100)
}

# print mean over the five data sets
print(round(apply(res,c(2,3),mean),2))

print(round(apply(res,c(2,3),median),2))


print(round(apply(res,c(2,3),mean),2))

# efficiency for pmmmh-100

R = 60000
burn_in = 10000
nbr_na = 0

sim_res <- matrix(nrow = 5, ncol = 4)

N_time = "200"
M_subjects = "40"


nbr_particles = "100"
rho = "0.0"

job = paste(100, "_", M_subjects, "_", N_time,"_", nbr_particles,"_", rho, sep = "")

simdata_pmmh <- read_csv(paste("data/SDEMEM OU/cpmmh/sim_data_", job, ".csv", sep = ""))
run_time_pmmh = as.numeric((simdata_pmmh)[1,1])/60


chain_sigma_epsilon <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_sigma_epsilon_", job, ".csv", sep = ""))
chain = chain_sigma_epsilon[(burn_in+1):R,1]
ess_epsilon = ess(chain)

chain_eta <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_eta_", job, ".csv", sep = ""))
chain = chain_eta[(burn_in+1):R,1:6]
ess_eta = ess(chain)

chain_phi <- read_csv(paste("data/SDEMEM OU/cpmmh/chain_phi_", job, ".csv", sep = ""))
chain = chain_phi[(burn_in+1):R,1:120]
ess_phi = ess(chain)

ess = c(ess_epsilon, ess_eta, ess_phi)
nbr_na = nbr_na + sum(is.na(ess))
ess_min_per_min_pmmh = min(na.omit(ess))/run_time_pmmh

print("Runtime (min) ")
print(paste0("PMMH: ", run_time_pmmh))
print(paste0("PMMH: ", ess_min_per_min_pmmh*run_time_pmmh))
print(paste0("PMMH: ", ess_min_per_min_pmmh))
