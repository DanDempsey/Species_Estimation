##### Run simulation analysis
##### Daniel Dempsey

### Load in simulation analysis functions
source( 'Code/Simulations/Simulation_Function.R' )

### Run simulation study
# Set working directory
sim_res_wd <- 'Output/Simulations'
dir.create( sim_res_wd )
setwd( sim_res_wd )

### Set parameters
sim_tau <- 250
pits <- 1e4
fits <- 1e6
between_loops <- 6
within_loops <- 2
m <- 10
n_cores <- detectCores() - 1

worldA <- list( C = 1e4, theta_N = c(4, 2.5), theta_N_sim = c(4, 2.5), theta_l = c(5, 0.1), theta_l_sim = c(5, 0.1),
                theta_x = log(30), theta_x_sim = log(30), Tau = sim_tau, hyperprior = FALSE, sim_name = 'World A' )

worldB <- list( C = 5e4/2, theta_N = c(6, 2), theta_N_sim = c(6, 2), theta_l = c(7, 0.08), theta_l_sim = c(7, 0.08),
                theta_x = log(50), theta_x_sim = log(50), Tau = sim_tau, hyperprior = FALSE, sim_name = 'World B' )

param_list <- list( worldA, worldB )

### Run simulation study
tic()
RNGkind("L'Ecuyer-CMRG") # So that parallel computing results are reproducible
sim_run( param_list = param_list, between_loops = between_loops, within_loops = within_loops, 
         pilot_its = pits, actual_its = fits, dat_seed = 2024, n_cores = n_cores )
toc() # 84263.485 sec elapsed


