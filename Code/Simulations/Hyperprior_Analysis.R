##### Simulation hyperprior analysis
##### Daniel Dempsey

### Load in simulation analysis functions
source( 'Code/Simulations/Simulation_Function.R' )

### Run simulation study
# Set working directory
sim_res_wd <- 'Output/Simulations/World_A/Hyperprior_Analysis'
dir.create( sim_res_wd )
setwd( sim_res_wd )

### Set parameters
sim_tau <- 250
pits <- 1e4
fits <- 1e6
between_loops <- 6
within_loops <- 1
m <- 10
C <- 1e4
n_cores <- detectCores() - 1

sim_vals <- list( theta_l = c(5, 0.1), theta_x = log(30), C = C, 
                  theta_N_sim = c(4, 2.5), theta_l_sim = c(5, 0.1), theta_x_sim = log(30),
                  Tau = sim_tau, hyperprior = TRUE )

lower_bound_vals <- seq( 1, 0.5, length.out = 5 )
upper_bound_vals <- seq( 3, 9, length.out = 5 )

param_list <- list()
for ( i in 1:length(lower_bound_vals) ) {
  sim_vals$theta_N <- c( 4, lower_bound_vals[i], upper_bound_vals[i] )
  sim_vals$sim_name <- paste0( 'Hyperprior_', i )
  param_list <- c( param_list, list(sim_vals) )
}

load( '../Simulation_data.Rdata' )
dat_list <- list()
for ( i in 1:length(dat) ) {
  dat_list[[i]] <- dat[[i]]$dat
}
rm( dat )
gc()

### Run simulation study
tic()
RNGkind("L'Ecuyer-CMRG") # So that parallel computing results are reproducible
set.seed( 1608 )
sim_run( param_list = param_list, dat_list = dat_list, between_loops = between_loops, 
         within_loops = within_loops, pilot_its = pits, actual_its = fits, n_cores = n_cores )
toc() # 40409.857 sec elapsed

