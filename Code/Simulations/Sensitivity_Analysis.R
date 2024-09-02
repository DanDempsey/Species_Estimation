##### Simulation sensitivity analysis
##### Daniel Dempsey

### Load in simulation analysis functions
source( 'Code/Simulations/Simulation_Function.R' )

### Run simulation study
# Set working directory
sim_res_wd <- 'Output/Simulations/World_A/Sensitivity_Analysis'
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

sim_vals <- list( theta_N = c(4, 2.5), theta_l = c(5, 0.1), theta_x = log(30), C = C, 
                  theta_N_sim = c(4, 2.5), theta_l_sim = c(5, 0.1), theta_x_sim = log(30),
                  Tau = sim_tau, hyperprior = FALSE )

param_list <- list()
for ( i in 1:3 ) {
  bump_up <- sim_vals
  bump_down <- sim_vals
  nm <- names( sim_vals )[i]
  for ( j in 1:length(sim_vals[[i]]) ) {
    bump_up[[i]][j] <- sim_vals[[i]][j] * 1.2
    bump_up$sim_name <- paste0( nm, j, '_Up' )
    bump_down[[i]][j] <- sim_vals[[i]][j] * 0.8
    bump_down$sim_name <- paste0( nm, j, '_Down' )
    param_list <- c( param_list, list(bump_up), list(bump_down) )
  }
}

load( '../Simulation_data.Rdata' )
dat_list <- list()
for ( i in 1:length(dat) ) {
  dat_list[[i]] <- dat[[i]]$dat
}
rm( dat )
gc()

### Run simulation study
# Split into two runs to avoid memory issues
tic()
RNGkind("L'Ecuyer-CMRG") # So that parallel computing results are reproducible
set.seed( 1234 )
sim_run( param_list = param_list[1:5], dat_list = dat_list, between_loops = between_loops, 
         within_loops = within_loops, pilot_its = pits, actual_its = fits, n_cores = n_cores )
gc()
sim_run( param_list = param_list[6:10], dat_list = dat_list, between_loops = between_loops, 
         within_loops = within_loops, pilot_its = pits, actual_its = fits, n_cores = n_cores )
toc() # 75050.707 sec elapsed

### Extract results
file_names <- c( '..', substr(list.dirs(), 3, nchar(list.dirs()))[-1] )
sensitivity_performance <- list()
xrange <- numeric()
num_sims <- length( file_names )

for ( i in 1:num_sims ) {
  comp_list <- list()
  load( paste0(file_names[i], '/Simulation_data.Rdata') )
  sensitivity_performance[[file_names[i]]] <- sapply( dat, function(x) { mean( x$res$Run_1$C ) } )
  xrange <- range( c(sensitivity_performance[[file_names[i]]]), xrange )
}

par_names <- substr( file_names[-1], 1, 8 ) %>% unique

suffix <- substr( file_names, nchar(file_names)-1, nchar(file_names) )
up_inds <- which( suffix == "Up" )
down_inds <- which( suffix == "wn" )
up_list <- sensitivity_performance[up_inds]
down_list <- sensitivity_performance[down_inds]

n_reps <- length( sensitivity_performance[[1]] )
baseline_range <- range( sensitivity_performance$.. )

png( 'Mean_Distibution.png' )
plot( 0, type = 'n', xlim = xrange, ylim = c(0, (num_sims-1)/2) + 0.5, ylab = '', 
      yaxt = 'n', xlab = 'C', main = 'Distribution of Means across Sensitivity Analysis' )
axis( 2, at = 1:5, labels = par_names, las = 2 )
for( i in 1:length(up_list) ) {
  points( x = up_list[[i]], y = rep(i, n_reps) + 0.1, col = 'dodgerblue', pch = 20 )
  line_range <- range( up_list[[i]] )
  lines( x = line_range, y = c(i, i) + 0.1, col = 'dodgerblue' )
}
for( i in 1:length(down_list) ) {
  points( x = down_list[[i]], y = rep(i, n_reps) - 0.1, col = 'orange', pch = 20 )
  line_range <- range( down_list[[i]] )
  lines( x = line_range, y = c(i, i) - 0.1, col = 'orange' )
}
abline( v = baseline_range, lty = 2, col = 'darkgrey' )
dev.off()

