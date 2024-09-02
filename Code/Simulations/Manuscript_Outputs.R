##### Simulation Outputs
##### Daniel Dempsey

library( magrittr )

### World A and B simulations with proper header
worlds <- list()
load( 'Output/Simulations/World_A/Simulation_data.Rdata' )
worlds[[1]] <- dat
rm( dat )
gc()
load( 'Output/Simulations/World_B/Simulation_data.Rdata' )
worlds[[2]] <- dat
rm( dat )
gc()

main_nm <- c( 'World A', 'World B' )
true_C <- c( 10000, 25000 )

png( 'Output/Simulations/Number_of_Species.png', height = 500, width = 800 )
par( mfrow = c(1, 2), mar = c(2.3, 2, 2, 2) )
for ( j in 1:length(worlds) ) {
  dat <- worlds[[j]]
  C_extract <- function( x ) {
    lapply( x$res[-1], '[[', 'C' )
  }
  C_sim <- lapply( dat, C_extract )
  C_means <- sapply( C_sim, function(x) { sapply( x, mean ) } ) %>% as.vector
  C_means_mat <- sapply( C_sim, function(x) { sapply( x, mean ) } )
  if ( is.matrix(C_means_mat) ) { 
    C_means_overall <- apply( C_means_mat, 2, mean ) %>% mean 
  }
  else {
    C_means_overall <- mean( C_means_mat )
  }
  C_dens <- lapply( C_sim, function(x) { lapply( x, density ) } ) %>% do.call( what = 'c' )
  C_dens_x <- lapply( C_dens, '[[', 'x' )
  C_dens_y <- lapply( C_dens, '[[', 'y' )
  xlim <- lapply( C_dens_x, range ) %>% do.call( what = 'rbind' ) %>% range
  ylim <- lapply( C_dens_y, range ) %>% do.call( what = 'rbind' ) %>% range
  col_code <- rep( 1:length(C_sim), each = length(C_sim[[1]]) ) + 1
  
  plot( 0, type = 'n', xlim = xlim, ylim = ylim, yaxt = 'n', ylab = '', xlab = '',
        main = main_nm[j] )
  for ( i in 1:length(C_dens_x) ) {
    lines( C_dens_x[[i]], C_dens_y[[i]], col = col_code[i] )
  }
  p_pos <- rep( seq(0, 0.25, length.out = length(C_sim)), length(C_sim[[1]]) )
  y_pos <- quantile( ylim, p_pos )
  abline( v = true_C[j], lty = 2, lwd = 2 )
  points( C_means, y_pos, pch = 20, col = col_code )
  abline( v = C_means_overall, lty = 2, col = 'grey', lwd = 2 )
}
dev.off()

### Nicer hyperprior plots
worlds <- list()
setwd( 'Output/Simulations/World_A/Hyperprior_Analysis/' )
sets <- paste0( 'Hyperprior_', 1:5 )
for ( i in 1:length(sets) ) {
  load( paste0(sets[i], '/Simulation_data.Rdata') )
  worlds[[i]] <- dat
  rm( dat )
  gc()
}

lower_bound_vals <- seq( 1, 0.5, length.out = 5 )
upper_bound_vals <- seq( 3, 9, length.out = 5 )

main_nm <- paste0( lower_bound_vals, ' to ', upper_bound_vals )

png( 'All_Hyperpriors.png', height = 500, width = 800 )
par( mfrow = c(2, 3), mar = c(2.3, 2, 2, 2) )
for ( j in 1:length(worlds) ) {
  dat <- worlds[[j]]
  C_extract <- function( x ) {
    lapply( x$res[-1], '[[', 'C' )
  }
  C_sim <- lapply( dat, C_extract )
  C_means <- sapply( C_sim, function(x) { sapply( x, mean ) } ) %>% as.vector
  C_means_mat <- sapply( C_sim, function(x) { sapply( x, mean ) } )
  if ( is.matrix(C_means_mat) ) { 
    C_means_overall <- apply( C_means_mat, 2, mean ) %>% mean 
  }
  else {
    C_means_overall <- mean( C_means_mat )
  }
  C_dens <- lapply( C_sim, function(x) { lapply( x, density ) } ) %>% do.call( what = 'c' )
  C_dens_x <- lapply( C_dens, '[[', 'x' )
  C_dens_y <- lapply( C_dens, '[[', 'y' )
  xlim <- lapply( C_dens_x, range ) %>% do.call( what = 'rbind' ) %>% range
  ylim <- lapply( C_dens_y, range ) %>% do.call( what = 'rbind' ) %>% range
  col_code <- rep( 1:length(C_sim), each = length(C_sim[[1]]) ) + 1
  
  plot( 0, type = 'n', xlim = xlim, ylim = ylim, yaxt = 'n', ylab = '', xlab = '',
        main = main_nm[j] )
  for ( i in 1:length(C_dens_x) ) {
    lines( C_dens_x[[i]], C_dens_y[[i]], col = col_code[i] )
  }
  #p_pos <- rep( seq(0, 0.25, length.out = length(C_sim)), length(C_sim[[1]]) )
  #y_pos <- quantile( ylim, p_pos )
  abline( v = 10000, lty = 2, lwd = 2 )
  #points( C_means, y_pos, pch = 20, col = col_code )
  #abline( v = C_means_overall, lty = 2, col = 'grey', lwd = 2 )
}
dev.off()

### Sensitivity analysis
worlds <- list()
setwd( 'Output/Simulations/World_A/Sensitivity_Analysis/' )

for ( i in 1:length(sets) ) {
  load( paste0(sets[i], '/Simulation_data.Rdata') )
  worlds[[i]] <- dat
  rm( dat )
  gc()
}

file_names <- c( '..', list.dirs( full.names = FALSE )[-1] )
sensitivity_performance <- list()
xrange <- numeric()
num_sims <- length( file_names )

for ( i in 1:num_sims ) {
  load( paste0(file_names[i], '/Simulation_data.Rdata') )
  sensitivity_performance[[file_names[i]]] <- lapply( dat, function(x) { 
    quantile( x$res$Run_1$C, c(0.025, 0.975) ) } )
  xrange <- range( c(sensitivity_performance[[file_names[i]]]), xrange )
}

par_names <- substr( file_names[-1], 1, 8 ) %>% unique

suffix <- substr( file_names, nchar(file_names)-1, nchar(file_names) )
up_inds <- which( suffix == "Up" )
down_inds <- which( suffix == "wn" )
up_list <- sensitivity_performance[up_inds]
down_list <- sensitivity_performance[down_inds]

n_reps <- length( sensitivity_performance[[1]] )
baseline_range <- sensitivity_performance[[1]]

png( 'Sensitivity_CIs.png', height = 500, width = 800 )
par( mfrow = c(2, 3), mar = c(2.3, 2, 2, 2) )
for ( i in 1:n_reps ) {
  use_up_dat <- lapply( up_list, '[[', i )
  use_down_dat <- lapply( down_list, '[[', i )
  baseline_dat <- baseline_range[[i]]
  xrange <- range( use_up_dat, use_down_dat, baseline_dat )
  plot( 0, type = 'n', xlim = xrange, ylim = c(0, (num_sims-1)/2) + 0.5, ylab = '', 
        yaxt = 'n', xlab = 'C', main = '' )
  axis( 2, at = 1:5, labels = expression(mu[L[1]], sigma[L]^2, mu[N], sigma[N]^2, theta[x]), 
        las = 2 )
  abline( v = 1e4, lty = 2 )
  for( j in 1:length(use_up_dat) ) {
    points( x = use_up_dat[[j]], y = rep(j, 2) + 0.1, col = 'dodgerblue', pch = 20 )
    line_range <- range( use_up_dat[[j]] )
    lines( x = line_range, y = c(j, j) + 0.1, col = 'dodgerblue' )
  }
  for( j in 1:length(use_down_dat) ) {
    points( x = use_down_dat[[j]], y = rep(j, 2) - 0.1, col = 'orange', pch = 20 )
    line_range <- range( use_down_dat[[j]] )
    lines( x = line_range, y = c(j, j) - 0.1, col = 'orange' )
  }
  abline( v = baseline_dat, lty = 2, col = 'darkgrey' )
}
dev.off()
