##### Functions for running simulation analyses
##### Daniel Dempsey

### Load required code
source( 'Code/Main_Software/Class_Discovery_ABC.R' )
library( scales )
library( tictoc )

### Utility function for computing discovery times
t_val <- function( d, n_sum ) {
  sum( n_sum < d ) + 1
}

### Simulation function 
sim_dat <- function( C, theta_N, theta_l, theta_x, Tau, graph = FALSE ) {
  
  # Latent effort L
  L <- numeric( Tau )
  L[1] <- rnorm( 1, theta_l[1], theta_l[2] )
  for ( i in 2:Tau ) {
    L[i] <- rnorm( 1, L[i-1], theta_l[2] )
  }
  
  # Class abundance N
  N <- rlnorm( C, theta_N[1], theta_N[2] ) %>% sort( decreasing = TRUE )
  N_inds <- sample( 1:C, C, prob = N )
  N <- N[N_inds]
  N_sum_prop <- cumsum( N[1:(C-1)] ) / sum( N )
  
  # Inter-discovery times d
  d <- rep( 1, C )
  d[2:C] <- rgeom( n = C-1, 1 - N_sum_prop )
  d_sum <- cumsum( d )
  
  # Sample size n
  n <- rpois( Tau, exp(L) )
  n_sum <- cumsum( as.numeric(n) )
  
  # Simulate observed data x and t
  x <- rpois( Tau, exp(L - theta_x) )
  t <- sapply( d_sum, t_val, n_sum )
  x_sim <- x[1:Tau]
  t_sim <- as.numeric( table(t[t <= Tau] %>% factor(levels = 1:Tau)) )
  
  # Plot the data if desired
  if ( graph ) {
    par( mfrow = c(2, 1) )
    plot( x_sim, type = 'l', main = 'Effort Proxy over Time', xlab = 'Time', ylab = 'Effort Proxy' )
    plot( t_sim, type = 'l', main = 'Number of Discoveries over Time', xlab = 'Time', 
          ylab = 'Number of Discoveries' )
    par( mfrow = c(1, 1) )
  }
  
  # Return simulated data
  list( x = x_sim, t = t_sim )
  
}

### Function that runs ABC
sim_analysis <- function( dat, theta_N, theta_l, theta_x, m, hyperprior, 
                          pilot_its, actual_its, loops, n_cores ) {
  
  # Extract observed data
  x_sim <- dat$x
  t_sim <- dat$t
  
  # Run pilot
  pilot_res <- CD_ABC( x = x_sim, t = t_sim, theta_N = theta_N, theta_x = theta_x, 
                       theta_l = theta_l, m = m, iters = pilot_its, epsilon = Inf, 
                       abundance_hyperprior = hyperprior, cores = n_cores )
  
  sim_epsilon <- quantile( pilot_res$distance, 0.1 )
  
  # Full simulation
  all_res <- list()
  all_res$Pilot <- pilot_res
  seeds <- sample( 1:1e6, loops )
  for ( i in 1:loops + 1 ) {
    
    set.seed( seeds[i - 1] )  # We need to do this since L'Ecuyer-CMRG does not reset seed automatically
    all_res[[i]] <- CD_ABC( x = x_sim, t = t_sim, theta_N = theta_N, theta_x = theta_x, 
                            theta_l = theta_l, m = m, iters = actual_its, epsilon = sim_epsilon, 
                            abundance_hyperprior = hyperprior, cores = n_cores )
    
  }
  
  names( all_res )[1:loops + 1] <- paste0( 'Run_', 1:loops )
  all_res
    
}

### Function that creates outputs
sim_output <- function( dat, param_list ) {
  
  n_reps <- length( dat )
  
  # Create Directory for output
  current_wd <- getwd()
  wd <- gsub( ' ', '_', param_list$sim_name )
  suppressWarnings( dir.create(wd) )
  setwd( wd )
  
  # Save data
  save( dat, file = 'Simulation_data.Rdata' )
  
  # Plot C
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
  
  png( 'Number_of_Species.png' )
  plot( 0, type = 'n', xlim = xlim, ylim = ylim, yaxt = 'n', ylab = '', xlab = 'C' )
  for ( i in 1:length(C_dens_x) ) {
    lines( C_dens_x[[i]], C_dens_y[[i]], col = col_code[i] )
  }
  p_pos <- rep( seq(0, 0.25, length.out = length(C_sim)), length(C_sim[[1]]) )
  y_pos <- quantile( ylim, p_pos )
  abline( v = param_list$C, lty = 2, lwd = 2 )
  points( C_means, y_pos, pch = 20, col = col_code )
  abline( v = C_means_overall, lty = 2, col = 'grey', lwd = 2 )
  dev.off()
  
  # Plot effort proxy
  Tau <- param_list$Tau
  for ( i in 1:length(C_sim) ) {
    png( paste0('Effort_Proxy_Sim', i, '.png') )
    x_res <- dat[[i]]$res$Run_1$x_hat
    x_range <- range( x_res )
    plot( 0, type = 'n', xlim = c(1, Tau), ylim = x_range, xlab = 'Time', ylab = 'Effort Proxy' )
    for ( j in 1:nrow(x_res) ) {
      lines( x_res[j, ], col = alpha('grey', 0.5) )
    }
    lines( dat[[i]]$dat$x, lwd = 1.5 )
    dev.off()
  }
  
  # Plot discovery times
  for ( i in 1:length(C_sim) ) {
    png( paste0('Discovery_Time_Sim', i, '.png') )
    t_res <- dat[[i]]$res$Run_1$t_hat
    t_range <- range( t_res )
    plot( 0, type = 'n', xlim = c(1, Tau), ylim = t_range, xlab = 'Time', ylab = 'Effort Proxy' )
    for ( j in 1:nrow(t_res) ) {
      lines( t_res[j, ], col = alpha('grey', 0.5) )
    }
    lines( dat[[i]]$dat$t, lwd = 1.5 )
    dev.off()
  }
  
  setwd( current_wd )
  return( NULL )
  
}

### Combine all functions
sim_run <- function( param_list, dat_list, between_loops, within_loops, pilot_its, actual_its, dat_seed, n_cores ) {
  
  all_sims <- list()
  analysis_len <- length( param_list )
  for ( i in 1:analysis_len ) {
    
    cat( paste0('Now running analysis for set ', i, ' of ', analysis_len, '.\n') )
    all_sims[[i]] <- list()
    
    for ( j in 1:between_loops ) {
      cat( paste0('Repetition ', j, ' of ', between_loops, '...\n') )
      list2env( param_list[[i]], environment() )
      
      # Generate simulation data
      if ( missing(dat_list) ) {
        if ( !missing(dat_seed) ) { set.seed( dat_seed + j ) }
        dat <- sim_dat( C = C, theta_N = theta_N_sim, theta_l = theta_l_sim, 
                        theta_x = theta_x_sim, Tau = Tau )
      }
      else {
        dat <- dat_list[[j]]
      }
      
      # Fit ABC
      res <- sim_analysis( dat, theta_N = theta_N, theta_l = theta_l, 
                           theta_x = theta_x, m = m, hyperprior = hyperprior, 
                           pilot_its = pilot_its, actual_its = actual_its, 
                           loops = within_loops, n_cores = n_cores )
      
      all_sims[[i]][[j]] <- list( dat = dat, res = res )
      gc()
    }
    
  }
  
  gc()
  cat( 'All simulation analyses complete! Now producing output... ' )
  
  # Produce output
  Map( sim_output, dat = all_sims, param_list = param_list )
  
  cat( 'complete!\n' )
  return( invisible(NULL) )
  
}

