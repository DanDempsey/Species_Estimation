##### Code for Running Analysis on Real Data
##### Daniel Dempsey

### Load in required libraries
library( tictoc )
library( readxl )
library( dplyr )
library( scales )
library( stringr )
source( 'Code/Main_Software/Class_Discovery_ABC.R' )

### Read in taxonomy data
read_taxa <- function( z ) {
  
  unzip( zipfile = paste0('Data/', z), files = 'data.xlsx', exdir = 'Data' )
  taxa_dat <- read_excel( 'Data/data.xlsx' ) %>% filter( .data[['col:rank']] == 'species' ) %>% 
    select( 'col:status', 'col:authorship' ) %>% suppressWarnings
  file.remove( 'Data/data.xlsx' )
  
  taxa_dat
  
}

### Utility function for replacing author synonyms
synonomise <- function( x ) {
  for ( i in seq_along(name_syn_list) ) {
    inds <- which( x %in% name_syn_list[[i]] )
    x[inds] <- name_syn_list[[i]][1]
  }
  x
}

### Clean taxonomy data
construct_ABC_data <- function( dat, synonym_list ) {
  
  # Remove missing values
  dat_complete <- na.omit( dat )
  
  # Clean authorship
  status_dat <- dat_complete[[1]]
  author_dat <- dat_complete[[2]]
  no_brackets <- gsub( "[()]", "", author_dat )
  
  name_vec <- gsub( '[[:digit:]]+', '', no_brackets )
  year_list <- str_extract_all( no_brackets, "\\d{4}" )
  drop_inds <- which( sapply(year_list, length) == 0 ) # Check for missing years
  
  if ( length(drop_inds) > 0 ) {
    status_dat <- status_dat[-drop_inds]
    name_vec <- name_vec[-drop_inds]
    year_list[drop_inds] <- NULL
  }
  
  year_vec_num <- sapply( year_list, '[', 1 ) %>% as.numeric # Take first year when two years are provided
  drop_inds2 <- which( (year_vec_num < 1758) | (year_vec_num > 2019) ) # Remove exceptions
  
  if( length(drop_inds2) > 0 ) {
    status_dat <- status_dat[-drop_inds2]
    name_vec <- name_vec[-drop_inds2]
    year_vec_num <- year_vec_num[-drop_inds2]
  }
  
  accept_inds <- which( status_dat == 'accepted' )
  
  year_vec_num_accepted <- year_vec_num[ accept_inds ]
  year_vec <- factor( year_vec_num_accepted, 
                      levels = min(year_vec_num_accepted):max(year_vec_num_accepted) )
  
  names( name_vec ) <- year_vec_num_accepted
  clean_step1 <- strsplit( name_vec, "[,&]+" ) # Split out commas and ampersands
  clean_step2 <- lapply( clean_step1, gsub, pattern = "\\ in .*", replacement = "" ) # Remove everything after "in"
  clean_step3 <- lapply( clean_step2, function(x) { gsub(' ', '', x) } ) # Remove whitespace
  clean_step_final <- lapply( clean_step3, function(x) { x[!(x=="")] } ) # Remove elements of pure whitespace
  
  if ( !missing(synonym_list) ) {
    clean_step_final <- lapply( clean_step_final, synonomise )
  }
  
  # Group the names by year
  name_yearSplit_raw <- with(stack(clean_step_final), split(values, ind)) %>% lapply( unique ) # Group list by year
  missing_years <- levels( year_vec )[ !(levels( year_vec ) %in% names( name_yearSplit_raw )) ]
  missing_year_list <- list()
  for ( i in seq_along(missing_years) ) {
    missing_year_list[[i]] <- character(0)
  }
  names( missing_year_list ) <- missing_years
  name_yearSplit_complete <- c( name_yearSplit_raw, missing_year_list ) # Add in missing years
  name_yearSplit <- name_yearSplit_complete[order(names(name_yearSplit_complete), decreasing = TRUE)] %>% 
    rev # Sort
  name_yearSplit <- lapply( name_yearSplit, sort ) # Just for readability when debugging
  
  ### Compile ABC inputs
  x <- sapply( name_yearSplit, length )
  t <- table( year_vec ) %>% as.vector
  
  cat( paste0(sum(t), ' species.\n') ) 
  list( x = x, t = t, name_list = name_yearSplit )
  
}

### Calculate time taken
time_taken <- function( x ) {
  round( x$toc - x$tic, 2 )
}

### Fit the data via ABC
run_ABC <- function( dat, param_list, metric = 'm', p = 2, pilot_iterations, 
                     fit_iterations, fixed_run = TRUE, hyperprior_run = FALSE,
                     x_weight = 1, t_weight = 1, n_cores = detectCores() - 1 ) {
  
  # Extract parameters
  x <- dat$x
  t <- dat$t
  
  theta_N <- param_list$theta_N 
  theta_N_hyperprior <- param_list$theta_N_hyperprior
  theta_l <- param_list$theta_l 
  theta_x <- param_list$theta_x 
  m <- param_list$m 
  
  res_list <- list()
  
  # Fixed prior run
  if ( fixed_run ) {
    
    # Pilot run
    tic()
    pilot_run <- CD_ABC( x = x, t = t, theta_N = theta_N, theta_x = theta_x, theta_l = theta_l, 
                         m = m, iters = pilot_iterations, epsilon = Inf, metric = metric, p = p,
                         x_weight = x_weight, t_weight = t_weight, cores = n_cores )
    pilot_speed <- toc( quiet = TRUE )
    cat( paste0('Pilot: ', time_taken(pilot_speed), ' seconds.\n') )
    
    sim_epsilon <- quantile( pilot_run$distance, 0.1 )
    
    # Actual fit
    tic()
    res <- CD_ABC( x = x, t = t, theta_N = theta_N, theta_x = theta_x, theta_l = theta_l, 
                   m = m, iters = fit_iterations, epsilon = sim_epsilon, metric = metric, p = p,
                   cores = n_cores )
    res_speed <- toc( quiet = TRUE )
    cat( paste0('Actual fit: ', time_taken(res_speed), ' seconds.\n') )
    
    res_list$fixed <- res
    
  }
  
  # Hyperprior run
  if ( hyperprior_run ) {
    
    # Pilot run
    tic()
    pilot_run <- CD_ABC( x = x, t = t, theta_N = theta_N_hyperprior, theta_x = theta_x, theta_l = theta_l, 
                         m = m, iters = pilot_iterations, epsilon = Inf, abundance_hyperprior = TRUE,
                         metric = metric, p = p, cores = n_cores )
    pilot_speed <- toc( quiet = TRUE )
    cat( paste0('Hyperprior pilot: ', time_taken(pilot_speed), ' seconds.\n') )
    
    sim_epsilon <- quantile( pilot_run$distance, 0.1 )
    
    # Actual fit
    tic()
    res <- CD_ABC( x = x, t = t, theta_N = theta_N_hyperprior, theta_x = theta_x, theta_l = theta_l, 
                   m = m, iters = fit_iterations, epsilon = sim_epsilon, abundance_hyperprior = TRUE, 
                   metric = metric, p = p, cores = n_cores )
    res_speed <- toc( quiet = TRUE )
    cat( paste0('Hyperprior fit: ', time_taken(res_speed), ' seconds.\n') )
    
    res_list$hyperprior <- res
    
  }
  
  res_list
  
}

### Create graphs and save results
output_function <- function( res, dat, suffix_file, suffix_graph, z ) {
  
  # Save results
  complete_res <- list( res = res, dat = dat )
  save( complete_res, file = 'res.Rdata' )
  rm( complete_res )
  
  # Make graphs
  x <- dat$x
  t <- dat$t
  
  # Plot the data
  png( paste0('Author_and_Discovery_Data_', z, '.png') )
  par( mfrow = c(2, 1) )
  plot( names(x), x, type = 'l', main = 'Number of Authors', xlab = 'Year', ylab = '' )
  plot( names(x), t, type = 'l', main = 'Discoveries per Year', xlab = 'Year', ylab = '' )
  par( mfrow = c(1, 1) )
  dev.off()
  
  # Plot the approximate posterior of C
  png( paste0('Total_Species_Num_Approx_Posterior_', z, suffix_file, '.png') )
  C_dens <- density( res$C )
  approx_mode <- C_dens$x[which.max(C_dens$y)]
  remove_inds_raw <- which( C_dens$x < sum(dat$t) )
  remove_inds <- remove_inds_raw[-length(remove_inds_raw)]
  C_dens$x <- C_dens$x[-remove_inds]
  C_dens$y <- C_dens$y[-remove_inds]
  C_dens$y[1] <- 0
  plot( C_dens, main = paste0(z, suffix_graph),
        xlab = '', ylab = '', xlim = c(0, max(C_dens$x)), yaxt = 'n' )
  #abline( v = approx_mode, lty = 2, col = 'blue' )
  dev.off()
  
  # Plot the simulated number of authors per year
  png( paste0('Effort_Proxy_ABC_', z, suffix_file, '.png') )
  xhat <- res$x_hat
  plot( x, type = 'n', ylim = range(rbind(x, xhat)), ylab = 'Effort Proxy', 
        xlab = 'Time', main = paste0(z, suffix_graph) )
  grid()
  for ( i in 1:nrow(xhat) ) {
    lines( xhat[i, ], col = alpha('grey', 0.5) )
  }
  lines( x )
  dev.off()
  
  # Plot the simulated number of discoveries per year
  png( paste0('Discovery_Time_ABC_', z, suffix_file, '.png') )
  that <- res$t_hat
  plot( t, type = 'n', ylim = range(rbind(t, that)), ylab = 'Number of Discoveries', 
        xlab = 'Time', main = paste0(z, suffix_graph) )
  grid()
  for ( i in 1:nrow(that) ) {
    lines( that[i, ], col = alpha('grey', 0.5) )
  }
  lines( t )
  dev.off()
  
}

analysis_function <- function( z, synonym_list, param_list, metric = 'm', p = 2,
                               pilot_iterations, fit_iterations, fixed_run = TRUE,
                               hyperprior_run = FALSE, x_weight = 1, t_weight = 1,
                               n_cores = detectCores() - 1, wd = getwd() ) {
  
  taxa_data <- read_taxa( paste0(z, '_20240815.zip') )
  
  ABC_data <- construct_ABC_data( taxa_data, synonym_list )
  
  res <- run_ABC( dat = ABC_data, param_list = param_list,
                  metric = metric, p = p,
                  pilot_iterations = pilot_iterations, 
                  fit_iterations = fit_iterations, 
                  fixed_run = fixed_run, hyperprior_run = hyperprior_run, 
                  x_weight = x_weight, t_weight = t_weight, n_cores = n_cores )
  
  default_wd <- getwd( )
  if ( !dir.exists(wd) ) { dir.create( wd ) }
  setwd( wd )
  
  if ( fixed_run ) {
    output_function( res$fixed, ABC_data, '', '', z )
  }
  
  if ( hyperprior_run ) {
    output_function( res$hyperprior, ABC_data, '_hyperprior', ' Hyperprior', z )
  }
  
  setwd( default_wd )
  
  return( invisible(NULL) )
    
}

