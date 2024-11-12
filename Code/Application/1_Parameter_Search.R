##### Parameter Search
##### Daniel Dempsey

### Load in analysis code
library( dplyr )
source( 'Code/Application/0_Data_Analysis_Function.R' )

### RNG settings
RNGkind( "L'Ecuyer-CMRG" ) # Carries RNG stream across parallel computing, so that parallel code is reproducible

### Function used to find optimal parameters
search_fun <- function( x, param_grid, m = 50, weights = c(1, 1), n_cores = detectCores()-1, 
                        its = 1e4, seed = 1234 ) {
  
  # Load in and prepare data
  taxa_data <- read_taxa( paste0(x, '_20240815.zip') )
  ABC_data <- construct_ABC_data( taxa_data )
  
  # Infer parameters that we don't need to use a grid search for
  mu_N <- log( 1e9 )
  max_ind <- which.max( ABC_data$t[-(1:2)] ) + 2
  theta_x <- ( ABC_data$t[max_ind] / ABC_data$x[max_ind] ) %>% log
  
  # Run grid search
  res <- mutate( param_grid, mu_N = mu_N, theta_x = theta_x, distance = 0 )
  for ( i in 1:nrow(res) ) {
    
    cat( paste0(x, ': iteration ', i, '\n') )
    set.seed( seed ) # Stochastic consistency across runs
    theta_x <- res$theta_x[i]
    theta_N <- c( res$mu_N[i], res$sigma_N[i] )
    theta_l <- c( res$mu_l[i], res$sigma_l[i] )
    all_dist <- CD_ABC( x = ABC_data$x, t = ABC_data$t, theta_N = theta_N, 
                        theta_x = theta_x, theta_l = theta_l, m = m, iters = its, 
                        x_weight = weights[1], t_weight = weights[2], epsilon = Inf, 
                        cores = n_cores )$distance
    res$distance[i] <- sort( all_dist )[1:round(its/10)] %>% mean
    
  }
  
  # Return result
  wd <- paste0( 'Output/Real_Data_Analysis/', x )
  if ( !dir.exists(wd) ) { dir.create( wd ) }
  save( res, file = paste0(wd, '/param_search.Rdata') )
  res
  
}

### Set parameter grid
param_grid <- expand.grid( paste0(3, 'e', seq(-4, 1)) %>% as.numeric %>% log, 
                           seq(0.08, 0.16, length.out = 5), 
                           seq(1, 2, length.out = 5) )
names( param_grid ) <- c( 'mu_l', 'sigma_l', 'sigma_N' )

### Run parameter search for each class
classes <- c( 'Aves', 'Isopoda', 'Mammalia', 'Nematoda' ) 
tic()
for ( nm in classes ) {
  search_fun( nm, param_grid )
}

### Visualise results
use_cols <- c( 'dodgerblue', 'orange' )
for ( nm in classes ) {
  
  load( file = paste0('Output/Real_Data_Analysis/', nm, '/param_search.Rdata') )
  png( paste0('Output/Real_Data_Analysis/', nm, '/grid_search_results.png') )
  par( mfrow = c(2, 2) )
  
  min_ind <- which.min( res$distance )
  mu_l_min <- res$mu_l[min_ind]
  sigma_l_min <- res$sigma_l[min_ind]
  sigma_N_min <- res$sigma_N[min_ind]
  
  mu_l_col <- ifelse( sort(unique(res$mu_l)) == mu_l_min, 2, 1 )
  sigma_l_col <- ifelse( sort(unique(res$sigma_l)) == sigma_l_min, 2, 1 )
  sigma_N_col <- ifelse( sort(unique(res$sigma_N)) == sigma_N_min, 2, 1 )
  
  boxplot( distance ~ round(mu_l, 2), data = res, xlab = "mu_l", border = use_cols[mu_l_col], pch = 20 )
  boxplot( distance ~ sigma_l, data = res, border = use_cols[sigma_l_col], pch = 20 )
  boxplot( distance ~ sigma_N, data = res, border = use_cols[sigma_N_col], pch = 20 )
  dev.off()
  
}

