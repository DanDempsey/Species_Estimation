##### Fit Truncated Data and Extrapolate
##### Daniel Dempsey

### Load in analysis code
source( 'Code/Application/0_Data_Analysis_Function.R' )

### RNG settings
RNGkind( "L'Ecuyer-CMRG" ) # Carries RNG stream across parallel computing, so that parallel code is reproducible

### Utility function for computing discovery times
t_val <- function( d, n_sum ) {
  sum( n_sum < d ) + 1
}

### Extrapolation function
extrap <- function( C, D, x, theta_x, theta_l, theta_N, tstart, tend ) {
  
  Tau <- tend - tstart + 1
  
  # Latent effort L
  L <- numeric( Tau )
  L[1] <- log( x ) + theta_x
  for ( i in 2:length(L) ) {
    L[i] <- rnorm( 1, L[i-1], theta_l[2] )
  }
  
  # Class abundance N
  N <- rlnorm( C, theta_N[1], theta_N[2] )
  N_inds <- sample( 1:C, C, prob = N )
  N <- N[N_inds]
  N_sum_prop <- cumsum( N[1:(C-1)] ) / sum( N )
  
  # Inter-discovery times d
  d <- rep( 1, C )
  d[2:C] <- rgeom( n = C-1, 1 - N_sum_prop )
  d_sum <- cumsum( d )#[-(1:D)]
  
  # Sample size n
  n <- rep( d_sum[D], Tau )
  n[2:Tau] <- rpois( Tau - 1, exp(L[-1]) )
  n_sum <- cumsum( n )
  
  # Simulate t
  t <- sapply( d_sum, t_val, n_sum )
  t_sim <- as.numeric( table(t[t <= Tau] %>% factor(levels = 1:Tau)) )
  
  # Return simulated data
  t_sim[-1]
  
}

### Run truncated data
classes <- c( 'Aves', 'Isopoda', 'Mammalia', 'Nematoda' )
n_cores <- detectCores() - 1
set.seed( 1234 )
for ( i in 1:length(classes) ) {
  
  # Load in parameters
  par_list <- list()
  load( paste0('Output/Real_Data_Analysis/', classes[i], '/param_search.Rdata') )
  opt_row <- res[ which.min(res$distance), ]
  par_list$theta_N <- c( opt_row$mu_N, opt_row$sigma_N )
  par_list$theta_l <- c( opt_row$mu_l, opt_row$sigma_l )
  par_list$theta_x <- opt_row$theta_x
  par_list$m <- 50
  
  # Load in data
  taxa_data <- read_taxa( paste0(classes[i], '_20240815.zip') )
  ABC_data <- construct_ABC_data( taxa_data )
  train_inds <- names( ABC_data$x ) %>% as.numeric %>% '<='( 1920 )
  x_dat <- list( x_train = ABC_data$x[train_inds], x_test = ABC_data$x[!train_inds] )
  t_dat <- list( t_train = ABC_data$t[train_inds], t_test = ABC_data$t[!train_inds] )
  train_dat <- list( x = x_dat$x_train, t = t_dat$t_train )
  
  # Run code
  res1920 <- run_ABC( dat = train_dat, param_list = par_list,
                      pilot_iterations = 1e4, 
                      fit_iterations = 1e6, 
                      n_cores = n_cores )
  save( res, file = paste0('Output/Real_Data_Analysis/', classes[i], '/res_1920.Rdata') )
  
  # Plot density with original density
  load( file = paste0('Output/Real_Data_Analysis/', classes[i], '/res.Rdata') )
  dens_1920 <- density( res1920$fixed$C )
  dens_C <- density( complete_res$res$C )
  xrange <- c( dens_1920$x, dens_C$x ) %>% range
  yrange <- c( dens_1920$y, dens_C$y ) %>% range
  
  png( paste0('Output/Real_Data_Analysis/', classes[i], '/density_overlap.png') )
  plot( dens_1920, xlim = xrange, ylim = yrange, main = classes[i],
        xlab = 'C', ylab = '', yaxt = 'n', col = 'dodgerblue' )
  lines( dens_C, col = 'orange' )
  legend( 'topright', c('1920 data', 'Full data'), col = c('dodgerblue', 'orange'), 
          lty = 1, bty = 'n' )
  dev.off()
  
  # Perform extrapolation
  D <- sum(train_dat$t)
  C_boot <- sample( res1920$fixed$C, 1000, replace = TRUE )
  t_forecast <- mclapply( C_boot, extrap, D = D, x = train_dat$x[length(train_dat$x)], 
                          theta_x = par_list$theta_x, theta_l = par_list$theta_l, theta_N = par_list$theta_N,
                          tstart = length(train_dat$x) + 1, tend = length(ABC_data$x), mc.cores = n_cores )
  
  # Plot forecasts
  ind_1920 <- length(t_dat$t_train) + 1
  ind_2019 <- length(ABC_data$t)
  all_curves <- lapply( t_forecast, function(x) { cumsum(c(D, x)) } )
  yrange <- range( cumsum(c(t_dat$t_train, t_dat$t_test)), all_curves )
  
  png( paste0('Output/Real_Data_Analysis/', classes[i], '/model_assessment_', classes[i], '.png') )
  plot( 0, type = 'n', xaxt = 'n', xlab = '', xlim = c(0, ind_2019),
        ylab = 'Cumulative Discoveries', ylim = yrange, main = classes[i] )
  axis( 1, at = c(ind_1920 - 1, ind_2019), label = c(1920, 2019) )
  lines( 2:ind_1920 - 1, cumsum(t_dat$t_train) )
  for ( j in 1:length(all_curves) ) {
    lines( ind_1920:ind_2019, all_curves[[j]], col = 'grey' )
  }
  lines( ind_1920:ind_2019, cumsum(t_dat$t_test) + D, lty = 2 )
  dev.off()
  
}

