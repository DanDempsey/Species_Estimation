##### Run Data Analyses 
##### Daniel Dempsey

### Load in analysis code
source( 'Code/Application/0_Data_Analysis_Function.R' )

### RNG settings
RNGkind( "L'Ecuyer-CMRG" ) # Carries RNG stream across parallel computing, so that parallel code is reproducible

### Parameter values
m <- 20

# Aves
Aves_parameters <- list()
Aves_parameters$theta_N <- c( log(1e9), 1.1 ) 
Aves_parameters$theta_l <- c( log(30), 0.08 )
Aves_parameters$theta_x <- log( 2 )
Aves_parameters$m <- m

# Isopoda
Isopoda_parameters <- list()
Isopoda_parameters$theta_N <- c( log(1e9), 1.5 )
Isopoda_parameters$theta_l <- c( log(0.0001), 1 )
Isopoda_parameters$theta_x <- log( 0.6 )
Isopoda_parameters$m <- m

# Mammalia
Mammalia_parameters <- list()
Mammalia_parameters$theta_N <- c( log(1e9), 1.1 ) 
Mammalia_parameters$theta_l <- c( log(30), 0.15 )
Mammalia_parameters$theta_x <- log( 1 )
Mammalia_parameters$m <- m

# Nematoda
Nematoda_parameters <- list()
Nematoda_parameters$theta_N <- c( log(1e9), 1.5 )
Nematoda_parameters$theta_l <- c( log(0.0001), 1 )
Nematoda_parameters$theta_x <- log( 1.2 )
Nematoda_parameters$m <- m

# Combine all
par_list <- list( Aves = Aves_parameters, Isopoda = Isopoda_parameters, 
                  Mammalia = Mammalia_parameters, Nematoda = Nematoda_parameters )

### Run fit for all data
pits <- 1e4
fits <- 1e4

set.seed( 20240820 )
for ( i in 1:length(par_list) ) {
  
  nm <- names( par_list )[i]
  cat( paste0(nm, ':\n') )
  analysis_function( nm, param_list = par_list[[i]], 
                     pilot_iterations = pits, fit_iterations = fits, 
                     wd = paste0('Output/Real_Data_Analysis/', nm) )
  cat( '\n' )
  
}

### Collate graphs
wd <- 'Output/Real_Data_Analysis/All'
if ( !dir.exists(wd) ) { dir.create( wd ) }
setwd( wd )

all_res <- list()
for ( i in 1:length(par_list) ) {
  
  nm <- names( par_list )[i]
  load( paste0('../', nm, '/res.Rdata') )
  all_res[[nm]] <- complete_res
  
}

# Approximate posterior of C
png( 'Total_Species_Num_Approx_Posterior.png' )
par( mfrow = c(2, 2), mar = c(3, 1, 3, 2) + 0.1 )
for ( i in 1:length(par_list) ) {
  
  res <- all_res[[i]]$res
  dat <- all_res[[i]]$dat
  nm <- names( all_res )[i]
  
  x <- dat$x
  t <- dat$t
  
  C_dens <- density( res$C )
  approx_mode <- C_dens$x[which.max(C_dens$y)]
  remove_inds_raw <- which( C_dens$x < sum(dat$t) )
  remove_inds <- remove_inds_raw[-length(remove_inds_raw)]
  C_dens$x <- C_dens$x[-remove_inds]
  C_dens$y <- C_dens$y[-remove_inds]
  C_dens$y[1] <- 0
  plot( C_dens, main = nm, xlab = '', ylab = '', 
        xlim = c(0, max(C_dens$x)), yaxt = 'n' )
  #abline( v = approx_mode, lty = 2, col = 'blue' )
  
}
dev.off()

# Plot accepted simulated effort proxies
png( 'Effort_Proxy_ABC.png' )
par( mfrow = c(2, 2), mar = c(3, 2, 3, 2) + 0.1 )
for ( i in 1:length(par_list) ) {
  
  res <- all_res[[i]]$res
  dat <- all_res[[i]]$dat
  nm <- names( all_res )[i]
  
  x <- dat$x
  t <- dat$t
  
  xhat <- res$x_hat
  plot( x, type = 'n', ylim = range(rbind(x, xhat)), ylab = 'Effort Proxy', 
        xlab = 'Time', main = nm )
  grid()
  for ( j in 1:nrow(xhat) ) {
    lines( xhat[j, ], col = alpha('grey', 0.5) )
  }
  lines( x )
  
}
dev.off()

# Plot accepted simulated discovery times
png( 'Discovery_Time_ABC.png' )
par( mfrow = c(2, 2), mar = c(3, 2, 3, 2) + 0.1 )
for ( i in 1:length(par_list) ) {
  
  res <- all_res[[i]]$res
  dat <- all_res[[i]]$dat
  nm <- names( all_res )[i]
  
  x <- dat$x
  t <- dat$t
  
  that <- res$t_hat
  plot( t, type = 'n', ylim = range(rbind(t, that)), ylab = 'Number of Discoveries', 
        xlab = 'Time', main = nm )
  grid()
  for ( j in 1:nrow(that) ) {
    lines( that[j, ], col = alpha('grey', 0.5) )
  }
  lines( t )
}
dev.off()

setwd( '../../..' )

