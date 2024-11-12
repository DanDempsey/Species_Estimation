##### Run Data Analyses 
##### Daniel Dempsey

### Load in analysis code
source( 'Code/Application/0_Data_Analysis_Function.R' )

### RNG settings
RNGkind( "L'Ecuyer-CMRG" ) # Carries RNG stream across parallel computing, so that parallel code is reproducible

### Set parameter values
classes <- c( 'Aves', 'Isopoda', 'Mammalia', 'Nematoda' )
par_list <- list()
for ( i in 1:length(classes) ) {
  
  par_list[[i]] <- list()
  load( paste0('Output/Real_Data_Analysis/', classes[i], '/param_search.Rdata') )
  opt_row <- res[ which.min(res$distance), ]
  par_list[[i]]$theta_N <- c( opt_row$mu_N, opt_row$sigma_N )
  par_list[[i]]$theta_l <- c( opt_row$mu_l, opt_row$sigma_l )
  par_list[[i]]$theta_x <- opt_row$theta_x
  par_list[[i]]$m <- 50
  
}
names( par_list ) <- classes

### Run fit for all data
pits <- 1e4
fits <- 1e6

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

# Discovery curves
png( 'Discovery_Curves.png' )
discovery_curves <- lapply( all_res, function(x) { cumsum(c(0, x$dat$t)) } )
years <- 1757:2019
nms <- names( discovery_curves )
plot( years, discovery_curves[[1]], type = 'n', ylim = range(discovery_curves),
      xlab = 'Year', ylab = '', main = 'Cumulative Number of Discoveries')
grid()
lines( years, discovery_curves[[1]], col = 2 )
lines( years, discovery_curves[[2]], col = 3 )
lines( years, discovery_curves[[3]], col = 4 )
lines( years, discovery_curves[[4]], col = 7 )
text( 2010, 11900, nms[1], col = 2 )
text( 2010, 7500, nms[2], col = 3 )
text( 2010, 4500, nms[3], col = 4 )
text( 1990, 15500, nms[4], col = 7 )
dev.off()

# Approximate posterior of C
png( 'Total_Species_Num_Approx_Posterior.png', height = 600, width = 600 )
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
  abline( v = sum(dat$t), lty = 2 )
  
}
dev.off()

# Accepted simulated effort proxies
png( 'Effort_Proxy_ABC.png', height = 600, width = 600 )
par( mfrow = c(2, 2), mar = c(3, 2, 3, 2) + 0.1 )
for ( i in 1:length(par_list) ) {
  
  res <- all_res[[i]]$res
  dat <- all_res[[i]]$dat
  nm <- names( all_res )[i]
  
  x <- dat$x
  t <- dat$t
  
  xhat <- res$x_hat
  plot( years[-1], x, type = 'n', ylim = range(rbind(x, xhat)), 
        ylab = 'Effort Proxy', xlab = 'Time', main = nm )
  grid()
  for ( j in 1:nrow(xhat) ) {
    lines( years[-1], xhat[j, ], col = alpha('grey', 0.5) )
  }
  lines( years[-1], x )
  
}
dev.off()

# Accepted simulated discovery times
png( 'Discovery_Time_ABC.png', height = 600, width = 600 )
par( mfrow = c(2, 2), mar = c(3, 2, 3, 2) + 0.1 )
for ( i in 1:length(par_list) ) {
  
  res <- all_res[[i]]$res
  dat <- all_res[[i]]$dat
  nm <- names( all_res )[i]
  
  x <- dat$x
  t <- dat$t
  
  that <- res$t_hat
  plot( years[-1], t, type = 'n', ylim = range(rbind(t, that)), 
        ylab = 'Number of Discoveries', xlab = 'Time', main = nm )
  grid()
  for ( j in 1:nrow(that) ) {
    lines( years[-1], that[j, ], col = alpha('grey', 0.5) )
  }
  lines( years[-1], t )
}
dev.off()

# Summary statistics
C_list <- lapply( all_res, function(x) { x$res$C } )
summary_res <- as.data.frame( matrix(0, nrow = 4, ncol = 4) )
colnames( summary_res ) <- c( 'Mean', 'Mode', 'Median', '95% CI' )
rownames( summary_res ) <- classes

summary_res$Mean <- sapply( C_list, mean ) %>% round(-3)
summary_res$Mode <- sapply( C_list, function(x) { z <- density(x); z$x[which.max(z$y)] } ) %>% round(-3)
summary_res$Median <- sapply( C_list, quantile, probs = 0.5 ) %>% round(-3)
summary_res$`95% CI` <- sapply( C_list, function(x) { 
  paste0(round(quantile(x, probs = c(0.025, 0.975)), -3), collapse = ' - ') } )

xtable::xtable( summary_res )

setwd( '../../..' )

