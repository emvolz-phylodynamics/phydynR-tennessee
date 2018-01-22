# author : Erik Volz
# date : August 1 2017
# Bayesian MCMC fit model 
# Metropolis-hastings

# iterations 
ITER <- 10e3 

require(akima)
source('fit5.4.helpers.R') # load the model and data 


## PRIORS
ln.dprior0 <- function(theta, theta_clade)
{
	lnd<- 0
	for (n in names(theta)){
		if (n=='p12') {
			lnd <- lnd + dbeta( theta[n], 2, 2, log=T)
		}
		if (n=='p22') {
			lnd <- lnd + dbeta( theta[n], 2, 2, log=T)
		}
		if (n=='importRate') {
			lnd <- lnd + dexp( theta[n], 1/20, log=T)
		}
		if (n=='w2') {
			lnd <- lnd + dlnorm( theta[n], log(1), 1, log=T)
		}
		if (n=='wchron') {
			lnd <- lnd + dlnorm( theta[n], log(1), 2, log=T)
		}
	}
	
	for (n in rownames(theta_clade ) ){
		if (n == 'I0'){
			lnd <- lnd + sum( dexp( theta_clade[n,], 1, log=T) )
		}
		if (n == 'N'){
			lnd <- lnd + sum( dlnorm( theta_clade[n,], log(1e3), 2, log=T) )
		}
		if (n == 'beta'){
			lnd <- lnd + sum( dlnorm( theta_clade[n,], log(.1), 1, log=T) )
		}
	}
	lnd
}

## PROPOSALS
# parameters estimated for each clade: 
prop_sd_clade <- c( beta = .0025, N = 500, I0 = .1)
# parameters shared by all clades:
pnames_pooled <- c( 'wchron', 'w2', 'importRate', 'p12', 'p22') 

prop_sd_pooled <- c( wchron = .1
	, w2 = .5
	, importRate = .01
	, p12 = .01
	, p22 = .05
)

prop0 <- function(theta, theta_clade, iter, focus_clade = 5e3 )
{
	vn_clade <- rownames(theta_clade )
	np_clade <- length( vn_clade )
	
	.theta <- theta
	.theta_clade <- theta_clade
	if ( (iter %% 2 == 1) & (iter > focus_clade)){ 
		l <- 1 + (iter %% np_pooled )
		n <- pnames_pooled[ l] 
		.theta[n] <- rnorm( 1, .theta[n], prop_sd_pooled[n] )
	} else{ # update a clade var in order 
		k <- 1 + floor( iter/2 )
		l <- 1 + (k %% np_clade )
		n <- vn_clade[l] 
		.theta_clade[n,] <- theta_clade[n, ] + rnorm( ncol(theta_clade), 0, prop_sd_clade[ n ] )
	}
	list(.theta ,  .theta_clade )
}


## OBJFUN 
of.mcmc0 <- function( xtheta, xtheta_clade){
	lndp <- ln.dprior0( xtheta, xtheta_clade )
	if (is.infinite( lndp)) return(lndp)
	
	.theta<- theta
	.theta[names(xtheta)] <- unname( xtheta  )
	.theta[ c('wb', 'wc')] <- unname( xtheta['wchron'] )
	
	ll <- sum( unlist(  mclapply( 1:NCLADES, function(k) {
		bdt <- bdts[[k]]
		
		.theta[pnames_clade] <- unname(  xtheta_clade[pnames_clade,k] ) # set beta I0 and N for this clade
		if ( !.parm.limits( .theta)) return(-Inf)
		
		.parms <- as.list( .theta ) 
		.parms$ssbeta.fun <- ssbeta.fun
		
		.x0 <- x0 
		.x0[-m] <- .x0[-m] * .theta['I0'] 
		suppressWarnings( { ll <- colik.pik(bdt, .parms, dm, .x0, 1980 , maxHeight=MAXHEIGHT
			,  integrationMethod = 'lsoda'
			, res = 12*MAXHEIGHT) })
		print(c( date(), ll ) )
		ll
	}, mc.cores = mc_cores ) ) )
print('#####' )
print( xtheta  )
print( xtheta_clade  )
print(c( date(), ll ) )
print('#####' )
	if (is.na( ll)) return(-Inf)
	
	lnd <- ll + lndp
}

## MCMC 
MAXHEIGHT <- 15

LOGPO <- rep( NA, ITER)
THETA_POOLED <- matrix( NA, nrow=ITER, ncol = length(pnames_pooled ))
THETA_CLADE <- array( NA, dim = c( length(pnames_clade), NCLADES, ITER) )

.LOGPO <- rep( NA, ITER)
.THETA_POOLED <- matrix( NA, nrow=ITER, ncol = length(pnames_pooled ))
.THETA_CLADE <- array( NA, dim = c( length(pnames_clade), NCLADES, ITER) )



## START CONDITIONS
i_rand_start <- sample(1:ncol(starts), size=NCLADES)
theta_pooled <- rowMeans( starts[ pnames_pooled, i_rand_start ] )
theta_clade <- starts[ pnames_clade, i_rand_start ]
logpo <- of.mcmc0( theta_pooled, theta_clade )



for (iter in 1:ITER){
	cat('*****************\n')
	print( paste( "ITER" , iter ))
	
	prop <- prop0( theta_pooled, theta_clade, iter )
	.theta_pooled <- prop[[1]]
	.theta_clade <- prop[[2]]
	.logpo <- of.mcmc0( .theta_pooled, .theta_clade )
	
	.LOGPO[iter] <- .logpo
	.THETA_POOLED[iter,] <- .theta_pooled
	.THETA_CLADE[,,iter] <- .theta_clade
	
	if ( runif(1) < exp( .logpo - logpo ) ){
		logpo <- .logpo
		theta_pooled <- .theta_pooled
		theta_clade <- .theta_clade
	} 
	LOGPO[iter] <- logpo
	THETA_POOLED[iter,] <- theta_pooled
	THETA_CLADE[,,iter] <- theta_clade
}


## SAVE 
ofn <- paste0( 'f5.4.mcmc0/', Sys.getpid(), '.rds' )
saveRDS( list( logpo = LOGPO, clade = THETA_CLADE, pooled = THETA_POOLED ) , ofn)

