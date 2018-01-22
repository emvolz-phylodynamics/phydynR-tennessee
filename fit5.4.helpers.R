# author : Erik Volz
# date : August 1, 2017
# Loads the model and data for downstream ML or Bayesian analysis 
# makes starting conditions 
# defines parameter names and transforms 

source('m5.4.R') 
require(parallel)

NCLADES <- 4
MAXHEIGHT <- 30
mc_cores <- 4

## load tres
tres <- readRDS( 'vccc_dater_clades.rds' )

## load states of sampled lineages (age, stage etc )
X <- readRDS('mstates0.rds' )[[1]]
colnames( X) <- DEMES

# make dated trees 
.tre2bdt <- function(tre )
{
	tip2yrs <- setNames( sapply( strsplit( tre$tip, '_'), function(x) as.numeric( tail(x,1))), tre$tip) 
	mst <- max( tip2yrs[tre$tip]  , na.rm=T)  + 1
	n <- length( tre$tip )
	sts <- node.depth.edgelength( tre )[1:n] 
	sts <- setNames( sts - max(sts ) + mst , tre$tip )
	bdt <- DatedTree( tre, sts , sampleStates = X[tre$tip, ] , minEdgeLength = 1/12, tol = 1.5)
	bdt
}

bdts <-  lapply( tres[1:NCLADES], .tre2bdt ) 

# parameter names (shared by all clades)
pnames_pooled <- c( 'wchron', 'w2', 'importRate', 'p12', 'p22') 
log_vars_pooled <- c( 'wchron', 'w2', 'importRate')
logistic_vars_pooled <- c('p12', 'p22')

# parameter names for parameters estimated in each clade 
pnames_clade <- c('beta', 'N' , 'I0')
np_pooled <- length( pnames_pooled )

pnames0 <- c(pnames_clade, pnames_pooled )

##  make start conditions based on latin-hypercubes
require(lhs )
NSTART <- 20 
lhsmat <- improvedLHS( length(pnames0) , NSTART)
pbounds <- cbind( 
  c(beta =.04, N = 600, I0 = .1, wchron = .02, w2 = .5, importRate = 1/40, p12 = .01, p22 = .2)
  , c(beta =.35, N = 10e3, I0 = 5, wchron = 2, w2 = 50, importRate = 1/5, p12 = .5, p22 = .99)
)
starts <- ( pbounds[, 1] + (pbounds[,2] - pbounds[,1] ) * ( lhsmat ) ) # each col is a start 
rownames( starts ) <- pnames0

## box constraint 
.parm.limits <- function(theta)
{ # require parms to lie within bounds used for start conditions 
	pnames <- intersect( names(theta), rownames( pbounds))
	all( (theta[pnames] >= pbounds[pnames,1]) & (theta[pnames] <= pbounds[pnames,2]) )
}


# transform parameters to log or logistic scale 
nat2trscale <- function( .theta){
	.theta0 <- .theta
	pn <- names(.theta0)
	.theta0[ intersect(pn,  log_vars_pooled) ] <- log(  .theta[ intersect(pn,  log_vars_pooled) ] )
	.theta0[ intersect(pn,  pnames_clade) ] <- log(  .theta[ intersect(pn,  pnames_clade) ] )
	
	x <-  .theta[ intersect(pn,  logistic_vars_pooled) ]
	.theta0[ intersect(pn,  logistic_vars_pooled) ] <- log( x/ ( 1- x ) )
	.theta0
}
# transform back to natural scale 
tr2natscale <- function( trtheta ){
	.theta <- exp( trtheta ) 
	x <- .theta[ logistic_vars_pooled ]
	.theta[ logistic_vars_pooled ] <- x / (1 + x )
	.theta
}
