# author : Erik Volz
# date : August 1, 2017
# Model equations for interfacing with phydynR package 
# skyspline used for global reservoir with importation
# Includes assortativity parameters; p(recipient=young) is to be estimated for both old and young
# logistic growth, 
# 3 stages- a , b , c
# 2 age levels- 1 , 2

require(phydynR)
THRESHOLD_AGE <- 25.2 #

ssdata <- readRDS('srcGrowth0.rds')
GAMMA <- 1/10 # removal rate for SIR dynamics
ssbeta.fun <- with (ssdata$demo.history, approxfun( rev( times), (rev(reproduction.number)*GAMMA ), rule = 2 ) )

# MODEL

# stage progression parameters 
stage_prog_yrs <- c( .5 +3.32, 2.7, 5.50, 5.06 ) #cori AIDS
pstarts <- c( 
  pstartstage2 = 0.76
 , pstartstage3 = 0.19
 , pstartstage4 = 0.05
 , pstartstage5 = 0
)
time_in_stage <- unname( stage_prog_yrs* cumsum(pstarts) )


## parms

# parameter template : 
theta <- c(
 wb = .25 #rel to stage a
, wc = .25
, w2 = 2
, beta = .25 
, N = 5e3
, g1 = 1/time_in_stage[1] #  should refactor ga,gb, gc
, g2 = 1/sum(time_in_stage[2:3])
, g3 = 1/time_in_stage[4]
, p12 = .15 # prop recip in risk level 2
, p22 = .25 # prop recip in risk level 2
, mu21  = 1/ ( 26.40326  - 23.4661) # see mstates..R, rate young to old
, mu1 = 1/(78 - 40) # should depend on life expectancy and age in older group; see mstates..R
, importRate = 1/20 
, I0 = 1
, GAMMA = 1/10 # for source deme 
)

DEMES <- c( 'a1', 'b1', 'c1',  'a2', 'b2', 'c2' , 'src')
m <- length(DEMES)

## births (transmissions )
b <- matrix('0', nrow = m, ncol = m)
rownames(b) = colnames(b) <- DEMES 

b['a1', 'a1'] <- 'with(parms, (1-p12) * (a1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'
b['a1', 'a2'] <-  'with(parms, p12 * (a1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

b['b1', 'a1'] <- 'with(parms, (1-p12) * (wb*b1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'
b['b1', 'a2'] <-  'with(parms, p12 * (wb*b1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

b['c1', 'a1'] <- 'with(parms, (1-p12) * (wc*c1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'
b['c1', 'a2'] <-  'with(parms, p12 * (wc*c1/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

## risk 2 (young)
b['a2', 'a1'] <- 'with(parms, (1-p22) * (w2*a2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N)'
b['a2', 'a2'] <-  'with(parms, p22 * (w2*a2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

b['b2', 'a1'] <- 'with(parms, (1-p22) * (wb*w2*b2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'
b['b2', 'a2'] <-  'with(parms, p22 * (wb*w2*b2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

b['c2', 'a1'] <- 'with(parms, (1-p22) * (wc*w2*c2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'
b['c2', 'a2'] <-  'with( parms, p22 * (wc*w2*c2/((a1+w2*a2+wb*b1+wb*w2*b2+wc*c1+wc*w2*c2))) * (a1+a2+b1+b2+c1+c2)*beta*(N-((a1+a2+b1+b2+c1+c2))) / N )'

b['src', 'src'] <- 'parms$ssbeta.fun( t) * src'

## migrations (stage progression. aging and importation)
g <- matrix('0', nrow = m, ncol = m )
rownames(g) = colnames(g) <- DEMES 
# stage prog
g['a1', 'b1'] <- 'parms$g1 * a1'
g['a2', 'b2'] <- 'parms$g1 * a2'
g['b1', 'c1'] <- 'parms$g2 * b1'
g['b2', 'c2'] <- 'parms$g2 * b2'
# aging
g['a2', 'a1'] <- 'parms$mu21*a2' # young to old
g['b2', 'b1'] <- 'parms$mu21*b2'
g['c2', 'c1'] <- 'parms$mu21*c2'
# importation
for ( x in DEMES[-m]) g[m,x] <- paste0('parms$importRate * ' , x)
for ( x in DEMES[-m]) g[x,m] <- paste0('parms$importRate * ' , x)

## deaths
d <- setNames( rep('0', m ), DEMES)
d['a1'] = 'parms$mu1 * a1'
d['b1'] = 'parms$mu1 * b1'
d['c1'] = '(parms$mu1 + parms$g3) * c1' # young have approx no nat mort; including disease
d['c2'] = 'parms$g3 * c2' 
d['src'] = 'parms$GAMMA * src'


parms <- as.list(theta )
parms$ssbeta.fun <- ssbeta.fun
dm <- build.demographic.process( b, mig = g, deaths= d, parameterNames = names(parms), rcpp=F)

x0 <- c( c(a1 = 1, b1 = 1, c1 = 1, a2 = 1, b2 = 1, c2 = 1)  / (m-1)
 , src= unname( exp( ssdata$par['lny0'] )  )  # initial source size ; see srcGrowthRate
)
