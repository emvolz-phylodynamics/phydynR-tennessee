# author: Erik Volz
# skyspline package available at https://github.com/emvolz-phylodynamics/skyspline
require(skyspline)

TREEFN <- 'vbdata/vccc_dater.tre'
tr <- read.tree(TREEFN)
tips <- tr$tip
SOURCE <- substr( tips , 1, 1) != '9'
btre <- drop.tip( tr, tips[!SOURCE] )

yrs <- sapply( strsplit( btre$tip, '_'), function(x) as.numeric( tail(x,1)))
mst <- max( yrs , na.rm=T) 
ndepth <- node.depth.edgelength( btre )
maxheight <- max( ndepth )
n <- length( btre$tip )
sts <- setNames( mst - maxheight + ndepth[1:n], btre$tip )
bdt <- DatedTree( btre, sts , tol = .1)

ss1 <- skyspline( bdt
 , death_rate_guess = 1/10
 , R0guess = 5
 , t0 = 1980 #mst - maxheight
 , np_range = 3, hessian=F, trace=T)

##
ss <- ss1 
with (ss$demo.history, plot( times, reproduction.number * GAMMA / pop.size  ) ) 

saveRDS( ss, 'srcGrowth0.rds') 
