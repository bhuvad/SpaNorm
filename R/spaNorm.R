#' spatially-dependent normalisation for spatial transcriptomics datas
#'
#' Performs normalisation of spatial transcriptomics data using spatially-dependent spot-specific size factor.
#'
#' @param spe SpatialExperiment object with the raw data stored in assays(spe)$counts
#' @param pCells.touse the (maximum) proportion of cells used for model fitting. The default is 0.25.
#' @param lambda.a smoothing parameter for regularizing regression coefficients. The default is 0.0001*ncol(spe)
#' @param step.fac multiplicative factor to decrease IRLS step by when log-likelihood diverges. Default is 0.5
#' @param inner.maxit the maximum number of IRLS iteration for estimating mean parameters for a given dispersion parameter. Default is 50
#' @param outer.maxit the maximum number of IRLS iteration. Default is 25 
  
#' @return SpatialExperiment object with the adjusted data stored in assays(spe)$logPAC
#' @export

spaNorm <- function(spe,pCells.touse=0.25,lambda.a=0.0001,step.fac=0.5,inner.maxit=50,outer.maxit=25) {
 lambda.a <- lambda.a*ncol(spe)
 cl <- scran::quickCluster(spe)
 spe <- scran::computeSumFactors(spe, clusters = cl)
 spe$sizeFactor <- pmax(1e-08,spe$sizeFactor)
 spe <- subset(spe,matrixStats::rowMeans2(as.matrix(assays(spe)$counts))>0.1)
 spe$logLS <- log(spe$sizeFactor)
 x <- spatialCoords(spe)[,1]
 y <- spatialCoords(spe)[,2]
 bs.x <- splines::ns(x,df=6)
 bs.y <- splines::ns(y,df=6)
 bs.xy<- matrix(0,nrow=nrow(bs.x),ncol=ncol(bs.x)*ncol(bs.y))
 for(i in 1:ncol(bs.x)) {
  for(j in 1:ncol(bs.y)) {
    bs.xy[,(i-1)*ncol(bs.x)+j] <- bs.x[,i]*bs.y[,j]
  }
 }
 bs.xy <- scale(bs.xy,scale=F)

 Y   <- assays(spe)$counts
 idx <- sample(ncol(Y),size=min(3000,round(pCells.touse*ncol(Y))))
 Ysub<- as.matrix(Y[,idx]) 
 nsub<- ncol(Ysub)
 W   <- model.matrix(~spe$logLS*bs.xy)[,-1]
 nW  <- ncol(W)
 Wsub<- W[idx,]

 # initial estimate
 lmu.hat<- Z <- log(Ysub+1) 
 gmean  <- rowMeans(Z)
 alpha  <- matrix(0,nrow(Ysub),ncol(W))
 # alpha for logLS
 alpha[,1] <- 1

 # estimate initial psi
 psi.idx <- sample(ncol(Ysub),size=min(50,round(0.1*ncol(Y))))
 bef <- Sys.time() 
 offs   <- Wa <- Matrix::tcrossprod(alpha,Wsub) 
 offs.psi <- offs[,psi.idx]
 Ysub.psi <- Ysub[,psi.idx]
 nsub.psi  <- ncol(Ysub.psi)
 psi <- edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion

# initiate
conv <- FALSE
iter.outer <- 0
logl.outer <- NULL
step <- rep(1,nsub)
print('Start...')
while(!conv) {
 iter.outer <- iter.outer + 1
 logl.beta <- NULL
 
 ### calculate initial loglik
 lmu.hat    <- gmean + Wa
 sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))

 temp <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
 loglik <- colSums(temp)

 #############################
 step <- rep(1,nsub)
 conv.beta <- FALSE
 # saving temporary best estimate
 best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.gmean <- gmean
 
 iter <- 1
 halving <- 0
 while(!conv.beta) {
 #print(paste0('Inner iter:', iter))
 # calc working vector Z for a subset of cells 
 lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 

 sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))
 Z      <- lmu.hat + sweep( ((Ysub+0.01)/(exp(lmu.hat)+0.01) - 1),2,step,'*')

 # store current alpha est
 alpha.old <- alpha

 # get alpha 
 bef <- Sys.time()
 wt.cell  <- matrixStats::colMeans2(sig.inv)
 # prevent outlier large weight
 p98      <- quantile(wt.cell,probs=0.98)
 wt.cell[which(wt.cell>=p98)] <- p98
 b         <- (Z-gmean) %*% diag(wt.cell) %*% Wsub
 alpha<- b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub*wt.cell,Wsub)))
 
#print(dim(alpha))
#print(dim(b))
#print(dim(Wsub))
#print(length(wt.cell))

 alpha[,1] <- mean(alpha[,1])
 b         <- (Z-gmean-Matrix::tcrossprod(alpha[,1,drop=FALSE],Wsub[,1,drop=FALSE])) %*% diag(wt.cell) %*% Wsub[,-1]
 alpha[,-1]<- b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub[,-1]*wt.cell,Wsub[,-1])+lambda.a*diag(nW-1)))

 aft <- Sys.time()

 # if new alpha est has missing/inf values, revert to previous estimate
 if(any(is.na(alpha) | is.infinite(alpha) )) alpha <- alpha.old

 # reduce outliers
 a.med <- matrixStats::colMedians(alpha)
 a.mad <- matrixStats::colMads(alpha)
 for(i in 1:nW) {
   alpha[which(alpha[,i]> (a.med[i]+4*a.mad[i])),i] <- a.med[i]+4*a.mad[i]
   alpha[which(alpha[,i]< (a.med[i]-4*a.mad[i])),i] <- a.med[i]-4*a.mad[i]
 } 

 #step2: update gmean
 Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  
 bef <- Sys.time()
 gmean<- matrixStats::rowSums2(Z.res*sig.inv)/matrixStats::rowSums2(sig.inv)
 aft <- Sys.time()

 # calculate current logl
 lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 
 temp <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
 loglik.tmp <- colSums(temp) 

 # check degenerate case
 if(iter>=1) {
  degener<- sum(loglik)>sum(loglik.tmp)
  degener[is.na(degener)] <- TRUE
  if(degener) {
    check.gene <- loglik>loglik.tmp
    step[check.gene] <- step[check.gene]*step.fac
    Wsub <- best.W; alpha <- best.a ; psi <- best.psi ; gmean <- best.gmean 
    halving <- halving + 1
    if(halving>=3) {
     loglik.tmp <- loglik
     degener <- !degener
    }
  }
  if(!degener) {
   loglik     <- loglik.tmp
   logl.beta  <- c(logl.beta,sum(loglik))
   #print(paste0('Outer Iter ',iter.outer, ', Inner iter ', iter, ' logl-likelihood:', logl.beta[length(logl.beta)]))
   # saving temporary best estimate
   best.W <- Wsub; best.psi <- psi ; best.a <- alpha; best.gmean <- gmean 
   iter <- iter + 1
   halving <- 0
  }
 }

conv.logl <- FALSE
if(iter>2) {
 conv.logl <- ifelse( (logl.beta[length(logl.beta)] - logl.beta[length(logl.beta)-1])/abs(logl.beta[length(logl.beta)-1]) < 1e-04,TRUE,FALSE)
}
conv.beta <- iter>=inner.maxit  | conv.logl
} # end of IRLS inner loop

#update psi
logl.outer.tmp <- sum(loglik,na.rm=T)
updt.psi <- TRUE
if(iter.outer>1)
   updt.psi <- abs(logl.outer[length(logl.outer)]-logl.outer.tmp)/abs(logl.outer[length(logl.outer)]) > 1e-08  

Wa <- offs <- Matrix::tcrossprod(best.a,best.W) 
offs.psi <- offs[,psi.idx]
if(updt.psi) 
  psi.new <- tryCatch({ edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {NA})
# update psi
psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]


# record outer logl
Wsub <- best.W
alpha <- best.a
gmean <- best.gmean
# new changes
best.psi <- psi

# recalculate logl
lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 
temp       <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
print(paste0('Outer Iter ',iter.outer, ' logl-likelihood:', logl.outer[length(logl.outer)]))

if(iter.outer>1) {
  conv <- ifelse( (logl.outer[iter.outer] - logl.outer[iter.outer-1])/abs(logl.outer[iter.outer-1]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
}
}
# calculate residual
mu      <- exp( gmean + Matrix::tcrossprod(alpha,W))
pred2.mu<- mean(exp(spe$logLS)) * exp(gmean + Matrix::tcrossprod(alpha[,2:37,drop=FALSE],W[,2:37,drop=FALSE]))
#logPAC
lb <- pnbinom(as.matrix(Y)-1,mu=mu, size= 1/psi)
p  <- 0.5*lb + 0.5*dnbinom(as.matrix(Y),mu=mu, size=1/psi)
p[which(p>0.99)] <- 0.99
p[which(p<0.01)] <- 0.01
PAC <- qnbinom(p,mu=pred2.mu,size=1/psi)
assays(spe,withDimnames=FALSE)$logPAC <- log(PAC+1)
return(spe)
}

