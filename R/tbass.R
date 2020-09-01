#######################################################
# Authors: Andy Shen, Devin Francom, Los Alamos National Laboratory
# github.com/aashen12/TBASS
#######################################################

########################################################################
## main TBASS function
########################################################################

#' Includes auxiliary functions used in tbass() function.
#' pos() takes a vector x and returns max(0, x);
#' getd() uses Cholesky decomposition for TBASS likelihood to return V matrix and beta_hat;
#' makeBasis() creates a basis function given signs, variables, knot locations and data or X matrix.
#' @title t-Distributed Bayesian Adaptive Spline Surfaces (TBASS)
#' @name tbass
#' @description Robust BMARS function for Student's t-distributed likelihoods. See bass function in package BASS. Returns the estimated model parameters.
#' @usage tbass(X,y)
#' @param X a dataframe matrix of predictor values
#' @param y vector of response data, ideally with noise
#' @param max.int maximum degree of interaction, default is 3
#' @param max.basis maximum number of basis functions permitted by function, default is 50
#' @param tau2 prior for regression coefficients
#' @param nu degrees of freedom for Student's t distribution, default 10
#' @param nmcmc number of MCMC iterations, default 10000
#' @param g1 shape for IG prior on sigma^2, default 0
#' @param g2 scale for IG prior on sigma^2, default 0
#' @param h1 shape for gamma prior on lambda, default 10
#' @param h2 scale for gamma prior on lambda, default 10
#' @param verbose prints running output every ticker iterations if TRUE, can be changed
#' @return a list object with estimated model parameters
#' @export
#' @import mnormt

#Rcpp::sourceCpp("/Users/andyshen/Desktop/Git/TBASS/TBASS/getdC.cpp")

pos<-function(vec){
  (abs(vec)+vec)/2
}

makeBasis<-function(signs,vars,knots,datat,degree=1){
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1))
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2)
  }
}

getd<-function(X,v,s2,tau2,y){
  Vinv<-t(X)%*%diag(v)%*%X/s2+diag(ncol(X))/tau2
  Vinv.chol<-chol(Vinv)
  V<-chol2inv(Vinv.chol)#solve(Vinv)
  Vinv.ldet<-sum(log(diag(Vinv.chol)))#determinant(Vinv)$mod
  bhat<-V%*%t(X)%*%diag(v)%*%y/s2 #regression eqn
  d <- -.5*Vinv.ldet + .5*t(bhat)%*%Vinv%*%bhat
  return(list(d=d,bhat=bhat,Vinv.ldet=Vinv.ldet,V=V,Vinv.chol=Vinv.chol,Vinv=Vinv))
}

tbass <- function(X,y,max.int=3,max.basis=50,tau2=10^4,nu=10,nmcmc=10000,g1=0,g2=0,h1=10,h2=10,verbose=FALSE){
  ticker = nmcmc/10
  # see parameter verbose, should be a multiple of 100, 500, or 1000, default is 1000 for 10000 nmcmc iterations
  Xt<-t(X)
  n<-length(y)
  p<-ncol(X)
  ssy<-sum(y^2)
  knots<-signs<-vars<-array(dim=c(nmcmc,max.basis,max.int))
  nint<-matrix(nrow=nmcmc,ncol=max.basis) #J
  beta<-matrix(nrow=nmcmc,ncol=max.basis+1) # +1 for intercept
  s2<-lam<-nbasis<-rep(NA,nmcmc)
  v<-matrix(nrow=nmcmc,ncol=n)
  v[1,]<-1
  nbasis[1]<-0
  s2[1]<-1
  lam[1]<-1
  X.curr<-matrix(rep(1,n))

  d.curr<-getd(X.curr,v[1,],s2[1],tau2,y)

  count<-c(0,0,0) # count how many times we accept birth, death, change
  beta[1,1]<-d.curr$bhat

  for(i in 2:nmcmc){

    ## Reversible jump step

    move.type<-sample(c('birth','death','change'),1)
    if(nbasis[i-1]==0)
      move.type<-'birth'
    if(nbasis[i-1]==max.basis)
      move.type<-sample(c('death','change'),1)

    # set all of this iterations values to last iteration values...we'll change them if we accept a move below
    nbasis[i]<-nbasis[i-1]
    nint[i,]<-nint[i-1,]
    knots[i,,]<-knots[i-1,,]
    signs[i,,]<-signs[i-1,,]
    vars[i,,]<-vars[i-1,,]

    if(move.type=='birth'){
      nint.cand<-sample(max.int,1) # sample degree of interaction for new basis function
      knots.cand<-runif(nint.cand) # sample knots for new basis function
      signs.cand<-sample(c(-1,1),nint.cand,replace = T) # signs for new basis function
      vars.cand<-sample(p,nint.cand,replace = F) # variables to use in new basis function
      basis.cand<-makeBasis(signs.cand,vars.cand,knots.cand,Xt) # make the new basis function
      X.cand<-cbind(X.curr,basis.cand) # add the new basis function to the basis functions we already have
      d.cand<-getd(X.cand,v[i-1,],s2[i-1],tau2,y)

      llik.alpha <- .5*log(1/tau2) + d.cand$d - d.curr$d # calculate the log likelihood ratio
      lprior.alpha <- ( # log prior ratio
        log(lam[i-1])-log(nbasis[i-1]+1) # nbasis
        + log(1/max.int)  # nint
        + nint.cand*log(1/2)  # signs
        + log(1/choose(p,nint.cand)) # vars
        + log(nbasis[i-1]+1) # ordering
      )
      lprop.alpha <- ( # log proposal ratio
        # probability of going from proposed to current (death)
        log(1/3) + # probability of selection a death step
          + log(1/(nbasis[i-1]+1)) # probability that this basis function is selected to kill
        # probability of going from current to proposed
        - (
          log(1/3) # probability of selecting a birth step
          + log(1/max.int) # probability of nint.cand
          + nint.cand*log(1/2) # probability of the signs we got
          + log(1/choose(p,nint.cand)) # probability of vars.cand
        )
      )

      alpha <- llik.alpha + lprior.alpha + lprop.alpha

      if(log(runif(1))<alpha){
        X.curr<-X.cand
        d.curr<-d.cand
        nbasis[i]<-nbasis[i-1]+1
        nint[i,nbasis[i]]<-nint.cand
        knots[i,nbasis[i],1:nint.cand]<-knots.cand
        signs[i,nbasis[i],1:nint.cand]<-signs.cand
        vars[i,nbasis[i],1:nint.cand]<-vars.cand
        count[1]<-count[1]+1
      }

    } else if(move.type=='death'){
      tokill<-sample(nbasis[i-1],1) # which basis function we will delete
      X.cand<-X.curr[,-(tokill+1),drop=F] # +1 to skip the intercept
      d.cand <- getd(X.cand,v[i-1,],s2[i-1],tau2,y)

      llik.alpha <- -.5*log(1/tau2) + d.cand$d - d.curr$d
      lprior.alpha <- (
        -log(lam[i-1])+log(nbasis[i-1]) # nbasis
        - log(1/max.int)  # nint
        - nint[i-1,tokill]*log(1/2)  # signs
        - log(1/choose(p,nint[i-1,tokill])) # vars
        - log(nbasis[i-1]) # ordering
      )
      lprop.alpha <- (
        # probability of going from proposed to current (birth)
        log(1/3) # probability of selecting a birth step
        + log(1/max.int) # probability of nint
        + nint[i-1,tokill]*log(1/2) # probability of the signs we got
        + log(1/choose(p,nint[i-1,tokill])) # probability of vars
        # probability of going from current to proposed
        - (
          log(1/3) # probability of selection a death step
          + log(1/nbasis[i-1]) # probability that this basis function is selected to kill
        )
      )

      alpha <- llik.alpha + lprior.alpha + lprop.alpha

      if(log(runif(1))<alpha){
        X.curr<-X.cand
        d.curr<-d.cand
        nbasis[i]<-nbasis[i-1]-1
        nint[i,]<-NA
        knots[i,,]<-signs[i,,]<-vars[i,,]<-NA
        if(nbasis[i]==0){
          nint[i,]<-NA
        } else{
          nint[i,1:nbasis[i]]<-nint[i-1,(1:nbasis[i-1])[-tokill]]
          knots[i,1:nbasis[i],]<-knots[i-1,(1:nbasis[i-1])[-tokill],]
          signs[i,1:nbasis[i],]<-signs[i-1,(1:nbasis[i-1])[-tokill],]
          vars[i,1:nbasis[i],]<-vars[i-1,(1:nbasis[i-1])[-tokill],]
        }
        count[2]<-count[2]+1
      }

    } else{

      tochange<-sample(nbasis[i-1],1) # which basis function we will change
      tochange2<-sample(nint[i-1,tochange],1) # which element in the basis function tensor product we will change
      knots.cand<-knots[i-1,tochange,1:nint[i-1,tochange]] # copy previous
      knots.cand[tochange2]<-runif(1) # change one element
      signs.cand<-signs[i-1,tochange,1:nint[i-1,tochange]]
      signs.cand[tochange2]<-sample(c(-1,1),1)
      basis<-makeBasis(signs.cand,vars[i-1,tochange,1:nint[i-1,tochange]],knots.cand,Xt)
      X.cand<-X.curr
      X.cand[,tochange+1]<-basis # +1 for intercept

      d.cand <- getd(X.cand,v[i-1,],s2[i-1],tau2,y)

      llik.alpha <- d.cand$d - d.curr$d

      alpha <- llik.alpha

      if(log(runif(1))<alpha){
        X.curr<-X.cand
        d.curr<-d.cand
        knots[i,tochange,1:nint[i,tochange]]<-knots.cand
        signs[i,tochange,1:nint[i,tochange]]<-signs.cand
        count[3]<-count[3]+1
      }

    }

    ## Gibbs steps

    lam[i]<-rgamma(1,h1+nbasis[i],h2+1)
    #curr$beta<-curr$bhat/(1+curr$beta.prec)+curr$R.inv.t%*%rnorm(curr$nc)*sqrt(curr$s2/(1+curr$beta.prec)/data$itemp.ladder[curr$temp.ind])
    beta[i,1:(nbasis[i]+1)]<-mnormt::rmnorm(1,d.curr$bhat,d.curr$V)
    res<-y-X.curr%*%t(beta[i,1:(nbasis[i]+1),drop=F])
    v[i,]<-rgamma(n,(nu+1)/2,nu/2 + .5/s2[i-1]*res^2)
    s2[i]<-1/rgamma(1,n/2+g1,rate=g2+.5*sum(v[i,]*res^2))

    d.curr<-getd(X.curr,v[i,],s2[i],tau2,y)

    if(verbose == TRUE) {
      if(i%%ticker==0) {
        cat(timestamp(quiet = T),' nmcmc: ',i,' nbasis: ',nbasis[i],'\n')
      }
    }
  }
  return(list(X=X.curr,b=d.curr$bhat,count=count,knots=knots,signs=signs,vars=vars,nint=nint,nbasis=nbasis,beta=beta,s2=s2,lam=lam,v=v))
}
