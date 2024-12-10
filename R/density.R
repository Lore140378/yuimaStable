dsymStable<-function(x,beta0,nodes,logw, N=100, laguerreW=NULL,dens=numeric(length(x)),pos=75,posInt=0.5){
  if(beta0==1){
    return(1/(pi*(1+x^2)))
  }

  pneg0 <- x < -pos
  xneg0 <- x[x < -pos]
  pneg <- x< -posInt & x >= -pos
  xneg<- x[pneg]
  ppos <-x> posInt &  x <= pos
  xpos<-  x[ppos]
  ppos0<-x > pos
  xpos0 <- x[ppos0]
  pzero <- x>= -posInt & x<= posInt
  xzero <- x[pzero]

  if(beta0==1){
    return(1/(pi*(1+x^2)))
  }
  #N <- 2000
  theta<- pi/4+pi/4*nodes
  #theta<- seq(0,pi/2,length.out=N+1)+pi/(4*N)
  costheta <- cos(theta)
  exp1<-(beta0/(beta0-1))
  V <- (costheta/sin(beta0*theta))^(exp1)
  V <- V*cos((beta0-1)*theta)/costheta
  cons <- beta0/(pi*abs(beta0-1))
  
  #rep(log(pi/(2*N)),N+1)
  size <- length(V)
  ext <- rep(exp1,size)
  #correction <- log(beta0*sin(beta0*pi/(2))*gamma(beta0)/(abs(pos)^(1+beta0)*pi)) - log(mytail(1,pos,ext,V,logw)*cons/(pos)*pi/4)
  if(length(xneg0)>0){
    dens[pneg0]<-beta0*sin(beta0*pi/(2))*gamma(beta0)/(abs(xneg0)^(1+beta0)*pi)#*exp(-correction)
    #dens[pneg0]<- libstableR::stable_pdf(xneg0,pars=c(beta0,0,1,0))
  }
  
  if (length(xneg)>0){
    dens[pneg] <- mytail(length(xneg),-xneg,ext,V,logw)*cons/(-xneg)*pi/4
    #dens[pneg] <- libstableR::stable_pdf(-xneg,pars=c(beta0,0,1,0))
  }
  if(length(xzero)>0){
    dens[pzero] <- stable_internal_density(xzero, beta0, N, laguerreW = laguerreW)
    # dens[pzero] <- libstableR::stable_pdf(xzero,pars=c(beta0,0,1,0))
    #dens[pzero] <- stabledist::dstable(xzero,beta0,0)
  }
  if(length(xpos)>0){
    dens[ppos] <- mytail(length(xpos),xpos,ext,V,logw)*cons/xpos*pi/4
    #dens[ppos] <- libstableR::stable_pdf(xpos,pars=c(beta0,0,1,0))
  }
  if(length(xpos0)>0){
    dens[ppos0]<-beta0*sin(beta0*pi/(2))*gamma(beta0)/(abs(xpos0)^(1+beta0)*pi)#*exp(-correction)
    #dens[ppos0]<-libstableR::stable_pdf(xpos0,pars=c(beta0,0,1,0))
  }
  
  return(dens)
}

f_betal1<-function(x, beta0, nodes, weights){
  a<-nodes^(1/beta0-1)/beta0*weights
  b <- nodes^(1/beta0)
  lengthx<-length(x)
  return(mysum(lengthx, x, a, b, pi)) 
}

f_betag1<-function(x, beta0, nodes, weights){
  b <- nodes
  lengthx<-length(x)
  a <- matrix(exp(-nodes^beta0+nodes+log(weights)),1,length(weights))
  #d <- matrix(1,1, length(weights))
  
  return(mysum(lengthx, x, a, b, pi))
  #return(mysum(length(x), x, t(weights), nodes^(1/beta0), pi*beta0))
  #return(mysum(length(x), x, t(weights), nodes^(1/beta0), pi*beta0))
}


stable_internal_density <- function(x, beta0, N=2000, laguerreW = NULL){
  if (beta0==2){
    return(dnorm(x))
  }
  if(beta0==1){
    return(dcauchy(x))
  }
  if(beta0<1 ){
    if(is.null(laguerreW)){
      aa<-gauss.quad(N,"laguerre",0)
      cond <- aa$weights>0
      nodes <-aa$nodes[cond]
      weights <- aa$weights[cond] 
    }else{
      aa <- laguerreW
      cond <- aa$weights>0
      nodes <-aa$nodes[cond]
      weights <- aa$weights[cond] 
    }
    return(f_betal1(x, beta0, nodes, weights))
  }else{
    if(is.null(laguerreW)){
      aa<-gauss.quad(N,"laguerre",0)
      cond <- aa$weights>0
      nodes <-aa$nodes[cond]
      weights <- aa$weights[cond] 
    }else{
      aa <- laguerreW
      cond <- aa$weights>0
      nodes <-aa$nodes[cond]
      weights <- aa$weights[cond] 
    }
    # aa<-gauss.quad(N,"laguerre",1/beta0-1)
    # cond <- aa$weights>0
    # nodes <-aa$nodes[cond]
    # weights <- aa$weights[cond]
    return(f_betag1(x, beta0, nodes, weights))
  }
}

