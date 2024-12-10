#load("data/mydata.RData")
internalfunC1<-function(v){
  mytab <- NULL
  data("mydata",package="yuimaStable")
  myfun<-approxfun(x=mytab[,1],y=mytab[,2])
  return(myfun(v))
}
internalfunC2<-function(v){
  mytab <- NULL
  data("mydata",package="yuimaStable")
  myfun<-approxfun(x=mytab[,1],y=mytab[,3])
  return(myfun(v))
}


internal_minuslog_stable <- function(par, Xt, A_coef, C_coef, h, Nint, 
      N, N_al, N_gam, aa=NULL,nodes,logw, W=W,parallel = F,dens=NULL,
      pos=75, posInt=1, num_of_cores=1, aa_alter=NULL){
  # alpha beta gamma
  alpha <- par[1:N_al]
  #
  beta0 <- par[N_al+1]
  #beta0 <- 1.5
  
  gam <- par[1:N_gam+N_al+1]
  # gam <- par[1:N_gam+N_al]
  dum <-   Xt[1:Nint]
  A_drift <- A_coef(alpha, dum)
  C_jump <- C_coef(gam, dum)
  Epsilon <- (diff(Xt)-h*A_drift)/(h^(1/beta0)*C_jump)
  if(parallel){
    
    Epsilon <- matrix(Epsilon,nrow = num_of_cores, ncol=Nint/num_of_cores)
    if(beta0>1.07){
      logStable<-log(ParStable_pdf(x=Epsilon, N=N,beta0=beta0, aa=aa, nodes=nodes, logw=logw,pos=pos,posInt=posInt))
    }else{
      logStable<-log(ParStable_pdf(x=Epsilon, N=N,beta0=beta0, aa=aa_alter, nodes=nodes, logw=logw,pos=pos,posInt=pos))
    }
  }else{
    logStable<- log(dsymStable(Epsilon,beta0,nodes=nodes,logw=logw, N=N, laguerreW=aa,dens=dens,pos=pos,posInt=posInt))
    #logStable <-log(pmax(stable_pdf(Epsilon,pars=c(beta0, 0,1,0)),0))
  }
  
  loglik <- sum((- 1/beta0*log(h) - log(C_jump[1:Nint]) + logStable)*W) 
  
  #cat("\n", c(loglik, par))
  return(-loglik)
}

internal_qsmle_code <- function(Data, start, A_coef, C_coef, N_al, N_gam,NumericalHess=FALSE, Deriv_A,logDeriv_C,  method = "Nelder-Mead", 
                                lower = NULL, upper = NULL, info=NULL, joint=TRUE, Todorov=FALSE, N=50, 
                                num_of_cores= detectCores(),parallel = F,pos=75, posInt=1,
                                r0=0.9,  maxcount=100,scale=2, onlypar = FALSE,W=1, Residual = T,aa_alt=NULL){
  #timeT <- index(Data)
  timeT <- time(Data)
  Nint <- length(timeT)-1
  Deltat <- (tail(timeT,1)-timeT[1])/(Nint)
  Xt <- as.numeric(Data)
  #Nlaguerre<-max(N,2500)
  Nlaguerre<-N
  aa<- gauss.quad(Nlaguerre,"laguerre",0)
  cond<-aa$weights>0
  aa$nodes<-aa$nodes[cond]
  aa$weights<-aa$weights[cond]
  if(is.null(aa_alt)){
    aa_alt<- gauss.quad(max(c(Nlaguerre,10000)),"laguerre",0)
    cond_alt<-aa_alt$weights>0
    aa_alt$nodes<-aa_alt$nodes[cond_alt]
    aa_alt$weights<-aa_alt$weights[cond_alt]
  }
  poly <- gauss.quad(N)
  nodes <- poly$nodes
  logw <-log(poly$weights)
  if(joint){
    if(parallel){
      cl <- makePSOCKcluster(num_of_cores)
      registerDoParallel(cl,cores=num_of_cores)
      if(method!="L-BFGS-B" & method!="Brent"){
        res<-optim(start, internal_minuslog_stable, method = method,
                   Xt = Xt, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel, pos=pos, posInt=posInt,num_of_cores=num_of_cores,
                   aa_alter = aa_alt)
      }else{
        res<-optim(start, internal_minuslog_stable, method = method,lower =lower, upper=upper,
                   Xt= Xt, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,
                   aa=aa,W=W,nodes=nodes,logw=logw, parallel = parallel,
                   pos=pos, posInt=posInt,num_of_cores=num_of_cores,aa_alter = aa_alt)
      }
      stopCluster(cl)
    }else{
      if(method!="L-BFGS-B" & method!="Brent"){
        res<-optim(start, internal_minuslog_stable, method = method,
                   Xt = Xt, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel)
      }else{
        res<-optim(start, internal_minuslog_stable, method = method,lower =lower, upper=upper,
                   Xt= Xt, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,
                   aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel)
      }
    }
    res_par <- res$par
    #return(res)
  }else{
    if(!Todorov){
      gam <- start[1:N_gam+N_al+1]
      beta0 <- start[N_al+1]
      startg_b <- c(beta0,gam)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_first<-optim(startg_b, internal_minuslog_stable_g_b, method = method,
                         Xt = Xt, C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
        #Xt, C_coef, h, Nint,N, N_gam,
      }else{
        res_first<-optim(startg_b, internal_minuslog_stable_g_b, method = method,lower =lower, upper=upper,
                         Xt= Xt,  C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
      }
      start_a <- start[1:N_al]
      gam_hat<- res_first$par[1:N_gam+1]
      C_jump_hat <- C_coef(gam_hat,Xt)
      beta_hat <- tail(res_first$par,1L)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,
                          gam=gam_hat,beta0=beta_hat,
                          Xt = Xt, A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  aa=aa) 
        #gam, beta0, Xt, A_coef, C_jump, h, Nint, N, N_al,
      }else{
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,lower =lower, upper=upper,
                          gam=gam_hat,beta0=beta_hat,
                          Xt= Xt,  A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  
                          aa=aa) 
      }
      res_par <- c(res_second$par,res_first$par)
    }else{
      gam <- start[1:N_gam+N_al+1]
      #
      beta0 <- start[N_gam+N_al+1]
      beta_hat <- Todorov_beta_est(DeltaX = diff(Xt),
                                   r0=r0, maxcount=maxcount, scale=scale)
      # internal_minuslog_stable_g
      # (par=gam, beta0= ,Xt, C_coef, h, Nint,N, N_gam, aa=NULL)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_first<-optim(gam, internal_minuslog_stable_g, method = method,beta0=beta_hat,
                         Xt = Xt, C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,aa=aa) 
        #Xt, C_coef, h, Nint,N, N_gam,
      }else{
        res_first<-optim(gam, internal_minuslog_stable_g, method = method, beta0=beta_hat,
                         lower =lower, upper=upper,
                         Xt= Xt,  C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
      }
      start_a <- start[1:N_al]
      gam_hat<- res_first$par
      C_jump_hat <- C_coef(gam_hat,Xt)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,
                          gam=gam_hat,beta0=beta_hat,
                          Xt = Xt, A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  aa=aa) 
        #gam, beta0, Xt, A_coef, C_jump, h, Nint, N, N_al,
      }else{
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,lower =lower, upper=upper,
                          gam=gam_hat,beta0=beta_hat,
                          Xt= Xt,  A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  
                          aa=aa) 
      }
      # res_par<- c(res_second$par,res_first$par,beta_hat)
      res_par<- c(res_second$par,beta_hat,res_first$par)
    }
  }
  ############# Computation of VCOV matrix ######################
  if(onlypar){
    res <-list(par <- res_par)  
  }else{
    beta_hat <- res_par[N_al+1]
    #beta_hat <- 1.5
    alpha_hat<- res_par[1:N_al]
    gam_hat <- res_par[1:N_gam+N_al+1] 
    #gam_hat <- res_par[1:N_gam+N_al] 
    estim_drift <- A_coef(alpha_hat,Xt)[-length(Xt)]
    estim_jump <- C_coef(gam_hat,Xt)[-length(Xt)]
    
    u_n_y <- diag(c(rep(sqrt(Nint)*h^(1-1/beta_hat),length(alpha_hat)),
                    sqrt(Nint)*log(1/h)/beta_hat^2,
                    rep(sqrt(Nint),length(gam_hat))))
    
    # u_n_y <- diag(c(rep(sqrt(Nint)*h^(1-1/beta_hat),length(alpha_hat)),
    #                 rep(sqrt(Nint),length(gam_hat))))
    
    if(NumericalHess){
      
      Hess<-optimHess(par=res_par, fn=internal_minuslog_stable, gr=NULL,
                      Xt=Xt,A_coef=A_coef,C_coef=C_coef,
                      h=Deltat,N =Nlaguerre,
                      Nint=Nint,N_al=N_al,N_gam=N_gam, aa=aa, W=W)
      Hess[1:N_al,1:(N_gam+1)+N_al]<-0
      Hess[1:(N_gam+1)+N_al,1:N_al]<-0
      Fish <- solve(u_n_y)%*%Hess%*%solve(u_n_y)
      vcov <- solve(Fish)
    }else{
      Deriv_A <- Deriv_A(alpha_hat, Xt)
      logDeriv_C <- logDeriv_C(gam_hat, Xt)
      pos <- 1:(length(Xt)-1)
#      first_comp<-matrix(0, N_al,N_al)
#      second_comp <- matrix(0,N_gam+1,N_gam+1)
      #second_comp <- matrix(0,N_gam,N_gam)
      #Third_comp <- matrix(0,1,1)
      #mixed_der <- matrix(0,length(gam),1)
      C1 <- internalfunC1(beta_hat)
      C2 <- internalfunC2(beta_hat)
      # C1 <- 0.5
      # C2 <-0.5
      
      
      if(length(W)==1){
        # for(i in pos){
        #   dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2
        #   
        #   dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
        #                 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))
        #   first_comp<- first_comp+dumA
        #   second_comp <-  second_comp+dumC
        # }
        provamy <- CompMatAsym(tail(pos,1L), Deriv_A, estim_jump, logDeriv_C)
        first_comp<-provamy$avec
        second_comp<-provamy$amat
        
        I_hat_al <- C1/Nint*first_comp
        I_hat_gam <- C2/Nint*second_comp
        Fish <- matrix(0,length(res_par),length(res_par))
        Fish[1:N_al,1:N_al]<-I_hat_al
        newpos <- 1:dim(I_hat_gam)[1]+N_al
        Fish[newpos,newpos]<-I_hat_gam 
        vcov <- solve(Fish) 
      }else{
        provamy <- CompMatAsym1(tail(pos,1L), Deriv_A, estim_jump, logDeriv_C,W)
        
        first_comp <-provamy$avec
        first_comp2 <- provamy$avec2
        second_comp<-provamy$amat
        second_comp2<-provamy$amat2
        
        #first_comp2<-matrix(0, N_al,N_al)
        #second_comp2 <- matrix(0,N_gam+1,N_gam+1)
        
        #second_comp2 <- matrix(0,N_gam,N_gam)
        # for(i in pos){
        #   dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2*W[i]
        #   dumA2 <- dumA*W[i] 
        #   dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
        #                 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))*W[i]
        #   
        #   #           dumC <- logDeriv_C[,i]%*%t(logDeriv_C[,i])*W[i] 
        #   dumC2 <- dumC*W[i]
        #   first_comp<- first_comp+dumA
        #   first_comp2 <- first_comp2+dumA2
        #   second_comp <-  second_comp+dumC
        #   second_comp2 <- second_comp2+dumC2
        # }
        InvGamma_al_n <- solve(C1/Nint*first_comp) # Inverse_of_Gamma_hat_alpha_n
        Sigma_al_n <- C1/Nint*first_comp2  # Sigma_alpha_hat_n
        InvGamma_b_g_n <- solve(C2/Nint*second_comp)  
        Sigma_b_g_n <- C2/Nint*second_comp2
        vcov <- matrix(0,length(res_par),length(res_par))
        vcov[1:N_al,1:N_al]<- InvGamma_al_n%*%Sigma_al_n%*%InvGamma_al_n #  Inverse_of_Gamma_hat_alpha_n %*% Sigma_alpha_hat_n %*% Inverse_of_Gamma_hat_alpha_n 
        newpos <- 1:dim(Sigma_b_g_n)[1]+N_al
        vcov[newpos,newpos]<-InvGamma_b_g_n%*%Sigma_b_g_n%*%InvGamma_b_g_n 
        Fish <- solve(vcov) # AsimptoticInfoMatrix =  (Inverse_of_Gamma_hat_alpha_n %*% Sigma_alpha_hat_n %*% Inverse_of_Gamma_hat_alpha_n)^-1 
        # InvF1 <- solve(first_comp)
        # I_hat_al <- C1/Nint*solve(InvF1%*%first_comp2%*%InvF1)
        # InvS1 <- solve(second_comp)
        # I_hat_gam <- C2/Nint*solve(InvS1%*%second_comp2%*%InvS1)
        #InvF1 <- solve(first_comp)
        # I_hat_al <- C1/Nint*(first_comp%*%solve(first_comp2)%*%first_comp)
        # #InvS1 <- solve(second_comp)
        # I_hat_gam <- C2/Nint*solve(second_comp%*%solve(second_comp2)%*%second_comp)
        
      }
      
      
    } 
    loglik = NULL
    if(joint){loglik <- -res$value}
    if(joint){
    res <- list(par = res_par, 
                PrecMat = Fish,
                vcov = vcov, 
                D_n =u_n_y,
                loglik = loglik,
                optim  = res)
    }else{
      res <- list(par = res_par, 
                  PrecMat = Fish,
                  vcov = vcov, 
                  D_n =u_n_y,
                  loglik = loglik
                  )
    }
    if(!Residual){
      return(res)
    }
  }
  ###################
  # Filtering Noise #
  ###################
  # alpha_hat <-res_par[1:N_al]
  # gam_hat <- res_par[1:N_gam+N_al]
  # beta_hat <- tail(res_par,1L)
  Epsilon_Delta <-  diff(Xt)
  
  resi <- (Epsilon_Delta-Deltat*estim_drift)/estim_jump
  
  #LevyDelta <- zoo(cumsum(c(0,resi)),order.by=index(Data))
  LevyDelta <- zoo(cumsum(c(0,resi)),order.by=time(Data))
  unitLevy <- as.numeric(diff(na.approx(LevyDelta, xout=seq(0,floor(tail(timeT,1))), method = "constant")))
  res$unitNoise <- unitLevy
  res$DeltaNoise <- resi
  return(res)
  
}



internal_minuslog_stable_fixed_beta0 <- function(par, Xt,beta0, A_coef, C_coef, h, Nint, 
                                     N, N_al, N_gam, aa=NULL,nodes,logw, W=W,parallel = F,dens=NULL,
                                     pos=75, posInt=1, num_of_cores=1, aa_alter=NULL){
  # alpha beta gamma
  alpha <- par[1:N_al]
  #
  #beta0 <- par[N_al+1]
  #beta0 <- 1.5
  
  gam <- par[1:N_gam+N_al]
  # gam <- par[1:N_gam+N_al]
  dum <-   Xt[1:Nint]
  A_drift <- A_coef(alpha, dum)
  C_jump <- C_coef(gam, dum)
  Epsilon <- (diff(Xt)-h*A_drift)/(h^(1/beta0)*C_jump)
  if(parallel){
    
    Epsilon <- matrix(Epsilon,nrow = num_of_cores, ncol=Nint/num_of_cores)
    if(beta0>1.07){
      logStable<-log(ParStable_pdf(x=Epsilon, N=N,beta0=beta0, aa=aa, nodes=nodes, logw=logw,pos=pos,posInt=posInt))
    }else{
      logStable<-log(ParStable_pdf(x=Epsilon, N=N,beta0=beta0, aa=aa_alter, nodes=nodes, logw=logw,pos=pos,posInt=pos))
    }
  }else{
    logStable<- log(dsymStable(Epsilon,beta0,nodes=nodes,logw=logw, N=N, laguerreW=aa,dens=dens,pos=pos,posInt=posInt))
    #logStable <-log(pmax(stable_pdf(Epsilon,pars=c(beta0, 0,1,0)),0))
  }
  
  loglik <- sum((- 1/beta0*log(h) - log(C_jump[1:Nint]) + logStable)*W) 
  
  #cat("\n", c(loglik, par))
  return(-loglik)
}



internal_qsmle_code_Fixed_beta0 <- function(Data, start, beta0, A_coef, C_coef, N_al, N_gam,NumericalHess=FALSE, Deriv_A,logDeriv_C,  method = "Nelder-Mead", 
                                lower = NULL, upper = NULL, info=NULL, joint=TRUE, Todorov=FALSE, N=50, 
                                num_of_cores= detectCores(),parallel = F,pos=75, posInt=1,
                                r0=0.9,  maxcount=100,scale=2, onlypar = FALSE,W=1, Residual = T,aa_alt=NULL){
  timeT <- index(Data)
  Nint <- length(timeT)-1
  Deltat <- (tail(timeT,1)-timeT[1])/(Nint)
  Xt <- as.numeric(Data)
  #Nlaguerre<-max(N,2500)
  Nlaguerre<-N
  aa<- gauss.quad(Nlaguerre,"laguerre",0)
  cond<-aa$weights>0
  aa$nodes<-aa$nodes[cond]
  aa$weights<-aa$weights[cond]
  if(is.null(aa_alt)){
    aa_alt<- gauss.quad(max(c(Nlaguerre,10000)),"laguerre",0)
    cond_alt<-aa_alt$weights>0
    aa_alt$nodes<-aa_alt$nodes[cond_alt]
    aa_alt$weights<-aa_alt$weights[cond_alt]
  }
  poly <- gauss.quad(N)
  nodes <- poly$nodes
  logw <-log(poly$weights)
  if(joint){
    if(parallel){
      cl <- makePSOCKcluster(num_of_cores)
      registerDoParallel(cl,cores=num_of_cores)
      if(method!="L-BFGS-B" & method!="Brent"){
        res<-optim(start, internal_minuslog_stable_fixed_beta0, method = method,
                   Xt = Xt, beta0=beta0, A_coef = A_coef, C_coef = C_coef, 
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel, pos=pos, posInt=posInt,num_of_cores=num_of_cores,
                   aa_alter = aa_alt)
      }else{
        res<-optim(start, internal_minuslog_stable_fixed_beta0, method = method,lower =lower, upper=upper,
                   Xt= Xt, beta0=beta0, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,
                   aa=aa,W=W,nodes=nodes,logw=logw, parallel = parallel,
                   pos=pos, posInt=posInt,num_of_cores=num_of_cores,aa_alter = aa_alt)
      }
      stopCluster(cl)
    }else{
      if(method!="L-BFGS-B" & method!="Brent"){
        res<-optim(start, internal_minuslog_stable_fixed_beta0, method = method,
                   Xt = Xt,beta0=beta0, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel)
      }else{
        res<-optim(start, internal_minuslog_stable_fixed_beta0, method = method,lower =lower, upper=upper,
                   Xt= Xt, beta0=beta0, A_coef = A_coef, C_coef = C_coef,
                   h = Deltat, Nint = Nint, N = Nlaguerre, N_al =N_al, N_gam = N_gam,
                   aa=aa,W=W,
                   nodes=nodes,logw=logw, parallel = parallel)
      }
    }
    res_par <- res$par
    #return(res)
  }else{
    if(!Todorov){
      gam <- start[1:N_gam+N_al+1]
      beta0 <- start[N_al+1]
      startg_b <- c(beta0,gam)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_first<-optim(startg_b, internal_minuslog_stable_g_b, method = method,
                         Xt = Xt, C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
        #Xt, C_coef, h, Nint,N, N_gam,
      }else{
        res_first<-optim(startg_b, internal_minuslog_stable_g_b, method = method,lower =lower, upper=upper,
                         Xt= Xt,  C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
      }
      start_a <- start[1:N_al]
      gam_hat<- res_first$par[1:N_gam+1]
      C_jump_hat <- C_coef(gam_hat,Xt)
      beta_hat <- tail(res_first$par,1L)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,
                          gam=gam_hat,beta0=beta_hat,
                          Xt = Xt, A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  aa=aa) 
        #gam, beta0, Xt, A_coef, C_jump, h, Nint, N, N_al,
      }else{
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,lower =lower, upper=upper,
                          gam=gam_hat,beta0=beta_hat,
                          Xt= Xt,  A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  
                          aa=aa) 
      }
      res_par <- c(res_second$par,res_first$par)
    }else{
      gam <- start[1:N_gam+N_al+1]
      #
      beta0 <- beta0
      beta_hat <- Todorov_beta_est(DeltaX = diff(Xt),
                                   r0=r0, maxcount=maxcount, scale=scale)
      # internal_minuslog_stable_g
      # (par=gam, beta0= ,Xt, C_coef, h, Nint,N, N_gam, aa=NULL)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_first<-optim(gam, internal_minuslog_stable_g, method = method,beta0=beta_hat,
                         Xt = Xt, C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,aa=aa) 
        #Xt, C_coef, h, Nint,N, N_gam,
      }else{
        res_first<-optim(gam, internal_minuslog_stable_g, method = method, beta0=beta_hat,
                         lower =lower, upper=upper,
                         Xt= Xt,  C_coef = C_coef, 
                         h = Deltat, Nint = Nint, N = Nlaguerre,  N_gam = N_gam,
                         aa=aa) 
      }
      start_a <- start[1:N_al]
      gam_hat<- res_first$par
      C_jump_hat <- C_coef(gam_hat,Xt)
      if(method!="L-BGFS-B" & method!="Brent"){
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,
                          gam=gam_hat,beta0=beta_hat,
                          Xt = Xt, A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  aa=aa) 
        #gam, beta0, Xt, A_coef, C_jump, h, Nint, N, N_al,
      }else{
        res_second<-optim(start_a, internal_minuslog_stable_alpha, method = method,lower =lower, upper=upper,
                          gam=gam_hat,beta0=beta_hat,
                          Xt= Xt,  A_coef=A_coef, C_jump=C_jump_hat, 
                          h = Deltat, Nint = Nint, N = Nlaguerre, N_al=N_al,  
                          aa=aa) 
      }
      # res_par<- c(res_second$par,res_first$par,beta_hat)
      res_par<- c(res_second$par,beta_hat,res_first$par)
    }
  }
  ############# Computation of VCOV matrix ######################
  if(onlypar){
    res <-list(par <- res_par)  
  }else{
    beta_hat <- res_par[N_al+1]
    #beta_hat <- 1.5
    alpha_hat<- res_par[1:N_al]
    gam_hat <- res_par[1:N_gam+(N_al+1)] 
    #gam_hat <- res_par[1:N_gam+N_al] 
    estim_drift <- A_coef(alpha_hat,Xt)[-length(Xt)]
    estim_jump <- C_coef(gam_hat,Xt)[-length(Xt)]
    
    # u_n_y <- diag(c(rep(sqrt(Nint)*h^(1-1/beta_hat),length(alpha_hat)),
    #                 sqrt(Nint)*log(1/h)/beta_hat^2,
    #                 rep(sqrt(Nint),length(gam_hat))))
    
    u_n_y <- diag(c(rep(sqrt(Nint)*h^(1-1/beta_hat),length(alpha_hat)),
                    rep(sqrt(Nint),length(gam_hat))))
    
    if(NumericalHess){
      
      Hess<-optimHess(par=res_par, fn=internal_minuslog_stable, gr=NULL,
                      Xt=Xt,A_coef=A_coef,C_coef=C_coef,
                      h=Deltat,N =Nlaguerre,
                      Nint=Nint,N_al=N_al,N_gam=N_gam, aa=aa, W=W)
      Hess[1:N_al,1:(N_gam+1)+N_al]<-0
      Hess[1:(N_gam+1)+N_al,1:N_al]<-0
      Fish <- solve(u_n_y)%*%Hess%*%solve(u_n_y)
      vcov <- solve(Fish)
    }else{
      Deriv_A <- Deriv_A(alpha_hat, Xt)
      logDeriv_C <- logDeriv_C(gam_hat, Xt)
      pos <- 1:(length(Xt)-1)
      #      first_comp<-matrix(0, N_al,N_al)
      #      second_comp <- matrix(0,N_gam+1,N_gam+1)
      #second_comp <- matrix(0,N_gam,N_gam)
      #Third_comp <- matrix(0,1,1)
      #mixed_der <- matrix(0,length(gam),1)
      C1 <- internalfunC1(beta0)
      C2 <- internalfunC2(beta0)
      # C1 <- 0.5
      # C2 <-0.5
      
      
      if(length(W)==1){
        # for(i in pos){
        #   dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2
        #   
        #   dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
        #                 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))
        #   first_comp<- first_comp+dumA
        #   second_comp <-  second_comp+dumC
        # }
        provamy <- CompMatAsym(tail(pos,1L), Deriv_A, estim_jump, logDeriv_C)
        first_comp<-provamy$avec
        second_comp<-provamy$amat[2:N_gam,2:N_gam]
        
        I_hat_al <- C1/Nint*first_comp
        I_hat_gam <- C2/Nint*second_comp
        Fish <- matrix(0,length(res_par),length(res_par))
        Fish[1:N_al,1:N_al]<-I_hat_al
        newpos <- 1:dim(I_hat_gam)[1]+N_al
        Fish[newpos,newpos]<-I_hat_gam 
        vcov <- solve(Fish) 
      }else{
        provamy <- CompMatAsym1(tail(pos,1L), Deriv_A, estim_jump, logDeriv_C,W)
        
        first_comp <-provamy$avec
        first_comp2 <- provamy$avec2
        second_comp<-provamy$amat[2:(N_gam+1),2:(N_gam+1)]
        second_comp2<-provamy$amat2[2:(N_gam+1),2:(N_gam+1)]
        
        #first_comp2<-matrix(0, N_al,N_al)
        #second_comp2 <- matrix(0,N_gam+1,N_gam+1)
        
        #second_comp2 <- matrix(0,N_gam,N_gam)
        # for(i in pos){
        #   dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2*W[i]
        #   dumA2 <- dumA*W[i] 
        #   dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
        #                 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))*W[i]
        #   
        #   #           dumC <- logDeriv_C[,i]%*%t(logDeriv_C[,i])*W[i] 
        #   dumC2 <- dumC*W[i]
        #   first_comp<- first_comp+dumA
        #   first_comp2 <- first_comp2+dumA2
        #   second_comp <-  second_comp+dumC
        #   second_comp2 <- second_comp2+dumC2
        # }
        InvGamma_al_n <- solve(C1/Nint*first_comp) # Inverse_of_Gamma_hat_alpha_n
        Sigma_al_n <- C1/Nint*first_comp2  # Sigma_alpha_hat_n
        InvGamma_b_g_n <- solve(C2/Nint*second_comp)  
        Sigma_b_g_n <- C2/Nint*second_comp2
        vcov <- matrix(0,length(res_par),length(res_par))
        vcov[1:N_al,1:N_al]<- InvGamma_al_n%*%Sigma_al_n%*%InvGamma_al_n #  Inverse_of_Gamma_hat_alpha_n %*% Sigma_alpha_hat_n %*% Inverse_of_Gamma_hat_alpha_n 
        newpos <- 1:dim(Sigma_b_g_n)[1]+N_al
        vcov[newpos,newpos]<-InvGamma_b_g_n%*%Sigma_b_g_n%*%InvGamma_b_g_n 
        Fish <- solve(vcov) # AsimptoticInfoMatrix =  (Inverse_of_Gamma_hat_alpha_n %*% Sigma_alpha_hat_n %*% Inverse_of_Gamma_hat_alpha_n)^-1 
        # InvF1 <- solve(first_comp)
        # I_hat_al <- C1/Nint*solve(InvF1%*%first_comp2%*%InvF1)
        # InvS1 <- solve(second_comp)
        # I_hat_gam <- C2/Nint*solve(InvS1%*%second_comp2%*%InvS1)
        #InvF1 <- solve(first_comp)
        # I_hat_al <- C1/Nint*(first_comp%*%solve(first_comp2)%*%first_comp)
        # #InvS1 <- solve(second_comp)
        # I_hat_gam <- C2/Nint*solve(second_comp%*%solve(second_comp2)%*%second_comp)
        
      }
      
      
    } 
    loglik = NULL
    if(joint){loglik <- -res$value}
    if(joint){
      res <- list(par = res_par, 
                  PrecMat = Fish,
                  vcov = vcov, 
                  D_n =u_n_y,
                  loglik = loglik,
                  optim  = res)
    }else{
      res <- list(par = res_par, 
                  PrecMat = Fish,
                  vcov = vcov, 
                  D_n =u_n_y,
                  loglik = loglik
      )
    }
    if(!Residual){
      return(res)
    }
  }
  ###################
  # Filtering Noise #
  ###################
  # alpha_hat <-res_par[1:N_al]
  # gam_hat <- res_par[1:N_gam+N_al]
  # beta_hat <- tail(res_par,1L)
  Epsilon_Delta <-  diff(Xt)
  
  resi <- (Epsilon_Delta-Deltat*estim_drift)/estim_jump
  
  LevyDelta <- zoo(cumsum(c(0,resi)),order.by=index(Data))
  unitLevy <- as.numeric(diff(na.approx(LevyDelta, xout=seq(0,floor(tail(timeT,1))), method = "constant")))
  res$unitNoise <- unitLevy
  res$DeltaNoise <- resi
  return(res)
  
}
qsle <- function(model, start, beta0=1.5, data, 
                 method = "Nelder-Mead", 
                 upper = NULL, 
                 lower = NULL, 
                 numericalinfo = NULL){
  if(is.zoo(data)){
    data <- data
    originalData<- setData(data)
  }else{
    if(is(data,"yuima.data")){
      originalData <- data
      data <- get.zoo.data(data)[[1]]
    }
  }
  if(is.list(start)){
    start<-unlist(start)
  }else{
    stop("start must be a list")
  }
  A_coef <- model@drift_func
  C_coef <- model@jump_func
  Deriv_A <- model@Deriv_drift
  logDeriv_C <- model@Deriv_logjump
  start <- start[c(model@parameter@drift,model@parameter@jump)]
  names(beta0)<-"beta0"
  N_al <- length(model@parameter@drift)
  N_gam <- length(model@parameter@jump)
  start <- c(start[1:N_al],beta0,start[c((N_al+1):length(start))])
  if(is.null(numericalinfo)){
      NumericalHess <-FALSE
      info <- NULL
      parallel <- T
      joint <- T
      Todorov <- F
      N <- 175
      num_of_cores <- detectCores()-2
      pos <- 7
      posInt <- 1
      r0 <- 0.9  
      maxcount <- 100 
      scale <- 2 
      onlypar <- FALSE
      data_num <- as.numeric(data)
      W <- 1/(1+data_num[1:(length(data_num)-1)]^2) 
      Residual <- T 
      aa_alt <- NULL
  }else{
      
    NumericalHess <- numericalinfo@NumericalHess 
    info <- numericalinfo@info
    parallel <- numericalinfo@parallel
    joint <- numericalinfo@joint
    Todorov <- numericalinfo@Todorov
    parallel <- numericalinfo@parallel
    N <- numericalinfo@N
    num_of_cores <- numericalinfo@num_of_cores
    pos <- numericalinfo@pos
    posInt <- numericalinfo@posInt
    r0 <- numericalinfo@r0 
    maxcount <- numericalinfo@maxcount 
    scale <- numericalinfo@scale 
    onlypar <- numericalinfo@onlypar
    data_num <- as.numeric(data)
    W <- numericalinfo@W(data_num) 
    Residual <- numericalinfo@Residual 
    aa_alt <- numericalinfo@aa_alt
   
  }
  res <- internal_qsmle_code(Data=data, start=start, A_coef = A_coef, C_coef = C_coef, 
                             N_al, N_gam, NumericalHess=NumericalHess, Deriv_A=Deriv_A,
                             logDeriv_C = logDeriv_C,  method = method, 
                             lower = lower, upper = upper, info=info, joint=joint, 
                             Todorov=Todorov, N=N, 
                             num_of_cores= num_of_cores, 
                             parallel = parallel, pos=pos, 
                             posInt=posInt,
                             r0=r0,  maxcount=maxcount, 
                             scale=scale, onlypar = onlypar, W=W, Residual = Residual, 
                             aa_alt=aa_alt)
  Incr<-setData(res$unitNoise) # Check
  res1 <- new("yuima.qmleLevy.incr", Incr.Lev=Incr,
              # data=data, 
              # model = model, 
              # coef = res$par, 
              # vcov = res$vcov,
              # min = -res$loglik,
              # nobs = length(data)
              )
  res1@Data <- originalData
  res1@model <- model
  res1@coef <- res$par
  res1@fullcoef <- res$par
  res1@vcov <- res$vcov
  res1@min <- -res$loglik
  res1@details <- res$optim
  res1@nobs <- length(data)
  res1@method <- method
  
  res2 <- list(AddInfo = res, yuima.qsle = res1)
  
  return(res2)
}