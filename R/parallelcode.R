# library(Prova)
# library(doParallel)
#library(iterators)
#library(tictoc)

ParStable_pdf <- function(x, N,beta0,aa=NULL,nodes,logw,pos=75,posInt=1){
  mydens <- matrix(0,nrow = nrow(x), ncol=ncol(x))
  i<-1
  r<-foreach(i=1:nrow(x),
             #.combine = rbind,.export = c("dsymStable","mytail","stable_internal_density","f_betag1","f_betal1","mysum"), .packages=c("libstableR","statmod","stabledist"),mydens)%dopar%{
             #.combine = rbind, .packages=c("yuimaStable","libstableR"),mydens)%dopar%{
             .combine = rbind, .packages=c("yuimaStable"),mydens)%dopar%{
              dsymStable(x[i,],beta0,nodes,logw,N=N, laguerreW=aa,pos=pos,posInt=posInt)
               #libstableR::stable_pdf(x[d,],pars=c(beta0,0,1,0))
             }
  mydens1 <- as.numeric(r)
}


