setClass("yuimaLocStableSDE",
         representation(drift_func = "function",
                        jump_func =  "function",
                        Deriv_drift = "function",
                        Deriv_logjump =  "function"),
         contains="yuima.model")

setMethod("initialize", "yuimaLocStableSDE",
          function(.Object,
                   drift_func = function(){NULL},
                   jump_func =  function(){NULL},
                   Deriv_drift = function(){NULL},
                   Deriv_logjump =  function(){NULL},
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = expression(),
                   measure=list(),
                   measure.type=character(),
                   parameter = new("model.parameter"),
                   state.variable = "lambda",
                   jump.variable = "N",
                   time.variable = "t",
                   noise.number = numeric(),
                   equation.number = numeric(),
                   dimension = numeric(),
                   solve.variable = character(),
                   xinit = expression(),
                   J.flag = logical()){
            .Object@drift_func <- drift_func
            .Object@jump_func <- jump_func
            .Object@Deriv_drift <- Deriv_drift
            .Object@Deriv_logjump <- Deriv_logjump
            .Object@drift <- drift
            .Object@diffusion <- diffusion
            .Object@hurst <- hurst
            .Object@jump.coeff <- jump.coeff
            .Object@measure <- measure
            .Object@measure.type <- measure.type
            .Object@parameter <- parameter
            .Object@state.variable <- state.variable
            .Object@jump.variable <- jump.variable
            .Object@time.variable <- time.variable
            .Object@noise.number <- noise.number
            .Object@equation.number <- equation.number
            .Object@dimension <- dimension
            .Object@solve.variable <- solve.variable
            .Object@xinit <- xinit
            .Object@J.flag <- J.flag
            return(.Object)
          })
function_int <- function(expr, par, variable,drift=TRUE){
  par1 <- paste0(c("par, ",variable), collapse="")
  inputs <- paste0("(",par1,")")
  if(drift){
   funct <- paste0("A_coeff <- function",inputs,"{ \n ")
  }else{
    funct <- paste0("Jump_coeff <- function",inputs,"{ \n ")
  }
  
  for(i  in c(1:length(par))){
    funct <- paste0(funct, paste0(par[i]," <- par[",i,"]"," \n "))
  }
  
  funct<- paste0(funct,
                 paste0("res <- ", expr),
                 " \n return(res)"," \n }"
                 )
  return(funct)
}

function_int2<-function(expr,
              par,
              variable,
              drift=TRUE){
  par1 <- paste0(c("par",variable), collapse=",")
  inputs <- paste0("(",par1,")")
  if(drift){
    funct <- paste0("Deriv_d <- function",inputs,"{ \n ")
  }else{
    funct <- paste0("Deriv_logjump <- function",inputs,"{ \n ")
  }
  for(i  in c(1:length(par))){
    funct <- paste0(funct, paste0(par[i]," <- par[",i,"]"," \n "))
  }
  expr<- as.character(expr)
  if(length(expr)>1){
     funct <- paste0(funct,"res <- rbind(",expr[1]," , ")
     if(length(expr)>2){
       for (j in c(2:(length(expr)-1))){
         funct <- paste0(funct,expr[j], " , ")
       }
     }
     funct <- paste0(funct,expr[length(expr)], " ) \n ")
  }else{
    funct <- paste0(funct,"res <- expr", " \n ")
  }
  funct <- paste0(funct,"return(res) \n }")
}

setLocStableModel<-function(drift,
                            jump.coeff, 
                            measure,
                            solve.variable = "x",
                            jump.variable = "z",
                            time.variable = "t",
                            xinit = NULL,
                            Deriv_drift = NULL,
                            Log_Deriv_Jump = NULL
){
  model <- setModel(drift=drift,
                    jump.coeff=jump.coeff,
                    measure =measure,
                    measure.type = "code",
                    state.variable = solve.variable,
                    jump.variable = jump.variable, 
                    time.variable = time.variable, 
                    solve.variable = solve.variable, 
                    xinit =xinit,
  )
  if(is(drift,"character")){
    A_coeff_int <- function_int(expr=drift,
                                par=model@parameter@drift,
                                variable=solve.variable)
    eval(parse(text = A_coeff_int))
  }else{
    A_coeff <- drift
  }
  if(is(jump.coeff,"character")){
    jump.coeff_int <- function_int(expr=jump.coeff,
                                   par=model@parameter@jump,
                                   variable=solve.variable,
                                   drift=FALSE)
    eval(parse(text = jump.coeff_int))
  }else{
    Jump_coeff <- jump.coeff
  }
  if(is.null(Deriv_drift)){
    Deriv_dr <- derivative(f= drift, 
                           var= model@parameter@drift)
    Deriv_d0 <- function_int2(expr=Deriv_dr,
                              par=model@parameter@drift,
                              variable=solve.variable,
                              drift=TRUE)
    eval(parse(text = Deriv_d0))
  }else{
    if(is(Deriv_drift,"function")){
      Deriv_d <- Deriv_drift
    }else{
      stop("missing Deriv_drift definition: either NULL or function")
    }
  }
  if(is.null(Log_Deriv_Jump)){
    Log_Deriv <- derivative(f= paste0("log(",jump.coeff,")"), 
                            var=model@parameter@jump)
    if(length(Log_Deriv)==1){
      Log_Deriv <- as.character(Simplify(Log_Deriv))
    }else{
      for(j in c(1:dim(Log_Deriv)[2])){
        Log_Deriv[1,j]<- as.character(Simplify(Log_Deriv[1,j]))
      }
    }
    Deriv_logjump0 <- function_int2(expr=Log_Deriv,
                                    par=model@parameter@jump,
                                    variable=solve.variable,
                                    drift=FALSE)
    eval(parse(text = Deriv_logjump0))
  }else{
    
  }
  
  model1 <-new("yuimaLocStableSDE",
               drift_func =A_coeff,
               jump_func = Jump_coeff,
               Deriv_drift = Deriv_d,
               Deriv_logjump = Deriv_logjump,
               drift  = model@drift,
               diffusion  = model@diffusion,
               hurst  = model@hurst,
               jump.coeff  = model@jump.coeff,
               measure  = model@measure,
               measure.type  = model@measure.type,
               parameter  = model@parameter,
               state.variable  = model@state.variable,
               jump.variable  = model@jump.variable,
               time.variable  = model@time.variable,
               noise.number  = model@noise.number,
               equation.number  = model@equation.number,
               dimension  = model@dimension,
               solve.variable  = model@solve.variable,
               xinit  = model@xinit,
               J.flag  = model@J.flag)
  
  return(model1)
}


setClass("yuimaStable.Info",
         representation(
           NumericalHess = "logical",
           info = "list",
           parallel = "logical",
           joint = "logical",
           Todorov = "logical",
           N = "numeric",
           num_of_cores = "numeric",
           pos = "numeric",
           posInt = "numeric",
           r0 = "numeric",  
           maxcount = "numeric", 
           scale = "numeric", 
           onlypar = "logical",
           W = "function", 
           Residual = "logical", 
           aa_alt = "list"
         ))

setMethod("initialize", "yuimaStable.Info",
          function(.Object,
                   NumericalHess = FALSE,
                   info = NULL,
                   parallel = T,
                   joint = T,
                   Todorov = F,
                   N = 175,
                   num_of_cores = detectCores()-2,
                   pos = 7,
                   posInt = 1,
                   r0 = 0.9,  
                   maxcount = 100, 
                   scale = 2, 
                   onlypar = FALSE,
                   W = function(x){1/(1+x[1:(length(x)-1)]^2)}, 
                   Residual = T, 
                   aa_alt = gauss.quad(10000,"laguerre",0)
                   ){
            .Object@NumericalHess <- NumericalHess
            .Object@info <- info
            .Object@parallel <- parallel
            .Object@joint <- joint
            .Object@Todorov <- Todorov
            .Object@N <- N
            .Object@num_of_cores <- num_of_cores
            .Object@pos <- pos
            .Object@posInt <- posInt
            .Object@r0 <- r0  
            .Object@maxcount <- maxcount 
            .Object@scale <- scale 
            .Object@onlypar <- onlypar
            .Object@W <- W 
            .Object@Residual <- Residual 
            .Object@aa_alt <- aa_alt
            return(.Object)
          })
setNumerical<- function(NumericalHess = FALSE,
                        info = NULL,
                        parallel = T,
                        joint = T,
                        Todorov = F,
                        N = 175,
                        num_of_cores = max(detectCores()-2,1),
                        pos = 7,
                        posInt = 1,
                        r0 = 0.9,  
                        maxcount = 100, 
                        scale = 2, 
                        onlypar = FALSE,
                        W = function(x){1/(1+x[1:(length(x)-1)]^2)}, 
                        Residual = T, 
                        aa_alt = gauss.quad(10000,"laguerre",0)){
  res<-new("yuimaStable.Info",
           NumericalHess = NumericalHess,
           info = info,
           parallel = parallel,
           joint = joint,
           Todorov = Todorov,
           N = N,
           num_of_cores = num_of_cores,
           pos = pos,
           posInt = posInt,
           r0 = r0,  
           maxcount = maxcount, 
           scale = scale,
           onlypar = onlypar,
           W = W, 
           Residual = Residual, 
           aa_alt = aa_alt)
  return(res)
}