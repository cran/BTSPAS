## 2011-05-15 CJS limited the etaU=log(U) to a maximum of 15 which corresponds to around 400,000,000 fish. 
## 2011-05-09 CJS subtle bug with initial values of epsilon where if fixed values for logitP at the end of the
##                experiment, then the initial values for epsilon must be truncated
## Function to generate initial values for each chain of the MCMC algorithm.
##
## Every model requires different initial values, though much of the
## code can be reused.

genInitsTTln <-
    function(n1,m2,u2){
        ## Generate initial parameters for log-normal travel time model
        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        m2dot1 <- apply(m2[,1:Nstrata.cap],1,sum)

        init.muLogTT <- rep(NA,Nstrata.rel)
        tmp1 <- (m2dot1>0)
        init.muLogTT[tmp1] <- log((m2[tmp1,1:Nstrata.cap] %*% 1:Nstrata.cap) /
                                  (m2dot1[tmp1]) - (1:Nstrata.rel)[tmp1])

        init.muLogTT[!tmp1] <- mean(init.muLogTT[tmp1])

        init.xiMu <- mean(init.muLogTT)
        init.tauMu <- 1/var(init.muLogTT)

        init.etasdLogTT <- log(rep(.5,Nstrata.rel))  # note that log (sd(log travel time)) is being modelled
        init.xiSd <- mean(init.etasdLogTT)
        init.tauSd <- 1

        return(list(muLogTT=init.muLogTT,
                    xiMu=init.xiMu,
                    tauMu=init.tauMu,
                    xiSd=init.xiSd,
                    tauSd=init.tauSd,
                    etasdLogTT=init.etasdLogTT))
    }

genInitsTTnp <-
    function(n1,m2,u2,Delta.max){
        ## Generate initial parameters for non-parametric travel time model

        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        ## Compute empirical theta matrix
        init.Theta <- t(sapply(1:Nstrata.rel,function(i){
            if(all(is.na(m2[i,])) || sum(m2[i,])==0)
                return(rep(NA,Delta.max+1))

            else{
                thetatmp <- pmax(.01,
                                 pmin(m2[i,-(Delta.max+2)]/sum(m2[i,-(Delta.max+2)],na.rm=TRUE),
                                      .99,na.rm=TRUE))  # CJS 2011-02-16
                return(thetatmp/sum(thetatmp))
            }
        }))

        ## Compute initial r
        init.delta <- t(apply(as.matrix(init.Theta[,-(Delta.max+1)]),1,   # CJS 2011-02-16 as.matrix added
                              function(theta){    # CJS fixed -(Delta.max+1)
                                  if(length(theta) == 1){theta}
                                  else {theta/(1-c(0,cumsum(theta[-Delta.max])))}
                              }))

        init.r <- log(init.delta)

        ## mean and standard deviation of transition probabilties
        init.muTT <- apply(logit(init.delta),2,mean,na.rm=TRUE)
        init.sdTT <- sd(as.vector(t(logit(init.delta)))-init.muTT,na.rm=TRUE)

        return(list(muTT=init.muTT,
                    tauTT=1/init.sdTT^2,
                    r=init.r,Theta=init.Theta))
    }

genInitValsChain <-
    function(model,
             n1,                          # Individuals marked per strata at first location
             m2,                          # Individuals recovered at second location
             u2,                          # (List of) unmarked individuals captured per strata with a single spline
             Delta.max=NULL,              # Max travel time for NP model
             logitP.cov,                  # Covariate matrix for capture probabilities
             logitP.fixed=NULL,
             SplineDesign,                # (List of) design matrix(ces) for splines
             hatch.after=NULL,            # Data of release for hatchery fish in model with two splines
             pScale=1){

        ## Generate initial values for a single chain
        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        inits <- list()                   # Create empty list of initial values

        ## 1) Travel time parameters (for non-diagonal models only)
        if(model %in% "TSPNDE"){
            inits <- append(inits,genInitsTTln(n1,m2,u2))
        }
        if(model %in% "TSPNDENP"){
            inits <- append(inits,genInitsTTnp(n1,m2,u2,Delta.max))
        }

        ## 2) Capture probabilities
        ## 2.1) Compute initial logit capture probabilities
        if(model %in% c("TSPDE","TSPDE-WHchinook","TSPDE-WHsteel")){
            init.P <- (m2+1)/(n1+2) * pScale
        }
        else if(model %in% c("TSPNDE")){
            ## Compute expected number of marked fish in each cell
            Theta <- t(sapply(1:Nstrata.rel,function(i){
                tmp <- pnorm(log(i:Nstrata.cap),inits$muLogTT[i],exp(inits$etasdLogTT[i]))
                c(rep(0,(i-1)),tmp - c(0,tmp[-(Nstrata.cap - (i-1))]))
            }))

            M <- Theta * n1

            m2dot2 <- apply(m2[,1:Nstrata.cap],2,sum)
            init.P <- (m2dot2 + 1)/(apply(M,2,sum) + m2dot2 + 1) * pScale
        }
        else if(model %in% c("TSPNDENP")){
            ## Compute expected number of marked fish in each cell
            N2 <- lapply(1:Nstrata.rel,function(i) inits$Theta[i,]*n1[i])

            ## Compute expected number of marked fish in each capture strata
            n2 <- sapply(1:Nstrata.cap,function(i){
                n2tmp <- 0
                for(j in max(i-Delta.max,1):min(i,Nstrata.rel))
                    n2tmp <- N2[[j]][i-j+1] + n2tmp
                n2tmp
            })

            m2dot2 <- sapply(1:Nstrata.cap,function(i){
                m2tmp <- 0
                for(j in max(i-Delta.max,1):min(i,Nstrata.rel))
                    m2tmp <- m2[j,i-j+1] + m2tmp
                m2tmp
            })

            init.P <- (m2dot2 + 1)/(n2 + m2dot2 + 1) * pScale

        }

                                        # Constrain p to the interval (.00001,.99999)
        init.P <- pmax(.00001,pmin(init.P,.99999))

                                        # Remove missing values
        init.P[is.na(init.P)] <- mean(init.P,na.rm=TRUE)

                                        # Remove fixed values from inital vector
        init.P[!is.na(logitP.fixed)] <- NA

                                        # Compute logit
        init.logitP <- logit(init.P)

        ## 2.2) Compute associated coefficients for design matrix
        init.beta.logitP <- as.vector(lm(init.logitP ~ logitP.cov - 1)$coeff)

        ## 2.3) Set variance for hierarchical model of capture probabilities
        if(length(init.beta.logitP)==1)
            init.tauP <- 1/var(init.logitP - logitP.cov*init.beta.logitP,na.rm=TRUE)
        else
            init.tauP <- 1/var(init.logitP - logitP.cov %*% init.beta.logitP,na.rm=TRUE)

        init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a vector

        init.beta.logitP[is.na(init.beta.logitP)] <- 0

        inits <- append(inits,list(beta.logitP=init.beta.logitP,tauP=as.numeric(init.tauP)))

        ## 3) Numbers of unmarked individuals per strata (where u2 observed)
        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            init.U <- ceiling((u2+1)/init.P)
        }

        ## Option 2: Chinook model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHchinook")){
            init.U.W <- ceiling((u2$W+1)/init.P)
            init.U.H <- ceiling((u2$H+1)/init.P)
            init.U.H[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery
        }

        ## Option 3: Steelhead model with separate splines for wild, wild YOY, and hatchery fish
        if(model %in% c("TSPDE-WHsteel")){
            init.U.W.YoY <- ceiling((u2$W.YoY+1)/init.P)
            init.U.W.1 <- ceiling((u2$W.1+1)/init.P)
            init.U.H.1 <- ceiling((u2$H.1+1)/init.P)
            init.U.H.1[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery
        }

        ## 4) Spline coefficients

        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            ## 4.1) Fit Spline to strata with u2 observed
            tmp1 <- !is.na(init.U)
            init.bU <- lm(log(init.U[tmp1]) ~ SplineDesign[tmp1,]-1)$coeff
            init.bU[is.na(init.bU)] <- mean(init.bU,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.2) Compute variance of second differences between coefficients
            tmp2 <- 3:length(init.bU)
            sigmaU <- sd(init.bU[tmp2]-2*init.bU[tmp2-1]+init.bU[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU=init.bU,tauU=init.tauU))
        }

        ## Option 2: Chinook model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHchinook")){
            ## 4.1.a) Fit spline to wild fish
            tmp1.W <- !is.na(init.U.W)
            init.bU.W <- lm(log(init.U.W[tmp1.W]) ~ SplineDesign$W[tmp1.W,]-1)$coeff
            init.bU.W[is.na(init.bU.W)] <- mean(init.bU.W,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.b) Fit spline to hatchery fish
            tmp1.H <- c(rep(FALSE,hatch.after),!is.na(init.U.H[-(1:hatch.after)]))
            init.bU.H <- lm(log(init.U.H[tmp1.H]) ~ SplineDesign$H[tmp1.H,]-1)$coeff
            init.bU.H[is.na(init.bU.H)] <- mean(init.bU.H,na.rm=TRUE)       # Fix any coefficients that can't be c

            ## 4.2) Variance of second differences between coefficients (use only wild fish to initialize)
            tmp2 <- 3:length(init.bU.W)
            sigmaU <- sd(init.bU.W[tmp2]-2*init.bU.W[tmp2-1]+init.bU.W[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU.W=init.bU.W,bU.H=init.bU.H,tauU=init.tauU))
        }

        ## Option 3: Steelhead model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHsteel")){
            ## 4.1.a) Fit spline to wild YoY fish
            tmp1.W.YoY <- !is.na(init.U.W.YoY)
            init.bU.W.YoY <- lm(log(init.U.W.YoY[tmp1.W.YoY]) ~ SplineDesign$W.YoY[tmp1.W.YoY,]-1)$coeff
            init.bU.W.YoY[is.na(init.bU.W.YoY)] <- mean(init.bU.W.YoY,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.b) Fit spline to wild 1 fish
            tmp1.W.1 <- !is.na(init.U.W.1)
            init.bU.W.1 <- lm(log(init.U.W.1[tmp1.W.1]) ~ SplineDesign$W.1[tmp1.W.1,]-1)$coeff
            init.bU.W.1[is.na(init.bU.W.1)] <- mean(init.bU.W.1,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.c) Fit spline to hatchery fish
            tmp1.H.1 <- c(rep(FALSE,hatch.after),!is.na(init.U.H.1[-(1:hatch.after)]))
            init.bU.H.1 <- lm(log(init.U.H.1[tmp1.H.1]) ~ SplineDesign$H.1[tmp1.H.1,]-1)$coeff
            init.bU.H.1[is.na(init.bU.H.1)] <- mean(init.bU.H.1,na.rm=TRUE)       # Fix any coefficients that can't be c

            ## 4.2) Variance of second differences between coefficients (use only wild YoY fish to initialize)
            tmp2 <- 3:length(init.bU.W.YoY)
            sigmaU <- sd(init.bU.W.YoY[tmp2]-2*init.bU.W.YoY[tmp2-1]+init.bU.W.YoY[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU.W.YoY=init.bU.W.YoY,bU.W.1=init.bU.W.1,
                                       bU.H.1=init.bU.H.1,tauU=init.tauU))
        }

        ## 5) Variance about spline
        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            sigmaeU <- sd(log(init.U+1) - SplineDesign %*% init.bU,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        ## Option 2: Chinook models with two splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHchinook")){
            sigmaeU <- sd(log(init.U.W+1) - SplineDesign$W %*% init.bU.W,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        ## Option 3: Steelhead models with three splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHsteel")){
            sigmaeU <- sd(log(init.U.W.YoY+1) - SplineDesign$W.YoY %*% init.bU.W.YoY,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        inits <- append(inits,list(taueU=init.taueU))

        ## 6) Initialize missing U values by fitting spline and generating errors
        ## Option 1: Models with only 1 spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            if(sum(!tmp1)>0)
                init.U[!tmp1] <- ceiling(exp(as.vector(SplineDesign[!tmp1,] %*% init.bU)
                                             + rnorm(sum(!tmp1),0,sigmaeU))) + 1
            init.etaU <- pmin(log(init.U),20)  # limit the initial values to reasonable values

            inits <- append(inits,list(etaU=init.etaU))
        }

        ## Option 2: Chinook models with two splines
        if(model %in% c("TSPDE-WHchinook")){
            ## Wild fish
            if(sum(!tmp1.W)>0)
                init.U[!tmp1.W] <- ceiling(exp(as.vector(SplineDesign$W[!tmp1.W,] %*% init.bU.W)
                                               + rnorm(sum(!tmp1.W),0,sigmaeU))) + 1
            init.etaU.W <- pmin(log(init.U.W), 15) # limit the initial values to reasonable values

            ## Hatchery fish
            tmp2.H <- tmp1.H[-(1:hatch.after)]
            if(sum(!tmp2.H)>0)
                init.U[!tmp2.H] <- ceiling(exp(as.vector(SplineDesign$H[!tmp2.H,] %*% init.bU.H)
                                               + rnorm(sum(!tmp2.H),0,sigmaeU))) + 1

            init.etaU.H <- c(rep(NA,hatch.after),pmin(20,log(init.U.H[tmp1.H])))
      

            inits <- append(inits,list(etaU.W=init.etaU.W,etaU.H=init.etaU.H))
        }

        ## Option 3: Steelhead models with three splines
        if(model %in% c("TSPDE-WHsteel")){
            ## Wild YoY fish
            if(sum(!tmp1.W.YoY)>0)
                init.U[!tmp1.W.YoY] <- ceiling(exp(as.vector(SplineDesign$W.YoY[!tmp1.W.YoY,] %*% init.bU.W.YoY)
                                                   + rnorm(sum(!tmp1.W.YoY),0,sigmaeU))) + 1
            init.etaU.W.YoY <- pmin(20,log(init.U.W.YoY))


            ## Wild 1 fish
            if(sum(!tmp1.W.1)>0)
                init.U[!tmp1.W.1] <- ceiling(exp(as.vector(SplineDesign$W.1[!tmp1.W.1,] %*% init.bU.W.1)
                                                 + rnorm(sum(!tmp1.W.1),0,sigmaeU))) + 1
            init.etaU.W.1 <- pmin(20,log(init.U.W.1))

            ## Hatchery 1 fish
            tmp2.H.1 <- tmp1.H.1[-(1:hatch.after)]
            if(sum(!tmp2.H.1)>0)
                init.U[!tmp2.H.1] <- ceiling(exp(as.vector(SplineDesign$H.1[!tmp2.H.1,] %*% init.bU.H.1)
                                                 + rnorm(sum(!tmp2.H.1),0,sigmaeU))) + 1

            init.etaU.H.1 <- c(rep(NA,hatch.after),pmin(20,log(init.U.H.1[tmp1.H.1])))

            inits <- append(inits,list(etaU.W.YoY=init.etaU.W.YoY,
                                       etaU.W.1=init.etaU.W.1,etaU.H.1=init.etaU.H.1))
        }

        ## 7) Transform initial values for logitP to initial values for epsilon
        ## Option 1: Models with only one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            init.epsilon <- init.logitP - log(u2 + 1) + inits$etaU
            # subtle problem. If the logitP.fixed include elements at the end of the experiment
            # then init.epsion needs to be truncated at the end, otherwise OPENBugs gets upset
            if(length(logitP.fixed)>0){
               for(i in length(logitP.fixed):1){
                  if( is.na(logitP.fixed[i])){break}
                  init.epsilon <- init.epsilon[-length(init.epsilon)] # drop last term
               }
            }
        }

        ## Option 2: Chinook models with two splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHchinook")){
            init.epsilon <- init.logitP - log(u2$W + 1) + inits$etaU.W
        }

        ## Option 3: Steelhead models with three splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHsteel")){
            init.epsilon <- init.logitP - log(u2$W.YoY + 1) + inits$etaU.W.YoY
        }

        inits <- append(inits,list(epsilon=c(init.epsilon)))

        ## Remove working objects from the initial values
        if(model=="TSPNDENP"){
            inits$Theta <- NULL
        }

        return(inits)

        return(inits)
    }

genInitVals <-
    function(model,
             n1,                          # Individuals marked per strata at first location
             m2,                          # Individuals recovered at second location
             u2=NULL,                     # (List of) unmarked individuals captured
             Delta.max=NULL,              # Max travel time for NP model
             logitP.cov,                  # Covariate matrix for capture probabilities
             logitP.fixed=NULL,
             SplineDesign,                # (List of) desgin matrices for spline for models
             hatch.after=NULL,            # Data of release for hatchery fish in model with two splines
             n.chains=3,
             pStep=5){                   # Relative change in p between chains

        ## Generate initial values for n.chains

        lapply(1:n.chains,function(i){
            pScale <- pStep ^(-(n.chains-1)/2 + (i-1))

            genInitValsChain(model,
                             n1,m2,u2,Delta.max,
                             logitP.cov,logitP.fixed,
                             SplineDesign,hatch.after,
                             pScale)
        })
    }
