################################
## Plasticity and bet-hedging ##
###### T. R. Haaland 2020 ######
################################

rm(list=ls())
library(RColorBrewer)


##cost(): Calculate costs of plasticity
## Arguments:
##   x: Reaction norm slope (for c_s) or magnitude of phenotypic update (for c_u)
##   steepness: Slope and exponent parameters (s_1,s_2 or u_1,u_2)
cost <- function(x,steepness=c(1,1)){
  return(steepness[1]*x^steepness[2])
}

## sim(): Individual-based simulation.
## Arguments:
##   n:      number of time steps/instances before reproduction, over which to gather resources.
##   alpha:  Between-year mortality (only used if type="Overlap")
##   m.rate: Mutation  rate. Default 0.005
##   m.size: Sd of the normal distribution around the old gene value from which new gene value is drawn.
##   K:      Carrying capacity. Default 5000. Density dependence is 'ceiling' type.
##   T:      Number of years/generations.
##   reps:   Number of replicate simulations. Default 1
##   plot:   Logical. Whether to plot population trajectories, or not (default not).
##   grain_env: Correlation among individual in environments experienced in a time step. [0,1]
##   r:      Intrinsic rate of increase Default 2.
##   costvar:How are costs of plasticity implemented? "slope" (only c_s, default), "update" (c_u) or "both"
##   steepness: Cost function parameters s_1,s_2,u_1,u_2. See cost() function above.
##   Title:  Plot title. Default blank
##   decayvar:How does fitness decay with mismatch? "gaussian" (default), "exponential" or "linear"
##   decay:  How steeply does fitness decay with mismatch? Only used if decayvar="exponential".
## Output: A list of two elements:
##   [[1]] nstorage: Matrix of population sizes for each replicate (row) at each time step (column)
##   [[2]] zstorage: Matrix of mean RN slope for each replicate (row) at each time step (column)
sim <- function(n,alpha=0.2,m.rate=0.005,m.size=0.05,K=5000,T=2000,reps=1,plot=FALSE,r=2,
                grain_env=0.1,costvar="slope",title="",steepness=c(1,1,1,1),decayvar="gaussian",decay=5){
  if(plot){
    plot(1:T,rep(1,T),type="n",ylim=c(0,1),ylab="Reaction norm slope",xlab="Year",main=title)
  }
  
  nstorage <- zstorage <- matrix(NA,reps,T) #Storage matrices for population size and phenotypes.
  
  for(rep in 1:reps){
    #Initiate population
    N <- K #Starting pop. size is at carrying capacity
    pop <- runif(N) #Gene for reaction norm slope. Initialize uniform.
    
    for(Time in 1:T){
      #Record population traits
      zstorage[rep,Time] <- mean(abs(pop)) #Record plasticities
      nstorage[rep,Time] <- N #Record population size.
      
      W <- rep(0,N)
      
      # Developmental environment determining starting phenotypes (only necessary if including updatecosts)
      if(costvar=="update"||costvar=="both"){
        env <- runif(1)
        lower <- grain_env*env
        upper <- env + (1-grain_env)*(1-env)
        envs <- runif(N,lower,upper)*2-1 #Environments are between -1 and 1
        phens <- pop*envs # Phenotype is Environment*RN Slope. RN Slope==1 -> phenotype==environment
      }
      
      for(Step in 1:n){
        if(costvar=="update"||costvar=="both"){prev.phens <- phens} # Record last time's phenotype (only necessary if including updatecosts)
        this.W <- rep(0,N)
        #Find environment(s)
        env <- runif(1)
        lower <- grain_env*env
        upper <- env + (1-grain_env)*(1-env)
        envs <- runif(N,lower,upper)*2-1 #Environments are between -1 and 1
        #Find individual phenotypes
        phens <- pop*envs # Phenotype is Environment*RN Slope. RN Slope==1 -> phenotype==environment
        #Calculate fitnesses dependent on phenotype-environment mismatch (states-conds)
        this.W <- switch(decayvar,
                         "exponential"= dexp(abs(phens-envs),decay)/decay,
                         "linear"= 1-abs(phens-envs),
                         "gaussian"= dnorm(phens-envs,0,0.4))
        W <- W + this.W #Sum fitnesses accumulated
        if(costvar=="update"||costvar=="both"){ #If cost depends on how much the phenotype changed
          W <- W - cost(abs(phens-prev.phens),steepness=steepness[3:4]) 
        }else if(costvar=="slope"||costvar=="both"){ #If cost depends on plasticity gene
          W <- W - cost(abs(pop),steepness=steepness[1:2]) 
        }
      }
      W <- r*(W/n) # Standardize fitness on n, multiply by intrinsic rate of increase
      W[W<0] <- 0 
      
      if(round(sum(W)==0)){
        print(c("Population extinct at time",Time))
        break
      }
      
      #Selection
      if(alpha==1){ # If discrete generations (all adults die)
        NextN <- ifelse(sum(W)>K,K,round(sum(W))) #  Find population size next year
        pop <- sample(pop,size=NextN,prob=W,replace=TRUE) # Choose offspring according to parent's fitness
        if(costvar=="update"||costvar=="both") phens <- rep(0,NextN) # All offspring are adapted to 'mean' environmental conditions.
        
        #Mutation
        mut <- which(runif(NextN)<m.rate) #Select mutated individuals
        pop[mut] <- rnorm(length(mut),pop[mut],m.size)
      }
      else{ # If overlapping generations (some adults survive)
        Alive <- runif(N)>alpha # Random survival
        Alive <- which(Alive==TRUE) # Alive is a sequence of numbers of the individuals who survived.
        
        #How many 'slots' should be filled (i.e. how many offspring should be produced?)
        Offspring <- min(round(sum(W)),K-length(Alive))
        
        #Reproduction. IDs of those who get offspring.
        if(sum(W[Alive])>0){ #In this scenario, offspring die if no surviving parents.
          IDs <- sample(Alive,size=Offspring,prob=W[Alive],replace=TRUE)
        } else {
          print(c("Extinct at time",Time))
          break
        }

        #Next generation is a combination of new and old individuals.
        pop <- c(pop[IDs],pop[Alive])  #IDs has length Offspring, so the Offspring first entries here are new.
        if(costvar[1]=="update") phens <- c(rep(0,Offspring),phens[Alive]) # All offspring are adapted to 'mean' env. conditions.
        #Mutation: Some of the offspring mutate
        mut <- which(runif(Offspring)<m.rate)
        pop[mut] <- rnorm(length(mut),pop[mut],m.size)
        
        NextN <- length(pop)
      }
      N <- NextN
    }
    #print(rep) #Ticker
    if(plot){
      lines(1:T,zstorage[rep,],col=rgb(1,0,0,0.1))
      lines(1:T,nstorage[rep,]/K,col=rgb(0,0,1,0.1))
    }
  }
  return(list(zstorage,nstorage))
}


#### Run scenarios ####
gs <- seq(0.4,1,by=0.1) #Environmental grain scenarios
ns <- c(1,2,5,10) #Number of time steps prior to reproduction
T <- 2000
K <- 5000

#Initialize storage
meanphens <- sdphens <- n <- sdn <- numeric(length(ns))
out <- array(0,dim=c(length(ns),4,length(gs)))
dimnames(out) <- list(ns,c("meanphens","sdphens","n","sdn"),gs)

# Not hard-coded - change the parameters you want in function call
for(g in gs){
  for(i in ns){
    #Run simulation
    test <- sim(n=i,alpha=0.5,reps=20,plot=TRUE,grain_env=g,decayvar="gaussian",costvar=c("update","linear"),steepness=1)
    #Record summary statistics
    repmeans <- apply(test[[1]][,(T-1000):T],2,mean,na.rm=TRUE) #Phenotypes in the last 1000 generations.
    meanphens[which(i==ns)] <- mean(repmeans,na.rm=TRUE) #Mean phenotypes across replicates
    sdphens[which(i==ns)] <- sd(repmeans,na.rm=TRUE) #Sd phenotypes across replicates
    rep_mean_n <- apply(test[[2]][,(T-1000):T],1,mean,na.rm=TRUE) #Population sizes
    n[which(i==ns)] <- mean(rep_mean_n,na.rm=TRUE)/K #Mean population size across replicates (Relative to K)
    sdn[which(i==ns)] <- sd(rep_mean_n,na.rm=TRUE) #Sd population size across replicates
    print(c(g,i)) # Ticker
  }
  out[,,which(g==gs)] <- c(meanphens,sdphens,n,sdn)
}
#Save output
saveRDS(out,file="Your directory and file name.Rdata")

#### Create results plots (Fig. 2 and S2) ####

#Read output
data <- readRDS(file="Your directory and file name.Rdata")

# Set diverging red-blue color scheme
pal <- brewer.pal(length(data[,1,1])+2,"RdBu") 
pal <- c(head(pal,length(pal)/2-1),tail(pal,length(pal)/2-1))
pchs <- c(1,2,4,5) # Point types for plotting - length(pchs) should be equal to length (ns)

plot(1:length(gs),data[1,1,1:length(gs)],ylim=c(0,0.9),type="n",xaxt="n",xlab="Environmental grain (g)",ylab="Reaction norm slope (beta)")
legend("bottomleft",legend=ns,col=pal,lty=1,title="Time steps prior to reproduction (n):",bty="n",pch=pchs,cex=0.9,pt.cex=1.1)
axis(side=1,at=1:length(gs),label=gs)
for(i in 1:length(ns)){
  points(1:length(gs),data[i,1,1:length(gs)],pch=pchs[i],cex=1,col=pal[i])
  lines(1:length(gs),data[i,1,1:length(gs)],col=pal[i])
  for(j in 1:length(gs)){
    arrows(j,data[i,1,j],j,data[i,1,j]+data[i,2,j],angle=90,length=0.1,col=pal[i])
    arrows(j,data[i,1,j],j,data[i,1,j]-data[i,2,j],angle=90,length=0.1,col=pal[i])
  }
}

###################################
## Calculate variance in fitness ##
###################################

## Gaussian W(m)
# b is the different plasticity slopes (should be a vector, e.g. seq(0,1,by=0.01))
# cost is whether costs depend on RN slope ("slope", default) or how much the phenotype is updated ("update") or "both"
# steepness[1:2] determines linear and exponent of slopecosts (if cost= "slope" or "both") - default 1, 2
# steepness[3:4] determines linear and exponent of updatecosts (if cost= "update" or "both") - default 0.5, 1
# Returns a vector of:
# [1] Variance in fitness across environments
# [2] Arithmetic mean fitness across environments
# [3] Geometric mean fitness across environments (approximation)
# [4] Geometric mean fitness across environments (exact - but doesn't work if some fitnesses are negative)
sim_gaus <- function(b,cost=c("slope"),steepness=c(1,2,0.5,1)){
  es <- runif(1000000,-1,1) # Simulate random environments
  zs <- b*es # Phenotypes produced in each env.
  ws <- dnorm(tail(abs(zs-es),length(zs)-1),0,0.4) # Fitness benefits - first env. is "developmental env." (in case of using updatecosts)
  slopecost <- cost(b,steepness=steepness[1:2]) # c_s
  updatecost <- cost(abs(head(zs,length(zs)-1)-tail(zs,length(zs)-1)),steepness=steepness[3:4]) # c_u
  
  if(cost=="slope"){
    ws <- ws-slopecost
  } else if(cost=="update"){
    ws <- ws-updatecost
  } else if(cost=="both"){
    ws <- ws-slopecost-updatecost
  }

  logws <- log(ws)
  logws[which(is.nan(logws)==TRUE)] <- -10000
  return(c(var(ws),mean(ws),mean(ws)-(2*var(ws))/mean(ws),exp(mean(logws))))
}

#### Create predictions figures (S4 and S5) ####
bs <- seq(0,1,by=0.01) #For any RN slope
vars <- ameans <- gmeansapprox <- gmeansexact <- numeric(length(bs))
for(i in 1:length(bs)){
  tmp <- sim_gaus(bs[i],cost="update",steepness=c(0,0,1,2))
  vars[i] <- tmp[1]
  ameans[i] <- tmp[2]
  gmeansapprox[i] <- tmp[3]
  gmeansexact[i] <- tmp[4]
  print(i)
}
plot(bs,gmeansapprox,main=expression("A) No "*italic(c)[slope]*", Quadratic "*italic(c)[update]),type="l",col="Red",xlab="",ylim=c(0,1),ylab="Fitness",cex.lab=1.3,lwd=2)
lines(bs,ameans,lwd=2,col=rgb(0,0,1,0.6))
lines(bs,sqrt(vars))
arrows(bs[which.max(ameans)],-0.05,bs[which.max(ameans)],max(ameans),length=0,col=rgb(0,0,1,0.6),lty=2)
arrows(bs[which.max(gmeansapprox)],-0.05,bs[which.max(gmeansapprox)],max(gmeansapprox),length=0,col=" Red",lty=2)
legend("topleft",legend=c("Arithmetic mean","Geometric mean","Standard deviation"),
       lwd=c(2,2,1),lty=c(1,1,1),col=c(rgb(0,0,1,0.6),"Red","Black"),bty="n")