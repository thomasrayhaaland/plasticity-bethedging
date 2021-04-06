################################
## Plasticity and bet-hedging ##
#### (C) T. R. Haaland 2019 ####
################################

#### Preamble ####
rm(list=ls())
library(RColorBrewer)
setwd("C:\\Users\\Thomas\\Documents\\Plasticity and bet-hedging\\Output\\")
# Set directory: Use find-and-replace-all on 'C:\\Users\\Thomas\\Documents\\Plasticity and bet-hedging\\'.

##cost(): Calculate costs of plasticity
## Arguments:
##   x: Reaction norm slope (for c_m) or magnitude of phenotypic update (for c_p)
##   steepness: Slope and exponent parameters (m_1,m_2 or p_1,p_2)
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
##   r:     Maximum reproductive rate. Default 5.
##   steepness: Cost function parameters, a vector of format c(m_1, m_2, p_1, p_2. See cost() function above.
##   Title:  Plot title. Default blank
##   decayvar:How does fitness decay with mismatch? "gaussian" (default), "exponential" or "linear"
##   decay:  How steeply does fitness decay with mismatch? Only used if decayvar="exponential".
## Output: A list of two elements:
##   [[1]] nstorage: Matrix of population sizes for each replicate (row) at each time step (column)
##   [[2]] zstorage: Matrix of mean RN slope for each replicate (row) at each time step (column)
sim <- function(n,alpha=0.2,m.rate=0.005,m.size=0.05,K=5000,T=2000,reps=1,plot=FALSE,r=5,
                grain_env=0.5,title="",decay=5,steepness=c(1,1,1,1),decayvar="gaussian"){

  if(plot){
    plot(1:T,rep(1,T),type="n",ylim=c(0,1),ylab="Reaction norm slope",xlab="Year",main=title)
  }
   
  nstorage <- zstorage <- matrix(NA,reps,T) #Storage matrices for population size and phenotypes.
  
  costvar <- ifelse(any(steepness[1:2]>0),"maintenance","production") # Determine the plasticity cost structure
  costvar <- ifelse(any(steepness[1:2]>0)&&any(steepness[3:4]>0),"both",costvar)
  
  for(rep in 1:reps){
    #Initiate population
    N <- K #Starting pop. size is at carrying capacity
    pop <- runif(N) #Gene for reaction norm slope. Initialize uniform.
    
    for(Time in 1:T){
      #Record population traits
      zstorage[rep,Time] <- mean(abs(pop)) #Record plasticities
      nstorage[rep,Time] <- N #Record population size.
      
      W <- rep(0,N)
      
      # Developmental environment determining starting phenotypes (only necessary if including production costs)
      if(costvar=="production"||costvar=="both"){
        env <- runif(1,-1,1)
        lower <- env - (1-grain_env)*(1+env)
        upper <- env + (1-grain_env)*(1-env)
        envs <- runif(N,lower,upper) 
        phens <- pop*envs # Phenotype is Environment*RN Slope. RN Slope==1 -> phenotype==environment
      }
      
      for(Step in 1:n){
        if(costvar=="production"||costvar=="both"){prev.phens <- phens} # Record last time's phenotype
        this.W <- rep(0,N)
        #Find environment(s)
        env <- runif(1,-1,1)
        lower <- env - (1-grain_env)*(1+env)
        upper <- env + (1-grain_env)*(1-env)
        envs <- runif(N,lower,upper) 
        #Find individual phenotypes
        phens <- pop*envs # Phenotype is Reaction norm Slope * Environment.
        #Calculate fitnesses dependent on phenotype-environment mismatch
        this.W <- switch(decayvar, # This is called V_t in the paper.
                         "exponential"= dexp(abs(phens-envs),decay)/decay,
                         "linear"= 1-abs(phens-envs),
                         "gaussian"= dnorm(phens-envs,0,0.4)) # This is the one used in the paper. 
        W <- W + this.W # Sum accumulated payoffs.
        if(costvar=="production"||costvar=="both"){ #If cost depends on how much the phenotype changed
          W <- W - cost(abs(phens-prev.phens),steepness=steepness[3:4]) 
        }else if(costvar=="maintenance"||costvar=="both"){ #If cost depends on plasticity gene
          W <- W - cost(abs(pop),steepness=steepness[1:2]) 
        }
      }
      W <- r*(W/n) # Standardize fitness on n, multiply by intrinsic rate of increase
      W[W<0] <- 0
      
      if(round(sum(W)==0)){
        print(c("Extinct at time",Time))
        break
      }
      
      #Selection
      if(alpha==1){ # If discrete generations (all adults die)
        NextN <- ifelse(sum(W)>K,K,round(sum(W))) # How many individuals will the next generation consist of?
        pop <- sample(pop,size=NextN,prob=W,replace=TRUE) # Draw offspring depending on parents' relative fitness.
        if(costvar=="production"||costvar=="both") phens <- rep(0,NextN) # Next generation starts with phenotypes 0 (only relevant if using production costs)
        
        #Mutation
        mut <- which(runif(NextN)<m.rate) #Select mutated individuals
        pop[mut] <- rnorm(length(mut),pop[mut],m.size)
      }
      else{ # If overlapping generations (some adults survive)
        lucky <- runif(N)>alpha
        Alive <- which(lucky==TRUE) # Alive is a sequence of numbers of the individuals who survived.
        
        #How many 'slots' should be filled (i.e. how many offspring should be recruited?)
        Offspring <- min(round(sum(W)),K-length(Alive))
        
        #Reproduction. IDs of those who get offspring. If mortality happens before breeding:
        if(sum(W[Alive])>0){ #In this scenario, offspring die if no surviving parents.
          if(length(Alive)>1){ #To avoid a syntax error in sample().
            IDs <- sample(Alive,size=Offspring,prob=W[Alive],replace=TRUE)
          } else {
            IDs <- rep(Alive,Offspring)
          }
        } else {
          print(c("Extinct at time",Time))
          break
        }
        #If mortality happens after breeding (hashtag away if using):
        #IDs <- sample(1:N,size=Offspring,prob=W,replace=TRUE)
        
        #Next generation is a combination of new and old individuals.
        pop <- c(pop[IDs],pop[Alive])  #IDs has length Offspring, so the Offspring first pops are new.
        if(costvar=="update") phens <- c(rep(0,Offspring),phens[Alive])
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
#####

## Generate the data shown in the paper ####
gs <- seq(0.5,1,by=0.1)
ns <- c(1,2,5,10)
set.seed(1)
for(gg in gs){ 
  for(nn in ns){
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(1,2,0,0))
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(1,2,0,0))
    saveRDS(overlap, paste0("steepness=c(1,2,0,0), alpha=0.5, n=",nn,", g=",gg,".R")) #2A
    saveRDS(discrete, paste0("steepness=c(1,2,0,0), alpha=1, n=",nn,", g=",gg,".R")) #S2A
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0,0,1,2))
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0,0,1,2))
    saveRDS(overlap, paste0("steepness=c(0,0,1,2), alpha=0.5, n=",nn,", g=",gg,".R")) #2B
    saveRDS(discrete, paste0("steepness=c(0,0,1,2), alpha=1, n=",nn,", g=",gg,".R")) #S2B
    print(nn)
  }
}

for(gg in gs){
  for(nn in ns){
    tmp <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0.8,1,0,0))
    saveRDS(tmp, paste0("steepness=c(0.8,1,0,0), alpha=0.5, n=",nn,", g=",gg,".R")) #4A
    tmp <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0,0,0.8,1),r=5)
    saveRDS(tmp, paste0("steepness=c(0,0,0.8,1), alpha=0.5, n=",nn,", g=",gg,".R")) #4B
    tmp <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0.8,1,0,0))
    saveRDS(tmp, paste0("steepness=c(0.8,1,0,0), alpha=1, n=",nn,", g=",gg,".R")) #4A
    tmp <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0,0,0.8,1),r=5)
    saveRDS(tmp, paste0("steepness=c(0,0,0.8,1), alpha=1, n=",nn,", g=",gg,".R")) #4B
    print(nn)
  }
}

for(gg in gs){ 
  for(nn in ns){
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0.4,1,0.4,1))
    saveRDS(discrete, paste0("steepness=c(0.4,1,0.4,1), alpha=1, n=",nn,", g=",gg,".R")) #5A
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0.4,2,0.4,1))
    saveRDS(discrete, paste0("steepness=c(0.4,2,0.4,1), alpha=1, n=",nn,", g=",gg,".R")) #5B
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0.4,1,0.4,1))
    saveRDS(overlap, paste0("steepness=c(0.4,1,0.4,1), alpha=0.5, n=",nn,", g=",gg,".R")) #5A
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0.4,2,0.4,1))
    saveRDS(overlap, paste0("steepness=c(0.4,2,0.4,1), alpha=0.5, n=",nn,", g=",gg,".R")) #5B
    print(nn)
  }
}
for(gg in gs){ 
  for(nn in ns){
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0.4,1,0.4,2))
    saveRDS(discrete, paste0("steepness=c(0.4,1,0.4,2), alpha=1, n=",nn,", g=",gg,".R")) # 5C
    discrete <- sim(n=nn,grain_env=gg,alpha=1,reps=100,plot=T,steepness=c(0.4,2,0.4,2))
    saveRDS(discrete, paste0("steepness=c(0.4,2,0.4,2), alpha=1, n=",nn,", g=",gg,".R")) # 5D
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0.4,1,0.4,2))
    saveRDS(overlap, paste0("steepness=c(0.4,1,0.4,2), alpha=0.5, n=",nn,", g=",gg,".R")) # 5C
    overlap <- sim(n=nn,grain_env=gg,alpha=0.5,reps=100,plot=T,steepness=c(0.4,2,0.4,2))
    saveRDS(overlap, paste0("steepness=c(0.4,2,0.4,2), alpha=0.5, n=",nn,", g=",gg,".R")) # 5D
    
    print(nn)
  }
}
#####

## Retrieve data ####

## Function to obtain summary stats from raw files ##
# Provide scenario with file name in format "steepness=c(m1,m2,p1,p2), alpha=a")
stats <- function(scenario,gs=seq(0.5,1,by=0.1),ns=c(1,2,5,10)){
  #Initialize storage
  mean_phens <- sd_phens <- mean_n <- sd_n <- numeric(length(ns))
  out <- array(0,dim=c(length(ns),4,length(gs)))
  dimnames(out) <- list(ns,c("mean_phens","sd_phens","mean_n","sd_n"),gs)
  
  for(g in gs){
    for(n in ns){
      #Retrieve file
      data <- readRDS(paste0(scenario,", n=",n,", g=",g,".R"))
      zdat <- data[[1]]
      ndat <- data[[2]]
      T <- dim(data[[1]])[2]
      reps <- dim(data[[1]])[1]
      K <- data[[2]][1,2]
      
      #Record summary statistics
      repmeans <- apply(zdat[,(T/2):T],1,mean,na.rm=TRUE) #Phenotypes in the last 50% of generations.
      mean_phens[which(ns==n)] <- mean(repmeans,na.rm=TRUE) #Mean phenotypes across replicates
      sd_phens[which(ns==n)] <- sd(repmeans,na.rm=TRUE) #Sd phenotypes across replicates
      rep_mean_n <- apply(ndat[,(T/2):T],1,mean,na.rm=TRUE) #Population sizes
      mean_n[which(ns==n)] <- mean(rep_mean_n,na.rm=TRUE)/5000 #Mean population size across replicates
      sd_n[which(ns==n)] <- sd(rep_mean_n,na.rm=TRUE) #Sd population size across replicates
    }
    out[,,which(gs==g)] <- c(mean_phens,sd_phens,mean_n,sd_n)
  }
  return(out)
}

#
saveRDS(stats("steepness=c(0.4,2,0.4,1), alpha=0.5",gs=gs),file="Summary steepness=c(0.4,2,0.4,1), alpha=0.5.R")
saveRDS(stats("steepness=c(0.4,2,0.4,1), alpha=1",gs=gs),file="Summary steepness=c(0.4,2,0.4,1), alpha=1.R")


## Plot evolutionary trajectories (as in Fig. S1) ##
# Provide scenario in format "steepness=c(m1,m2,p1,p2), alpha=a"
timeplots <- function(scenario,n,g,main="",axes=c("Bottom","Left")){
  data <- readRDS(paste0(scenario,", n=",n,", g=",g,".R"))
  zdat <- data[[1]]
  ndat <- data[[2]]
  T <- dim(data[[1]])[2]
  reps <- dim(data[[1]])[1]
  K <- data[[2]][1,2]
  
  plot(1:T,zdat[1,],type="n",ylim=c(0,1),main=main,xaxt="n",yaxt="n")
  for(i in sample(reps,size=min(20,reps))){
    lines(1:T,zdat[i,],col=rgb(1,0,0,0.1))
    lines(1:T,ndat[i,]/K,col=rgb(0,0,1,0.05))
  }
  if(any(axes=="Left")){
    axis(2,at=c(0,0.4,0.8))
    mtext(c("N/K","RN slope"),col=c("Blue","Red"),side=2,line=c(1.8,2.7),cex=c(0.8,0.9))
  }
  if(any(axes=="Bottom")){
    axis(1,at=seq(0,T,length=5))
    mtext("Year",side=1,line=2.5)
  }
  if(any(axes=="Top")){
    mtext(n)
  }
}

## Plot Fig. S1
set.seed(1)
pdf("Evolutionary trajectories steepness=c(1,2,1,1), alpha=0.5.pdf",width=10,height=7)
par(mfrow=c(4,length(ns)),mar=c(1,2,0,0),oma=c(3,5,3,1))
txt <- logical(3)
sides <- c("Bottom","Left","Top")
txt[3] <- TRUE
for(gg in tail(gs,4)){
  for(nn in ns){
    txt[2] <- ifelse(nn==1,TRUE,FALSE)
    txt[1] <- ifelse(gg==1,TRUE,FALSE)
    timeplots("steepness=c(0,0,0.8,1), alpha=1, mu=5",n=nn,g=gg,axes="sides[which(txt)]")
  }
  txt[3] <- FALSE
}

mtext(bquote("Number of time steps prior to reproduction, "~italic("n")),side=3,cex=1.2,outer=T,line=1.2)
mtext(bquote("Environmental grain,  "~italic("g")),side=2,outer=T,line=2.8,cex=1.2)
mtext("1                              0.9                              0.8                              0.7",side=2,outer=T,line=1.8)
dev.off()
####

## Create results summary plots (Fig. 2, S2, S4, S5) ##
# scenario: Provide simulation scenario with file name in format "steepness=c(s1,s2,u1,u2), alpha=a"
# ns:       Vector of number of time steps prior to reproduction (n), will be plotted as different lines
# gs:       Vector of environmental grains (g), will be plotted on x-axis.
# ydefault: Whether (T, default) or not to use default y-axis limits, which are from 0 to the highest RN slope produced in this scenario +0.05.
# lower:    Custom lower y-axis limit, only used if ydefault=F
# upper:    Custom upper y-axis limit, only used if ydefault=F
# legend:   Whether or not to print legend. Default T.
summaryplot <- function(scenario,ns=c(1,2,5,10),gs=seq(0.5,1,by=0.1),ydefault=T,lower,upper,legend=T,title=""){
  data <- readRDS(paste0("Summary ",scenario,".R"))
  gg <- as.numeric(dimnames(data)[[3]])
  if(length(gs)>length(gg)){ print("Too many x-axis points. Provide different a range.")}
  first <- which(gg==gs[1]) # Given that there are fewer x-axis points (gs) than simulation scenarios (gg), which simulation scenario corresponds to first x-axis point?
  
  # Set diverging red-blue color scheme
  require(RColorBrewer)
  pal <- brewer.pal(length(ns)+2,"RdBu") # Typically 6, then reduce down to 4.
  pal <- c(head(pal,length(pal)/2-1),tail(pal,length(pal)/2-1))
  pchs <- c(1,2,4,5,6,7,8) # Point types for plotting - length(pchs) should be equal to or longer than length(ns).
  
  plot(1:length(gs),data[1,1,first:length(gg)],ylim=c(ifelse(ydefault,0,lower),ifelse(ydefault,max(data[1,1,],na.rm=T)+0.05,upper)),type="n",xaxt="n",xlab=expression(paste("Environmental grain, ",italic("g"))),
       ylab=expression(paste("Reaction norm slope,  ",italic('\U03B2'))),cex.lab=1.1,main=title,cex.main=1.6)
  if(legend) legend("bottomright",legend=ns,col=pal,pch,lty=1,title=expression(paste("Time steps prior to reproduction, ",italic(n),":")),bty="n",pch=pchs)
  axis(side=1,at=1:length(gs),label=tail(seq(0,1,by=0.1),length(gs)))
  for(i in 1:length(ns)){
    points(1:length(gs),data[i,1,first:length(gg)],pch=pchs[i],cex=1.2,col=pal[i])
    lines(1:length(gs),data[i,1,first:length(gg)],col=pal[i])
    for(j in 1:length(gs)){
      arrows(j,data[i,1,first+j-1],j,data[i,1,first+j-1]+data[i,2,first+j-1],angle=90,length=0.1,col=pal[i])
      arrows(j,data[i,1,first+j-1],j,data[i,1,first+j-1]-data[i,2,first+j-1],angle=90,col=pal[i],length=0.1)
    }
  }
}

# Produce results figures

summaryplot("steepness=c(0,0,1,2), alpha=1") # 2B
summaryplot("steepness=c(0,0,1,2), alpha=0.5") # S2B
summaryplot("steepness=c(1,2,0,0), alpha=1") # 2A
summaryplot("steepness=c(1,2,0,0), alpha=0.5") # S2A

summaryplot("steepness=c(0.8,1,0,0), alpha=1",ydefault=F,lower=0,upper=0.4,legend=F) # S4A
summaryplot("steepness=c(0.8,1,0,0), alpha=0.5",ydefault=F,lower=0,upper=0.4,legend=F)
summaryplot("steepness=c(0,0,0.8,1), alpha=1",ydefault=F,lower=0.4,upper=0.7) # S4B
summaryplot("steepness=c(0,0,0.8,1), alpha=0.5",ydefault=F,lower=0.4,upper=0.7)


summaryplot("steepness=c(0.4,1,0.4,1), alpha=1",ydefault=F,lower=0.7,upper=0.9) # S5A
summaryplot("steepness=c(0.4,1,0.4,1), alpha=0.5",ydefault=F,lower=0.7,upper=0.9)
summaryplot("steepness=c(0.4,1,0.4,2), alpha=1",ydefault=F,lower=0.6,upper=0.9) # S5B
summaryplot("steepness=c(0.4,1,0.4,2), alpha=0.5",ydefault=F,lower=0.6,upper=0.9)
summaryplot("steepness=c(0.4,2,0.4,1), alpha=1",ydefault=F,lower=0.6,upper=0.9) # S5C
summaryplot("steepness=c(0.4,2,0.4,1), alpha=0.5",ydefault=F,lower=0.6,upper=0.9)
summaryplot("steepness=c(0.4,2,0.4,2), alpha=1",ydefault=F,lower=0.6,upper=0.9) # S5D
summaryplot("steepness=c(0.4,2,0.4,2), alpha=0.5",ydefault=F,lower=0.6,upper=0.9)


#### Calculate fitness means and variance (for predictions) ####

# Gaussian W(m)
# b is the  plasticity slope
# steepness[1:2] determines linear and exponent of maintenance costs
# steepness[3:4] determines linear and exponent of production costs
# Outputs: [1] Variance. [2] Arit.mean. [3] Geom.mean approx. [4] Geom.mean exact.
sim_gaus <- function(b,steepness=c(1,2,0,0)){
  es <- runif(2000000,-1,1)
  zs <- b*es
  ws <- dnorm(tail(abs(zs-es),length(zs)-1),0,0.4)
  maintenancecost <- cost(b,steepness=steepness[1:2])
  productioncost <- cost(abs(head(zs,length(zs)-1)-tail(zs,length(zs)-1)),steepness=steepness[3:4])
  ws <- ws-maintenancecost-productioncost
  logws <- log(ws)
  logws[which(is.nan(logws)==TRUE)] <- -10000
  return(c(var(ws),mean(ws),mean(ws)-(2*var(ws))/mean(ws),exp(mean(logws))))
}
#####

### Find exact maxima of predictions curves ####
maxima <- function(steepness=c(1,2,0,0),precision=3,type="arit"){
  bs <- seq(0,1,by=0.1)
  out <- matrix(NA,3,length(bs))
  for(rep in 1:3){
    for(i in 1:length(bs)){
      tmp <- sim_gaus(bs[i],steepness=steepness)
      out[rep,i] <- switch(type,var=tmp[1],arit=tmp[2],geomapprox=tmp[3],geomexact=tmp[4])
    }
  }
  p <- 1
  outmax <- apply(out,1,which.max)
  peak <- bs[outmax]
  #print(peak)
  peak <- round(mean(peak),digits=p)
  
  while(p < precision){
    p <- p+1
    bs <- seq(bs[min(outmax)-1],bs[max(outmax)+1],by=10^-p)
    out <- matrix(NA,3,length(bs))
    for(rep in 1:3){
      for(i in 1:length(bs)){
        tmp <- sim_gaus(bs[i],steepness=steepness)
        out[rep,i] <- switch(type,var=tmp[1],arit=tmp[2],geomapprox=tmp[3],geomexact=tmp[4])
      }
    }
    outmax <- apply(out,1,which.max)
    peak <- bs[outmax]
    #print(peak)
    peak <- round(mean(peak),digits=p)
  }
  return(peak)
}

### Plot predictions fig. 1 ####
par(mfrow=c(1,2),mar=c(4,4.1,2,0.3))
bs <- seq(0,1,by=0.01) #For any RN slope
vars <- ameans <- gmeansapprox <- gmeansexact <- numeric(length(bs))
for(i in 1:length(bs)){
  tmp <- sim_gaus(bs[i],steepness=c(1,2,0,0))
  vars[i] <- tmp[1]
  ameans[i] <- tmp[2]
  gmeansapprox[i] <- tmp[3] # Use this one if investigating scenarios that may have negative fitness values.
  gmeansexact[i] <- tmp[4]
  print(i)
}
arit1a <- maxima(steepness=c(1,2,0,0),precision=3,type="arit") #0.313
geomexact1a <- maxima(steepness=c(1,2,0,0),precision=3,type="geomexact") # 0.431
geomapprox1a <- maxima(steepness=c(1,2,0,0),precision=3,type="geomapprox") # 0.554

plot(bs,gmeansexact,type="l",col="Red",ylim=c(0,1),ylab="Fitness",xlab=expression(paste("Reaction norm slope, ",italic('\U03B2'))),
     cex.lab=1.3,lwd=2,main="A) Maintenance costs",cex.main=1.6)
lines(bs,ameans,lwd=2,col=rgb(0,0,1,0.6))
lines(bs,sqrt(vars))
arrows(arit1a,-0.05,arit1a,max(ameans),length=0,col=rgb(0,0,1,0.6),lty=2)
arrows(geomexact1a,-0.05,geomexact1a,max(gmeansexact),length=0,col=" Red",lty=2)
lines(bs,bs^2,lty=2)
legend("topleft",legend=c("Arithmetic mean","Geometric mean","Standard deviation","Cost of plasticity"),
       lwd=c(2,2,1,1),lty=c(1,1,1,2),col=c(rgb(0,0,1,0.6),"Red","Black","Black"),bty="n")
#mtext(text="A)",side=2,line=2.5,at=1,las=2,cex=1.4)

vars <- ameans <- gmeansapprox <- gmeansexact <- numeric(length(bs))
for(i in 1:length(bs)){
  tmp <- sim_gaus(bs[i],steepness=c(0,0,1,2))
  vars[i] <- tmp[1]
  ameans[i] <- tmp[2]
  gmeansapprox[i] <- tmp[3] # Use this one if investigating scenarios that may have negative fitness values.
  gmeansexact[i] <- tmp[4]
  print(i)
}

arit1b <- maxima(steepness=c(0,0,1,2),precision=3,type="arit") # 0.494
geomexact1b <- maxima(steepness=c(0,0,1,2),precision=3,type="geomexact") #0.17 (fail)
geomapprox1b <- maxima(steepness=c(0,0,1,2),precision=3,type="geomapprox") # 0.411

plot(bs,gmeansapprox,type="l",col="Red",ylim=c(0,1),ylab="Fitness",xlab=expression(paste("Reaction norm slope, ",italic('\U03B2'))),
     cex.lab=1.3,lwd=2,main="B) Production costs",cex.main=1.6)
#mtext(text="B)",side=2,line=2.5,at=1,las=2,cex=1.4)
lines(bs,ameans,lwd=2,col=rgb(0,0,1,0.6))
lines(bs,sqrt(vars))
arrows(arit1b,-0.05,arit1b,max(ameans),length=0,col=rgb(0,0,1,0.6),lty=2)
arrows(geomapprox1b,-0.05,geomapprox1b,max(gmeansapprox),length=0,col=" Red",lty=2)
dev.off()