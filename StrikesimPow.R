# Strike simulation experiment Power analysis: Power to Detect an Environmental Co-variate
rm(list = ls())
gc()
# USER SPECIFIED PARAMETERS -------------------------------------------------
RootDir =  "above1"  # Examples "current" or "above1" or "C:/Users/XXXX/BayesianMethods"
AnalysisFolder = 'Strikemod'  # Folder path within RootDir where analysis code is stored
# RunFile = 'BayesCalls2_10'       # Version of BayesCalls to run
DataFolder = 'simdata'  # Folder path within RootDir where raw data files stored
ResultsFolder = 'simresults'  # Folder path within RootDir where results files stored
Resultsfile = "Strikesim_Kauai_2018.rdata"
ProjectName =  'Strikesim'  # Name of project (for user reference)
ProjectLocation =  'Kauai'  # Island name for data analysis
ProjectYear =  2018  # Year of analysis project 
# Nsims = 100
JagsfileP = 'Jags_Strikesimpow.jags'
# Set parameters for Power Analysis
simreps = 10       # Number reps for Power analysis (recomend at least 10)
EVtype = 1       # Type of Environmental variable (EV), 1 = categorical, 2 = continuous
EVeff = -0.1    # True effect of EV on Strike Detection (-.1 = 10% lower prob per unit incr. in EV)
EVfrq = .5       # If EV is categorical (1-0 var, e.g. rain), expected frequencey that EV = 1 
EVsd = 1         # If EV is continuous (centered so mean=0), std dev of values (default=1, use "scale()")
NSiteA = 15      # Total Number of Sites Monitored 
RecPsite = c(20,200) # Number of experimental trials per site, c(min, max)
P_signif = 0.95     # Desired level of certainty for CI and P values

# END USER SPECIFIED PARAMETERS ---------------------------------------------
#
# Process user input and import data ----------------------------------------
existing_Packages<-as.list(installed.packages()[,1])
# Add new packages you might need to this line, to check if they're installed and install missing packages
required_Packages<-c("rstudioapi","stats","lubridate","stringr","gtools","coda","mcmcplots",
                     "lattice","rjags","jagsUI","parallel","doParallel","fitdistrplus","ggplot2")
missing_Packages<- required_Packages[!required_Packages%in% existing_Packages]
if(length(missing_Packages)>0)install.packages(pkgs =  missing_Packages)
invisible(lapply(required_Packages, require, character.only=T,quietly = T))
# Install libraries
library(rstudioapi)
library(stats)
library(lubridate)
library(stringr)
library(gtools)
library(coda)
library(lattice)
library(rjags)
library(jagsUI)
library(parallel)
library(doParallel)
library(fitdistrplus)
library(mcmcplots)
library(ggplot2)
rm('existing_Packages')
rm('missing_Packages')
rm('required_Packages')
# Process paths and directory structure 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (RootDir=='current'){
  RootDir = getwd()
} else if(RootDir=='above1'){
  tmp1 = getwd()
  setwd('..')
  RootDir = getwd()
  setwd(tmp1)
  rm('tmp1')
} else if(RootDir=='above2'){
  tmp1 = getwd()
  setwd('..')
  setwd('..')
  RootDir = getwd()
  setwd(tmp1)
  rm('tmp1')
} else {
  RootDir = RootDir
}
#
if (ResultsFolder=='') {
  loadfile1 = paste0(RootDir,"/",Resultsfile)
} else {
  loadfile1 = paste0(RootDir,"/",ResultsFolder,"/",Resultsfile)
}
#
ResultsFolder = paste0(RootDir,"/",ResultsFolder)
AnalysisFolder = paste0(RootDir,"/",AnalysisFolder)
setwd(AnalysisFolder)
# Load results
load(loadfile1)
B = s_stats$Mean[which(startsWith(vn,"B["))]
Bpar1 = B; 
Bpar2 = 1/(s_stats$SD[which(startsWith(vn,"B["))])^2
theta = s_stats$Mean[which(startsWith(vn,"theta"))]
Teff = s_stats$Mean[which(startsWith(vn,"Teff"))]
sigS = s_stats$Mean[which(startsWith(vn,"sigS"))]
sigT = s_stats$Mean[which(startsWith(vn,"sigT"))]
gam = s_stats$Mean[which(startsWith(vn,"gam"))]
sampsize = round(seq(RecPsite[1],RecPsite[2],length.out = simreps))
sampsizeT = sampsize*NSiteA
alph = 1-P_signif
Ncores<-detectCores()
Nchains = min(20,Ncores)
Nburnin =  500  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  100  # Number of adapting reps, default 100
Totalreps = 1000 # Total desired reps (ie # simulations making up posterior)
Nsim =  Totalreps/Nchains + Nburnin 
cores = min(Ncores, Nchains)
cl <- makeCluster(cores[1])
registerDoParallel(cl)
#
Signif = numeric()
for (r in 1:simreps){
  NObs = NSiteA*sampsize[r] 
  ii = sample(seq(1,Nobs),NObs,replace = TRUE)
  FSr = FS[ii]; LAr = LA[ii]; BUr = BU[ii]; CLr = CL[ii]
  Per = PeriodN[ii]  
  MoonR = Moon[ii]
  SiteR = numeric()
  Seffr= numeric()
  for (s in 1:NSiteA){
    Seffr[s] = rnorm(1,0,sigS)
    iii = seq((s-1)*sampsize[r]+1,(s-1)*sampsize[r]+sampsize[r])
    SiteR[iii] = s
  }
  # EV = Environmental Variable of interest (for which we want to examine power)
  if (EVtype==1){
    EV = rbinom(NObs,1,EVfrq)
  }else{
    EV = rnorm(NObs,0, EVsd)
  }  
  BaseProb_lgt = B[1]+B[2]*FSr+B[3]*FSr^2+B[4]*FSr^3+B[5]*LAr+B[6]*BUr+B[7]*CLr+B[8]*FSr*LAr
  Mn_NoEV = inv.logit(BaseProb_lgt + theta*MoonR + Teff[Per] + Seffr[SiteR]) 
  base = median(Mn_NoEV)
  new = base + EVeff
  # Compute true value of phi, the logit param resulting in specified % change in EV
  phiT = logit(new) - logit(base)
  # Generate simulated strike data using beta-binomial distribution
  Mn_EV = inv.logit(BaseProb_lgt + theta*MoonR + Teff[Per] + Seffr[SiteR] + phiT*EV) 
  a = Mn_EV/gam 
  b = (1-Mn_EV)/gam
  Prob = rbeta(NObs,a,b)
  # plot(Mn_EV,Prob)
  TruHts = rep(TrueNhits,NObs)
  Hits = rbinom(NObs,TruHts,Prob)
  # Now estimate phi using Bayesian model
  data <- list(Hits=Hits,True=TruHts,NObs=NObs,Site=SiteR,NSite=NSiteA,
               Tpr=Per,Nprs=NPeriod,moon=MoonR,EV=EV,
               FS=FSr,LA=LAr,CL=CLr,BU=BUr,Bpar1=Bpar1,Bpar2=Bpar2) # Yr=YearN,Nyrs=Nyrs,
  #
  # Inits: Best to generate initial values using function
  inits <- function(){
    list(sigT=runif(1,.3*sigT,1.8*sigT),sigS=runif(1,.3*sigS,1.8*sigS))
         #gam=runif(1,.9*gam,1.1*gam),theta=runif(1,.9*theta,1.1*theta),
         #phi=runif(1,-.5,.5))
  }
  # List of parameters to monitor:
  paramsP <- c('theta','gam','sigT','sigS','phi') # 
  outP <- jags.basic(data = data,
                    inits = inits,
                    parameters.to.save = paramsP,
                    model.file = JagsfileP,
                    modules=c('mix'),
                    n.chains = Nchains,
                    n.adapt = Nadapt,
                    n.iter = Nsim,
                    n.burnin = Nburnin,
                    n.thin = 1,
                    parallel=TRUE,
                    n.cores=cores)
  sP = summary(outP) 
  sP_stats = data.frame(sP$statistics)
  sP_quantiles = data.frame(sP$quantiles)
  vnP = varnames(outP)
  pp = which(vnP=="phi")
  outmatP = as.matrix(outP)
  outdfP = as.data.frame(outmatP); rm(outmatP)
  phiest = sP_stats$Mean[pp]
  phipost = outdfP$phi
  if (phiT<0 & phiest<0 & quantile(phipost,P_signif+alph/2)<0){
    Sigcriteria = 1
  }else if(phiT>0 & phiest>0 & quantile(phipost,alph/2)>0){
    Sigcriteria = 1
  }else{
    Sigcriteria = 0
  }
  if(Sigcriteria == 1 & phiest>sP_quantiles[pp,1] & phiest<sP_quantiles[pp,5]){
    Signif[r] = 1
  }else{
    Signif[r] = 0
  }
}
pfp = c(which(paramsP=='gam'),which(startsWith(paramsP,'sig')),which(paramsP=='phi')) #,which(params=='Yeff')
#
for (i in pfp){
  parnm = paramsP[i]
  traplot(outP,parnm)
  denplot(outP,parnm,ci=.9,collapse = TRUE)
}
if (sum(Signif)==0){
  print("The specified effect size cannot be reliably detected with the current range of sample sizes")
  print(" (there is insufficient Power). Either increase sample size or specify larger effect.")
}else{
  dfP = data.frame(Ssize = sampsize, SsizeT = sampsizeT, Signif = Signif)
  fitP = glm(Signif ~ Ssize, family = binomial, data = dfP)
  summary(fitP)
  dfP$Power = predict(fitP, newdata = dfP, type = "response")
  subtitl = paste0("N sites = ", NSiteA,", True effect size = ", 100*EVeff,
                   "% change in detection probability per unit increase in variable")
  pltP = ggplot( dfP, aes(x=Ssize, y=Signif)) +
    geom_point() +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    ggtitle("Power vs Sample Size: Detecting a Variables Effect on Strike Detection" ,
            subtitle = subtitl ) +
    xlab("Sample size (playback experiments per site)") +
    ylab("Power (probability of correctly estimating effect)")
  print(pltP)
}
  

