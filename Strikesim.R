# Strikesim model: determine detection probability of strikes using playback experiments
rm(list = ls())
gc()
# USER SPECIFIED PARAMETERS -------------------------------------------------
RootDir =  "above1"  # Examples "current" or "above1" or "C:/Users/XXXX/BayesianMethods"
AnalysisFolder = 'Strikemod'  # Folder path within RootDir where analysis code is stored
# RunFile = 'BayesCalls2_10'       # Version of BayesCalls to run
DataFolder = 'simdata'  # Folder path within RootDir where raw data files stored
ResultsFolder = 'simresults'  # Folder path within RootDir where results files stored
Datafile =  'StrikeZero_BRS_2017R7_2017R10_2018R01_18Apr18.csv'  #Name of data file for analysis 
Envdata = ''  #Name of file with environmental data
ProjectName =  'Strikesim'  # Name of project (for user reference)
ProjectLocation =  'Kauai'  # Island name for data analysis
ProjectYear =  2018  # Year of analysis project 
# Nsims = 100
Jagsfile = 'Jags_Strikesim.jags'
Nchains = 20
Nburnin =  1500  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  100  # Number of adapting reps, default 100
Totalreps = 5000 # Total desired reps (ie # simulations making up posterior)
TrueNhits = 10
# END USER SPECIFIED PARAMETERS ---------------------------------------------
#
# Process user input and import data ----------------------------------------
existing_Packages<-as.list(installed.packages()[,1])
# Add new packages you might need to this line, to check if they're installed and install missing packages
required_Packages<-c("rstudioapi","readxl","lubridate","stringr","gtools","coda",
                     "lattice","rjags","jagsUI","parallel","doParallel","fitdistrplus","ggplot2")
missing_Packages<- required_Packages[!required_Packages%in% existing_Packages]
if(length(missing_Packages)>0)install.packages(pkgs =  missing_Packages)
invisible(lapply(required_Packages, require, character.only=T,quietly = T))
# Install libraries
library(rstudioapi)
library(readxl)
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
if (DataFolder=='') {
  loadfile1 = paste0(RootDir,"/",Datafile)
} else {
  loadfile1 = paste0(RootDir,"/",DataFolder,"/",Datafile)
}
#
if (Envdata!=''){
  if (DataFolder=='') {
    loadfile2 =  paste0(RootDir,"/",Envdata)
  } else {
    loadfile2 =  paste0(RootDir,"/",DataFolder,"/",Envdata)
  }
}
SaveResults = paste0(ProjectName,'_',ProjectLocation,'_',ProjectYear,'.rdata')
if (ResultsFolder=='') {
  SaveResults = paste0(RootDir,"/",SaveResults)
} else {
  SaveResults = paste0(RootDir,"/",ResultsFolder,"/",SaveResults)
}
#
DataFolder = paste0(RootDir,"/",DataFolder)
ResultsFolder = paste0(RootDir,"/",ResultsFolder)
AnalysisFolder = paste0(RootDir,"/",AnalysisFolder)
setwd(AnalysisFolder)
# Load data
data = read.csv(file = loadfile1, header = TRUE, sep = ",")
if (Envdata!=''){
  Edat = read.csv(file = loadfile2, header = TRUE, sep = ",")
}
# Process variables  -----------------------------------------------------
attach(data)
#
Nobs = dim(data)[1]
Nhits = hits
Wts = integer(length=Nobs) + TrueNhits
NightN = NightOfYear
Year = year
Year1 = min(Year)
YearN = Year-Year1+1
Nyrs = max(YearN)
DayN = integer(length=Nobs)
for (y in 1:Nyrs){
  Yr = Year1+y-1
  ii = which(Year==Yr)
  if (y==1){
    NightStrt = min(NightN[ii])
    
    DayN[ii] = NightN[ii] - NightStrt + 1
  }else{
    DayN[ii] =(y-1)*365 + NightN[ii] - NightStrt + 1
  }
}
WeekN = ceiling(DayN/7)
Mins = minutes
Site = as.factor(UMPLocationID)
SiteN = as.numeric(Site)
Sensor = as.factor(Sensor_Name)
SensorN = as.numeric(Sensor)
Area = Area_Short
Pole = Pole_Num
Period = factor(BRSPeriod, levels=c('SS15', 'SS135', 'SS195', 'Midnight',
                                  'SR195', 'SR135', 'SR15'))
PeriodN = as.numeric(Period)
tmp = unique(data.frame(Site, SiteN))
Sitelist = tmp[order(tmp$SiteN,tmp$Site),1:2];
NSite = dim(Sitelist)[1]
tmp = unique(data.frame(Period, PeriodN))
Periodlist = tmp[order(tmp$PeriodN,tmp$Period),1:2];
NPeriod = dim(Periodlist)[1]
tmp = unique(data.frame(Sensor, SensorN))
Sensorlist = tmp[order(tmp$SensorN,tmp$Sensor),1:2];
# Moon = 2*(Illu*Moonup)-1
Moon = Illu*Moonup
FS = flux_sensitive   # Flux sensitive variable as-is
LA = level_absolute+5.5 # Level absolute variable made positive by adding 5 
CL = log(1+click)     # Click variable re-scaled, log of click+1 
BU = 10*log(1+burst)     # Burst variable re-scaled, log of burst+1 times 10

detach(data)

# Fit preliminary glm to hits using signal metrics ------------------------
i = seq(1,Nobs)
df = data.frame(hits=pmin(1,Nhits/TrueNhits),FS = FS[i], LA = LA[i], CL = CL[i],BU = BU[i],
                FS2 = FS[i]^2, FS3 = FS[i]^3, LA2 = LA[i]^2, LA3 = LA[i]^3, FSLA = FS[i]*LA[i])
model <- glm(hits ~ FS + FS2 + FS3 + LA + BU + CL + FSLA,
             family=binomial(link='logit'),weights=Wts, data=df)
summary(model)
df$prob = predict(model, newdata = df, type = "response")
plot(df$FS[df$hits<1],df$LA[df$hits<1],col="blue",ylim = c(0,4),xlim = c(0,100))
points(df$FS[df$hits==1],df$LA[df$hits==1],col="red",ylim = c(0,4),xlim = c(0,100))
newdat = df[1:800,]
Prob1 = matrix(nrow=20,ncol=20)
Prob2 = matrix(nrow=20,ncol=20)
cntr = 0
cntrC = 0
for(j in c(0.01,.3)){
  cntrx = 0
  cntrC = cntrC +1
  for(f in seq(10,100,length.out = 20)){
    cntrx = cntrx+1
    cntry = 0    
    for(l in seq(.5,4,length.out = 20)){
      cntry = cntry+1
      cntr = cntr+1
      newdat$hits[cntr]=0
      newdat$FS[cntr]=f
      newdat$FS2[cntr]=f^2
      newdat$FS3[cntr]=f^3
      newdat$LA[cntr]=l
      newdat$LA2[cntr]=l^2
      newdat$BU[cntr]=j
      newdat$CL[cntr]=.001
      newdat$FSLA[cntr]=f*l
      if (cntrC == 1){
        Prob1[cntrx,cntry] = predict(model, newdata = newdat[cntr,], type = "response")
      }else{
        Prob2[cntrx,cntry] = predict(model, newdata = newdat[cntr,], type = "response")
      }
    }
  }
}
thrdvar = 'Burst'
newdat$prob = predict(model, newdata = newdat, type = "response")
xx = seq(10,100,length.out = 20); yy = seq(.5,4,length.out = 20)
library(lattice)
levelplot(Prob1, data = NULL, aspect = "fill",
          xlim = c(0,100),ylim = c(0.5,4),
          row.values = xx, column.values = yy,
          xlab = 'Flux Sensitive', ylab = 'Level Absolute',
          main= paste0(thrdvar, '= low'))
levelplot(Prob2, data = NULL, aspect = "fill",
          xlim = c(0,100),ylim = c(0.5,4),
          row.values = xx, column.values = yy,
          xlab = 'Flux Sensitive', ylab = 'Level Absolute',
          main= paste0(thrdvar, '= high'))
# Extract coefficients and std erros as priors for Bayesian model
Bpar1 =  summary(model)$coefficients[, 1] #
Bpar2 =  summary(model)$coefficients[, 2]
# Convert standard erros of param estimates to precision values
Bpar2 = 1/Bpar2^2

# Set up JAGS -------------------------------------------------------------

set.seed(123)
# For 8 cores, set next lines to 8, 2750, 1500, 100 to get 10,000 reps
# Nsim = 2120      # 2750
# Nburnin = 1500  # 2000
# Nadapt = 100

# For parallel (comment out for serial)
cores<-detectCores()
cores = min(cores, Nchains)
cl <- makeCluster(cores[1])
registerDoParallel(cl)

data <- list(Hits=Nhits,True=Wts,NObs=NObs,Site=SiteN,NSite=NSite,
             Yr=YearN,Tpr=PeriodN,Nyrs=Nyrs,Nprs=NPeriod,moon=Moon,
             FS=FS,LA=LA,CL=CL,BU=BU,Bpar1=Bpar1,Bpar2=Bpar2) 

# Inits: Best to generate initial values using function
inits <- function(){
  list(sigT=runif(1,0.3,0.6),sigN=runif(1,1,5),sigS=runif(1,2,10),sigW=runif(1,.1,.5),
       theta=runif(1,-.1,.1),Dispers=runif(1,.15,.25))
}
# List of parameters to monitor:
params <- c('theta','sigT','sigS','sigN','sigW','Dispers','peakTemp',
            'C0','C','Cs','Csite','B','Temp','eps') # 


# Run JAGS ----------------------------------------------------------------

# For detailed output stats & DIC use jags(), but do not inlcude Temporal or E
#  in params since they take up too much memory and it crashes
# To get Temporal matrix, use jags.basic() to save just basic mcmc list
#  which uses much less memory but does not save associated stats & DIC
out <- jags.basic(data = data,
                  inits = inits,
                  parameters.to.save = params,
                  model.file = Jagsfile,
                  n.chains = Nchains,
                  n.adapt = Nadapt,
                  n.iter = Nsim,
                  n.burnin = Nburnin,
                  n.thin = 1,
                  parallel=TRUE,
                  n.cores=cores)
# Diagnostics -------------------------------------------------------------
stopCluster(cl)
vn = varnames(out)
outmat = as.matrix(out)
reps = dim(outmat); reps = reps[1]
outdf = as.data.frame(outmat); rm(outmat)
s = summary(out) 
s_stats = data.frame(s$statistics)
s_quantiles = data.frame(s$quantiles) 
# 
# rhat = gelman.diag(out); rhat = data.frame(rhat$psrf);
# Deviance = s_stats["deviance",1]
# Vardev = s_stats["deviance",2]^2; pd = Vardev/2;
# DIC = Deviance + 2*pd
save(list = ls(all.names = TRUE),file=SaveResults)


