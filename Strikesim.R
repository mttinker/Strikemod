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
required_Packages<-c("rstudioapi","readxl","lubridate","stringr","gtools","coda","mcmcplots",
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
dat = read.csv(file = loadfile1, header = TRUE, sep = ",")
if (Envdata!=''){
  Edat = read.csv(file = loadfile2, header = TRUE, sep = ",")
}
# Process variables  -----------------------------------------------------
Nsim =  Totalreps/Nchains + Nburnin  # Total # MCMS sims: Actual saved reps = (Nsim-Nburnin) * (num Cores)
attach(dat)
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
LA = level_absolute+5.5 # Level absolute variable made positive by adding 5.5 
CL = log(1+click)     # Click variable re-scaled, log of click+1 
BU = 10*log(1+burst)     # Burst variable re-scaled, log of burst+1 times 10

detach(dat)

# Fit preliminary glm to hits using signal metrics ------------------------
i = seq(1,Nobs)
df = data.frame(hits=pmin(1,Nhits/TrueNhits),FS = FS[i], LA = LA[i], CL = CL[i],BU = BU[i],
                FS2 = FS[i]^2, FS3 = FS[i]^3, LA2 = LA[i]^2, LA3 = LA[i]^3, FSLA = FS[i]*LA[i])
model <- glm(hits ~ FS + FS2 + FS3 + LA + BU + CL + FSLA,
             family=binomial(link='logit'),weights=Wts, data=df)
plotGLMmod = 0
if(plotGLMmod==1){
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
  newdat$prob = predict(model, newdata = newdat, type = "response")
  thrdvar = 'Burst'
  xx = seq(10,100,length.out = 20); yy = seq(.5,4,length.out = 20)
  # library(lattice)
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
}
# Extract coefficients and std erros as priors for Bayesian model
Bpar1 =  summary(model)$coefficients[, 1] # Mean param estimats
Bpar2 =  summary(model)$coefficients[, 2] # Std Err of param estimates
# Convert standard errors of param estimates to precision values
Bpar2 = 1/(2*Bpar2)^2
#
# Set up JAGS -------------------------------------------------------------
#
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
#
data <- list(Hits=Nhits,True=Wts,NObs=Nobs,Site=SiteN,NSite=NSite,
             Tpr=PeriodN,Nprs=NPeriod,moon=Moon,
             FS=FS,LA=LA,CL=CL,BU=BU,Bpar1=Bpar1,Bpar2=Bpar2) # Yr=YearN,Nyrs=Nyrs,
#
# Inits: Best to generate initial values using function
inits <- function(){
  list(sigT=runif(1,.1,.5),sigS=runif(1,.1,.5),
       theta=runif(1,-.1,.1)) # sigY=runif(1,0.1,0.5),
}
# List of parameters to monitor:
params <- c('theta','sigT','sigS',
            'Yeff','Teff','Seff','B') # 
#
# Run JAGS ----------------------------------------------------------------
#
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
#
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
#
pfp = c(which(params=='theta'),which(startsWith(params,'sig'))) #,which(params=='Yeff')
#
for (i in pfp){
  parnm = params[i]
  traplot(out,parnm)
  denplot(out,parnm,ci=.9,collapse = TRUE)
}
#
# Effect plots ------------------------------------------------------------
#
parnum = which(params=='Teff')
Levlabels = as.character(Periodlist$Period)
caterplot(out,params[parnum],denstrip = FALSE, reorder=FALSE,
          quantiles=list(outer=c(0.025,0.975),inner=c(0.1666,0.8333)),lwd=c(.1,4),
          labels=Levlabels, labels.loc = 'above',las = 0, cex.labels = .8)
title(main = "Effect of time period", font.main = 4,
      xlab = "Logit effect")

parnum = which(params=='Seff')
Levlabels = as.character(Sitelist$Site)
caterplot(out,params[parnum],denstrip = FALSE, reorder=FALSE,
          quantiles=list(outer=c(0.025,0.975),inner=c(0.1666,0.8333)),lwd=c(.1,4),
          labels=Levlabels, labels.loc = 'above',las = 0, cex.labels = .8)
title(main = "Effect of Site", font.main = 4,
      xlab = "Logit effect")

B = s_stats$Mean[which(startsWith(vn,"B["))]
intcpt = numeric(length = 800)
FS1 = intcpt; FS2 = intcpt; FS3 = intcpt; CL1 = intcpt; 
LA1 = intcpt; BU = intcpt; FSLA = intcpt; Prob = intcpt; 
intcpt = 1+intcpt 
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
      FS1[cntr]=f
      FS2[cntr]=f^2
      FS3[cntr]=f^3
      LA1[cntr]=l
      BU[cntr]=j
      CL1[cntr]=.001
      FSLA[cntr]=f*l
      Prob[cntr] = min(1,inv.logit(B[1]+B[2]*FS1[cntr]+B[3]*FS2[cntr]+B[4]*FS3[cntr]+
                                         B[5]*LA1[cntr]+B[6]*BU[cntr]+B[7]*CL1[cntr]+B[8]*FSLA[cntr]))
      if (cntrC == 1){
        Prob1[cntrx,cntry] = Prob[cntr]
      }else{
        Prob2[cntrx,cntry] = Prob[cntr]
      }
    }
  }
}
xx = seq(10,100,length.out = 20); yy = seq(.5,4,length.out = 20)
plt1 = levelplot(Prob1, data = NULL, aspect = "fill",
                xlim = c(10,100),ylim = c(0.5,4),
                row.values = xx, column.values = yy,
                xlab = 'Flux Sensitive', ylab = 'Level Absolute (+5.5)',
                main="Effect of Signal, Flux Sensitive X Level Absolute, Low Burst",
                sub="(log(1+Burst)*10=0.01)",font.sub = 2)
print(plt1)

plt2 = levelplot(Prob2, data = NULL, aspect = "fill",
                xlim = c(10,100),ylim = c(0.5,4),
                row.values = xx, column.values = yy,
                xlab = 'Flux Sensitive', ylab = 'Level Absolute (+5.5)',
                main="Effect of Signal, Flux Sensitive X Level Absolute, High Burst",
                sub="(log(1+Burst)*10=0.3)",font.sub = 2)
print(plt2)

FSi = mean(FS)
LAi = mean(LA)
BUi = mean(BU)
CLi = mean(CL)
B = s_stats$Mean[which(startsWith(vn,"B["))]
theta = s_stats$Mean[which(startsWith(vn,"theta"))]
iii = which(Periodlist$Period=="SS15" | Periodlist$Period=="Midnight" | Periodlist$Period=="SR15")
BaseProb_lgt = B[1]+B[2]*FSi+B[3]*FSi^2+B[4]*FSi^3+B[5]*LAi+B[6]*BUi+B[7]*CLi+B[8]*FSi*LAi
Seff = s_stats$Mean[which(startsWith(vn,"Teff"))]
SS15eff = Seff[iii[1]]
Midnteff = Seff[iii[2]]
SR15eff = Seff[iii[3]]
ProbDetect = numeric()
ProbDetect[1] = inv.logit(BaseProb_lgt + SS15eff)
ProbDetect[2] = inv.logit(BaseProb_lgt + Midnteff)
ProbDetect[3] = inv.logit(BaseProb_lgt + SR15eff)
ProbDetect[4] = inv.logit(BaseProb_lgt + SS15eff + theta)
ProbDetect[5] = inv.logit(BaseProb_lgt + Midnteff + theta)
ProbDetect[6] = inv.logit(BaseProb_lgt + SR15eff + theta)
PrbDtct = matrix(nrow = Totalreps, ncol = 6)
ii = which(startsWith(vn,"B["))
B = outdf[,ii]
theta = outdf[,startsWith(vn,"theta")]
BaseProb_lgt = B[,1]+B[,2]*FSi+B[,3]*FSi^2+B[,4]*FSi^3+B[,5]*LAi+B[,6]*BUi+B[,7]*CLi+B[,8]*FSi*LAi
Seff = outdf[,which(startsWith(vn,"Teff"))]
SS15eff = Seff[,iii[1]]
Midnteff = Seff[,iii[2]]
SR15eff = Seff[,iii[3]]
PrbDtct[,1] = inv.logit(BaseProb_lgt + SS15eff)
PrbDtct[,2] = inv.logit(BaseProb_lgt + Midnteff)
PrbDtct[,3] = inv.logit(BaseProb_lgt + SR15eff)
PrbDtct[,4] = inv.logit(BaseProb_lgt + SS15eff + theta)
PrbDtct[,5] = inv.logit(BaseProb_lgt + Midnteff + theta)
PrbDtct[,6] = inv.logit(BaseProb_lgt + SR15eff + theta)
PrbDtctCI = apply(PrbDtct,2,quantile,probs =c(.0275,.975),na.rm=TRUE) 
dfTper3 = data.frame(Period = rep((Periodlist$Period[iii]),2),
                     Moon = c(rep(c("No moon"),3),rep(c("Full moon"),3)),
                     Mean = ProbDetect, Lo = PrbDtctCI[1,], Hi = PrbDtctCI[2,])
pd <- position_dodge(0.1) # move errorbars for groups .05 to the left and right
plt3 = ggplot(dfTper3, aes(x=Period, y=Mean, colour=Moon, group=Moon)) + 
  geom_errorbar(aes(ymin=Lo, ymax=Hi), colour="black", width=.1, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=5) + # 21 is filled circle
  xlab("Time Period") +
  ylab("Probability of Detection") +
  scale_colour_hue(name="Moon Illumination",    # Legend label, use darker colors
                   breaks=c("No moon", "Full moon"),
                   labels=c("No moon", "Full moon"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("Effect of Time Period and Moon Illumination on Strike Detection Probability",
          subtitle = "For average site and mean values of signal quality metrics") +
  expand_limits(y=0) +                        # Expand y range
  # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))               # Position legend in bottom right
print(plt3)
