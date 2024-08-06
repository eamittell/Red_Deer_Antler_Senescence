#### Models used in manuscript: "Antler size in red deer: declining selection and increasing genetic variance with age, but little senescence"
### Red deer antler size traits and annual breeding success
## Updated 23rd July 2024
# Lizy

library(MCMCglmm)

# File path to data
file_path_data <- "/Users/umer/Documents/Work/Evol_Climate_Change/Analyses/RedDeer/AntlerPaper_DataAndScripts/"
file_path_data <- ""

# Read in and format the pedigree
pedigree <- read.csv(paste(file_path_data,"PedigreeRandomIds.csv",sep=""))
# Format
pedigree$animal <- as.factor(pedigree$animal)
pedigree$dam <- as.factor(pedigree$dam)
pedigree$sire <- as.factor(pedigree$sire)
str(pedigree)

#########################
###### Antler Form ######
#########################
# Read in data
DataForm <- read.csv(paste(file_path_data,"DataFormRandomIds.csv",sep=""))
DataForm$animal <- as.factor(DataForm$Code)
str(DataForm)

## Univariate animal model
prior_mod_form = list(G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*1000),G2 = list(V = diag(1), nu = 0.002),G3 = list(V = diag(1), nu = 0.002)), R = list(V = diag(1), nu = 1.002))
mod_form <- MCMCglmm(Form~Age + I(Age^2) + I(Age^3) + LastSeen, random = ~animal + Year + Code, rcov = ~units, family = "gaussian", prior = prior_mod_form, data = DataForm, pedigree = pedigree, nitt = 2600000, thin = 500, burnin = 50000)

## Trivariate animal model
prior_tri_form = list(G = list(G1 = list(V = diag(3), nu = 0.002, alpha.mu = rep(0,3), alpha.V = diag(3)*1000),G2 = list(V = diag(3), nu = 0.002),G3 = list(V = diag(3), nu = 0.002)), R = list(V = diag(3), nu = 1.002))
tri_mod_form <- MCMCglmm(cbind(y_form,p_form,e_form)~trait -1 + trait:Age + trait:LastSeen, random = ~us(trait):animal + us(trait):Year + us(trait):Code, rcov = ~idh(trait):units, family = c("gaussian","gaussian","gaussian"), prior = prior_tri_form, data = DataForm, pedigree = pedigree, nitt = 8600000, thin = 500, burnin = 50000)

#################
###### ABS ######
#################
# Read in data
DataForm <- read.csv(paste(file_path_data,"DataFormRandomIds.csv",sep=""))
DataForm$animal <- as.factor(DataForm$Code)
str(DataForm)

## All ages together univariate animal model
prior1_abs <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 1.002),G3 = list(V = 1, nu = 1.002)), R = list(V = 1, nu = 2.002))
model_abs <- MCMCglmm(ABS~ Age+I(Age**2)+I(Age**3)+LastSeen,random=~animal+StagYear+Code,rcov=~units,family="poisson",prior=prior1_abs,data=DataForm,pedigree=pedigree,nitt=1205000,thin=500,burnin=5000)

## All ages together univariate animal model -- age^3 dropped as insignificant; this doesn't change the results
prior1_abs <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 1.002),G3 = list(V = 1, nu = 1.002)), R = list(V = 1, nu = 2.002))
model_abs2 <- MCMCglmm(ABS~ Age+I(Age**2)+LastSeen,random=~animal+StagYear+Code,rcov=~units,family="poisson",prior=prior1_abs,data=DataForm,pedigree=pedigree,nitt=1205000,thin=500,burnin=5000)

## Separate univariate animal models for each age class
# Young
young2 <- DataForm[which(DataForm$AgeClass==1),]
prior2_abs <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 1.002),G3 = list(V = 1, nu = 1.002)), R = list(V = 1, nu = 2.002))
young_abs <- MCMCglmm(ABSY~ Age+LastSeen, random=~animal+StagYear+Code,rcov=~units,family="poisson",prior=prior2_abs,data=young2,pedigree=pedigree,nitt=1205000,thin=500,burnin=5000) 
# Prime
prime2 <- DataForm[which(DataForm$AgeClass==2),]
prime_abs <- MCMCglmm(ABSP~ Age+LastSeen, random=~animal+StagYear+Code,rcov=~units,family="poisson",prior=prior2_abs,data=prime2,pedigree=pedigree,nitt=5010000,thin=1000,burnin=10000)
# Old -- without animal term
old2 <- DataForm[which(DataForm$AgeClass==3),]
prior2_absO <- list(G = list(G1 = list(V = 1, nu = 1.002),G2 = list(V = 1, nu = 1.002)), R = list(V = 1, nu = 2.002))
old_abs2 <- MCMCglmm(ABSO~Age, random=~StagYear+Code,rcov=~units,family="poisson",prior=prior2_absO,data=old2,nitt=201000,thin=100,burnin=1000)

## Trivariate model with all age classes -- animal component only fitted for young and prime
prior3b_abs <- list(G = list(G1 = list(V = diag(2), nu = 1.002, alpha.mu = rep(0,2), alpha.V = diag(2)*1000),G2 = list(V = diag(3), nu = 1.002),G3 = list(V = diag(3), nu = 1.002)), R = list(V = diag(3), nu = 2.002))
tri_abs <- MCMCglmm(cbind(ABSY,ABSP,ABSO)~trait -1 + trait:Age + trait:LastSeen, random=~us(at.level(trait,c("ABSY","ABSP"))):animal+us(trait):StagYear+us(trait):Code,rcov=~idh(trait):units,family=c("poisson","poisson","poisson"),prior=prior3b_abs,data=DataForm,pedigree=pedigree,nitt=15005000,thin=1000,burnin=5000)

####################################
####### PHENOTYPIC SELECTION #######
####################################
# Read in data
DataForm <- read.csv(paste(file_path_data,"DataFormRandomIds.csv",sep=""))
DataForm$animal <- as.factor(DataForm$Code)
str(DataForm)
###### Format data
DataFormS <- DataForm[which(!is.na(DataForm$Form)),]
DataFormABS <- DataFormS[which(!is.na(DataFormS$ABS)),]

# Univariate form all ages together
prior1 <- list(G = list(G1 = list(V = 1, nu = 1.002),G2 = list(V = 1, nu = 1.002)), R = list(V = 1, nu = 2.002))
absF <- MCMCglmm(ABS~Age+I(Age**2) +LastSeen + Form,random=~Year+Code,rcov=~units,family="poisson",prior=prior1,data=DataFormABS,pedigree=pedigree,nitt=305000,thin=100,burnin=5000)

# Univariate with interaction between age class and antler form
absFCH11 <- MCMCglmm(ABS ~ Age + Form:AgeClass + LastSeen -1, random=~Year+Code,rcov=~units,family="poisson",prior=prior1,data=DataFormABS,pedigree=pedigree,nitt=305000,thin=100,burnin=5000)

#################################
####### GENETIC SELECTION #######
#################################
# All ages together
prior2 <- list(G = list(G1 = list(V = diag(2), nu = 1.002, alpha.mu = rep(0,2), alpha.V = diag(2)*1000),G2 = list(V = diag(2), nu = 1.002),G3 = list(V = diag(2), nu = 1.002)), R = list(V = diag(2), nu = 2.002))
bi_abs_form <- MCMCglmm(cbind(ABS,Form)~trait-1+trait:Age+trait:I(Age**2)+trait:I(Age**3)+trait:LastSeen, random=~us(trait):animal+us(trait):Year+us(trait):Code,rcov=~us(trait):units,family=c("poisson","gaussian"),prior=prior2,data=DataForm,pedigree=pedigree,,nitt=2010000,thin=200,burnin=10000)

# Young
bi_abs_formY <- MCMCglmm(cbind(ABSY,y_form)~trait -1 + trait:Age + trait:LastSeen, random=~us(trait):animal+us(trait):Year+us(trait):Code,rcov=~us(trait):units,family=c("poisson","gaussian"),prior=prior2,data=DataForm,pedigree=pedigree,nitt=2010000,thin=200,burnin=10000)
# Prime
prior2b <- list(G = list(G1 = list(V = diag(2), nu = 3.002, alpha.mu = rep(0,2), alpha.V = diag(2)*1000),G2 = list(V = diag(2), nu = 1.002),G3 = list(V = diag(2), nu = 1.002)), R = list(V = diag(2), nu = 2.002))
bi_abs_formP <- MCMCglmm(cbind(ABSP,p_form)~trait -1 + trait:Age+ trait:LastSeen, random=~us(trait):animal+us(trait):Year+us(trait):Code,rcov=~us(trait):units,family=c("poisson","gaussian"),prior=prior2b,data=DataForm,pedigree=pedigree,nitt=2010000,thin=200,burnin=10000)
# Old -- no Va as doesn't estimate it well for ABS throughout all models
prior3 <- list(G = list(G1 = list(V = diag(2), nu = 1.002),G2 = list(V = diag(2), nu = 1.002)), R = list(V = diag(2), nu = 2.002))
bi_abs_formO <- MCMCglmm(cbind(ABSO,e_form)~trait -1 + trait:Age+ trait:LastSeen, random=~us(trait):Year+us(trait):Code,rcov=~us(trait):units,family=c("poisson","gaussian"),prior=prior3,data=DataForm,pedigree=pedigree,nitt=50500000,thin=1000,burnin=50000)

#########################################
###### Models in the Supplementary ######
#########################################
# Read in data
DataForm <- read.csv(paste(file_path_data,"DataFormRandomIds.csv",sep=""))
DataForm$animal <- as.factor(DataForm$Code)
str(DataForm)

# All traits scaled below

#### Antler Mass ####
mass <- DataForm[which(!is.na(DataForm$Mass)),]
mass$MassS <- scale(mass$Mass)
prior1 <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 2.002),G3 = list(V = 1, nu = 2.002)), R = list(V = 1, nu = 2.002))
mass_modS <- MCMCglmm(MassS~ Age + I(Age^2), random=~animal+Year+Code,rcov=~units,family="gaussian",prior=prior1,data=mass,pedigree= pedigree,nitt=1005000,thin=500,burnin=5000)

young_mass <- mass[which(mass$Age<8),]
young_mass$MassS <- scale(young_mass$Mass)
prior1 <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 2.002),G3 = list(V = 1, nu = 2.002)), R = list(V = 1, nu = 2.002))
young_mass_modS <- MCMCglmm(MassS ~ Age, random=~animal+Year+Code,rcov=~units,family="gaussian",prior=prior1,data=young_mass,pedigree= pedigree,nitt=505000,thin=500,burnin=5000)

prime_mass <- mass[which(mass$Age<8),]
prime_mass <- prime_mass[which(prime_mass$Age<12),]
prime_mass$MassS <- scale(prime_mass$Mass)
prior1 <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 2.002),G3 = list(V = 1, nu = 2.002)), R = list(V = 1, nu = 2.002))
prime_mass_modS <- MCMCglmm(MassS ~ Age, random=~animal+Year+Code,rcov=~units,family="gaussian",prior=prior1,data=prime_mass,pedigree= pedigree,nitt=705000,thin=500,burnin=5000)

#### Antler Length ####
length <- DataForm[which(!is.na(DataForm$Length)),]
length$LengthS <- scale(length$Length)
prior1 <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 2.002),G3 = list(V = 1, nu = 2.002)), R = list(V = 1, nu = 2.002))
length_modS <- MCMCglmm(LengthS ~ Age + I(Age^2) + I(Age^3) + LastSeen, random=~animal+Year+Code,rcov=~units,family="gaussian",prior=prior1,data=length,pedigree= pedigree,nitt=1505000,thin=500,burnin=5000)

young_length <- length[which(length$Age<8),]
young_length$LengthS <- scale(young_length$Length)
prior1 <- list(G = list(G1 = list(V = 1, nu = 3.002, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, nu = 2.002),G3 = list(V = 1, nu = 2.002)), R = list(V = 1, nu = 2.002))
young_length_modS <- MCMCglmm(LengthS ~ Age, random=~animal+Year+Code,rcov=~units,family="gaussian",prior=prior1,data=young_length,pedigree= pedigree,nitt=505000,thin=500,burnin=5000)

## Trivariate scaled
prior3 <- list(G = list(G1 = list(V = diag(3), nu = 1.002, alpha.mu = rep(0,3), alpha.V = diag(3)*1000),G2 = list(V = diag(3), nu = 2.002),G3 = list(V = diag(3), nu = 2.002)), R = list(V = diag(3), nu = 2.002))
tri_modS <- MCMCglmm(cbind(MassS,LengthS,PointsS) ~ AgeAntlerGrown + I(AgeAntlerGrown^2) + I(AgeAntlerGrown^3), random = ~us(trait):animal + us(trait):YearAntlerGrown + us(trait):Code, rcov = ~idh(trait):units, family = c("gaussian","gaussian","gaussian"), prior=prior3, data=antlers, pedigree = ped_ordered,nitt=1505000,thin=500,burnin=5000)#)#,nitt=10,burnin=1,thin=1)#)#
