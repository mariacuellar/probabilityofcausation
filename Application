# # Estimating the probability of causation for the RCT described (and with data) here 
# # https://www.povertyactionlab.org/evaluation/cleaning-springs-kenya
# test


# -----------------------
# Packages
# -----------------------

library(readstata13)
expit = function(x){ return( exp(x)/(1+exp(x)) ) }


# -----------------------
# Load and preselect data
# -----------------------

# Working directory
setwd("/Users/mariacuellar/Documents/CMU/Papers/2nd Heinz paper:ADA/Shaken Baby Syndrome/THESIS/Data leads/JPAL/Springs/dta/")

# Load data
dat1 = read.dta13("reg_data_children_Aug2010.dta", missing.type = TRUE)
dim(dat1)

# Remove the cases Kremer et al. removed for their article
dat2 = dat1[which( dat1$multiusers_l_base==0 | is.na(dat1$multiusers_l_base==0) ),]
dat3 = dat2[which( dat2$height_outlier_severe==0 | is.na(dat2$height_outlier_severe) ),]
dat4 = dat3[which( dat3$problem_weight==0 | is.na(dat3$problem_weight) ),]
dat5 = dat4[which( dat4$problem_bmi==0 | is.na(dat4$problem_bmi) ),]
dat6 = dat5[which( dat5$flag_age==0 | is.na(dat5$flag_age) ),]
dat7 = dat6[which( dat6$problem_age==0 | is.na(dat6$problem_age) ),]
dat8 = dat7[which( dat7$base_age<=3 | is.na(dat7$base_age) ),]

dat = dat8
dim(dat)

# Make my variables
a = ifelse(dat$evertreat=="TREAT", 0, 1) # Treatment. What about this var? ba_tc
y = dat$c14_d_child_diarrhea # Diarrhea in past week. Diarrhea defined as three or more “looser than normal” stools within 24 hours at any time in the past week
x1 = dat$c13_c_child_gender # gender
x2 = dat$base_age # age at baseline
x3 = dat$momeduc_orig # Mother's years of education
x4.0 = dat$splnecmpn_base # Baseline water quality, ln(spring water E. coli MPN) # High quality water: MPN <=1, High or moderate quality: MPN < 126, water is poor quality: MPN = 126-1000
x5 = dat$e1_iron_roof_base # Home has iron roof indicator
x6.0 = dat$hygiene_know_base # Mother's hygiene knowledge at baseline. average of demeaned sum of number of correct responses given to the open-ended question “to your knowledge, what can be done to prevent diarrhea?
x7.0 = dat$latrine_density_base # Baseline latrine density
x8.0 = dat$numkids_base # Number of children under 12 living at home
samplesize = length(A)

# Undoing the de-meaned variables
x4.1 = x4.0 - (min(x4.0[x4.0<0], na.rm = TRUE) + max(x4.0[x4.0<0], na.rm = TRUE))
x4 = x4.1 - min(x4.1, na.rm = TRUE)

x6.1 = x6.0 - (min(x6.0[x6.0<0], na.rm = TRUE) + max(x6.0[x6.0<0], na.rm = TRUE))
x6 = x6.1 - min(x6.1, na.rm = TRUE) + 1

x7.1 = x7.0 - (min(x7.0[x7.0<0], na.rm = TRUE) + max(x7.0[x7.0<0], na.rm = TRUE))
x7 = x7.1 - min(x7.1, na.rm = TRUE)

x8.1 = x8.0 - (min(x8.0[x8.0<0], na.rm = TRUE) + max(x8.0[x8.0<0], na.rm = TRUE))
x8 = x8.1 - min(x8.1, na.rm = TRUE)

# Make data frame, delete all missing values
dff0 = as.data.frame(cbind(a, y, x1, x2, x3, x4, x5, x6, x7, x8))
dff = na.omit(dff0)

x = as.matrix(cbind(x1, x2, x3, x4, x5, x6, x7, x8))
n = dim(dff)[1]

nrow(dff) # [1] 2933





# -----------------------
# Checking results for ATE, from Ed's R package
# -----------------------

library0 = c("SL.earth","SL.gam","SL.glm","SL.glmnet", "SL.glm.interaction", "SL.mean","SL.ranger")
library1 = c("SL.gam","SL.glm")

ATE = ate(y=dff$y, a=dff$a, x=data.frame(cbind(dff[,3:10])))
ATE[1]$res
1 - ATE[1]$res$est[1]/ATE[1]$res$est[2] # [1] 1.465586



# -----------------------
# Example point for my analysis
# -----------------------

# selected point to evaluate PC
# 1) gender	2) age	3) mom educ	4) water quality	5) iron roof ind	6) hygiene knowledge	7) latrine density	8) num kids in householf
y = dff$y
a = dff$a
x=as.data.frame(cbind(dff[,3:10]))
xtest=x[1,] ; xtest[1,1:8] <- c(1, 6, 6, 4, 1, 3, .4, 5) 
nsplits = 1
start.list=c(rep(1,8))
tracetf=TRUE


# parametric estimation of nuisance parameters
out.p <- pcausation(y=y, a=a, x=x, xtest=xtest, nsplits=1, start.list=c(rep(1,8)), 
                  tracetf=FALSE, printres=TRUE, sl.lib = "SL.glm")

out.p$res.pi
# est.pi
# 1 0.6781117

out.p$res.if
# est.pc      se.if  ci.ll.if  ci.ul.if
# 2 0.1207903 0.04094849 0.1125228 0.1295765



# nonparametric estimation of nuisance parameters
out.np <- pcausation(y=y, a=a, x=x, xtest=xtest, nsplits=1, start.list=c(rep(1,8)), 
                       tracetf=FALSE, printres=TRUE, sl.lib = "SL.ranger")

out.np$res.pi
# est.pi
# 1 0.2382169

out.np$res.if
# est.pc      se.if  ci.ll.if  ci.ul.if
# 2 0.1207903 0.04094849 0.1125228 0.1295765







# -----------------------
# Exploratory data analysis
# -----------------------

# Checking crosstab for Y and A
table(A, useNA="ifany")
table(Y, useNA="ifany")
table(dff$A, dff$Y, useNA="ifany")

# Gender distribution
table(dff$X1)

# Iron roof indicator
table(dff$X5)

# Lots of histograms for conginuous (or categorical) variables
pdf(file="EDA.pdf", width=10, height=5)

par(mfrow=c(2,3))
# Age distribution
hist(dff$X2, col="gray", main="Age distribution", xlab="Age in months")

# Mother's education
hist(dff$X3, col="gray", main="Mother education distribution \n at baseline", xlab="Years of education")

# Water quality
hist(dff$X4, col="gray", main="Water quality at baseline", xlab="Water quality (ln(E.coli MPN))")

# Hygiene knowledge at baseline
hist(dff$X6, col="gray", main="Mother's hygiene knowledge \n at baseline", xlab="Hygiene knowledge score")

# Latrine density at baseline
hist(dff$X7, col="gray", main="Latrine density near household \n at baseline", xlab="Latrine density (out of 0.6)")

# Number of children at baseline
hist(dff$X8, col="gray", main="Number of children in household \n at baseline", xlab="Number of children")

par(mfrow=c(1,1))

dev.off()







# -----------------------
# Exploring the results for the coefficients
# -----------------------
# confidence intervals for betas 
coef.estimates <- round(out.np$res.coefs, 2)
coef.estimates
# Estimate Robust.SE z.val p.val ci.ll ci.ul
# b1     0.17      0.09  1.99  0.05  0.00  0.34
# b2    -0.21      0.04 -4.93  0.00 -0.30 -0.13
# b3    -0.03      0.01 -2.07  0.04 -0.05  0.00
# b4    -0.06      0.02 -2.81  0.00 -0.10 -0.02
# b5    -0.05      0.09 -0.57  0.57 -0.24  0.13
# b6    -0.08      0.02 -4.01  0.00 -0.12 -0.04
# b7    -0.28      0.24 -1.15  0.25 -0.76  0.20
# b8    -0.02      0.02 -1.15  0.25 -0.05  0.01

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


coef.est.tab <- coef.estimates[,1:4]
coef.est.tab$. <- c("", "***", "*", "***", "", "***", "*", "*")
rownames(coef.est.tab) <- c("Gender (m:0, f:1)", "Age (in months)", "Mother's years of educ.", 
                           "Water quality at spring", "Iron roof indicator", "Mother's hygiene knowledge", 
                           "Latrine density near household", "Diarrhea prevention score")

# produce table for latex
xtable(coef.est.tab)





# How does PC change for different values of covariates?

# 1) gender	2) age	3) mom educ	4) water quality	5) iron roof ind	6) hygiene knowledge	7) latrine density	8) num kids in householf

# PC vs. mother's years of education
varyinggammahat = c()
momeduc_range = seq(0,12,1)
for(i in momeduc_range){
  apoint = c(1, 1.5, i, 4, 1, 3, .4, 4)
  gammahat_apoint = expit(betahat[1]*apoint[1] + betahat[2]*apoint[2] + 
                            betahat[3]*apoint[3] + betahat[4]*apoint[4] +
                            betahat[5]*apoint[5] + betahat[6]*apoint[6] +
                            betahat[7]*apoint[7] + betahat[8]*apoint[8])
  varyinggammahat[i] = gammahat_apoint
}

setwd("/Users/mariacuellar/Desktop/")
pdf(file="pcvsmothereduc1.pdf", width=5, height=4)
plot(varyinggammahat, type="o", ylim=c(0,.3),
     main="Estimated PC v. mother's years of education", ylab="Estimated PC", xlab="Mother's years of education")
dev.off()





# PC vs. age
varyinggammahat1 = c()
age_range = seq(0,17,1)
for(i in age_range){
  apoint = c(1, i, 6, 4, 1, 3, .4, 4)
  gammahat_apoint = expit(betahat[1]*apoint[1] + betahat[2]*apoint[2] + 
                            betahat[3]*apoint[3] + betahat[4]*apoint[4] +
                            betahat[5]*apoint[5] + betahat[6]*apoint[6] +
                            betahat[7]*apoint[7] + betahat[8]*apoint[8])
  varyinggammahat1[i] = gammahat_apoint
}
varyinggammahat1

pdf(file="pcvsage1.pdf", width=5, height=4)
plot(varyinggammahat1, type="o", ylim=c(0,.3),
     main="Estimated PC v. age", ylab="Estimated PC", xlab="Age")
dev.off()



# PC v. gender

apoint = c(0, 1.5, 6, 4, 1, 3, .4, 4)
gammahat_male = expit(betahat[1]*apoint[1] + betahat[2]*apoint[2] + 
                            betahat[3]*apoint[3] + betahat[4]*apoint[4] +
                            betahat[5]*apoint[5] + betahat[6]*apoint[6] +
                            betahat[7]*apoint[7] + betahat[8]*apoint[8])

apoint = c(1, 1.5, 6, 4, 1, 3, .4, 4)
gammahat_female = expit(betahat[1]*apoint[1] + betahat[2]*apoint[2] + 
                            betahat[3]*apoint[3] + betahat[4]*apoint[4] +
                            betahat[5]*apoint[5] + betahat[6]*apoint[6] +
                            betahat[7]*apoint[7] + betahat[8]*apoint[8])

gammahat_male
gammahat_female










# Trying out new ways of fitting my IF model

hist(IF_N_X_ystar)

IF_N_X_ystar_norm = IF_N_X_ystar
IF_N_X_ystar_norm = (IF_N_X_ystar_norm-min(IF_N_X_ystar_norm))
IF_N_X_ystar_norm = IF_N_X_ystar_norm/max(IF_N_X_ystar_norm)
hist(IF_N_X_ystar_norm)

logit = function(x){log(x/(1-x))}

mod_norm = nls(formula = IF_N_X_ystar~logit(X1+X2+X3+X4+X5+X6+X7+X8), start = rep(0,8), data=dff)

IF_N_X_model = nls(IF_N_X_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
                   start=list(beta1=0.1, beta2=0.1, beta3=0.1, beta4=0.1),
                   lower=rep(-2,4), upper=rep(2,4), algorithm="port",
                   data=dff, nls.control(maxiter = 500))

IF_N_X_model
newyhat = predict(mod_norm, newdata = dff.test1, type="response")
newyhat

