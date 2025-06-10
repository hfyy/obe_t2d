install.packages("msm")
install.packages("elect")

library(msm)
library(elect)
library(haven)
library(dplyr)
library(ggplot2)
library(tibble)

data <- read_dta("/Users/yiyue/Library/CloudStorage/OneDrive-CUHK-Shenzhen/T2D_obesity/data/sample_t2dobe.dta")
data <- data %>% select(c("wave","raeduc","hhidpn","rbmi","rdiab_raw","race","age","gender","rdiab","robe","dead","hhidpn",
                          "PGS_BMI","age_base","ec2_std","PGS_T2D","wave","year"))
data <- data %>% mutate(age = age-45)
data <- data %>% mutate(edu = case_when(raeduc == 1 ~ 0,
                                        raeduc == 2 ~ 0, 
                                        raeduc == 3 ~1 ))

data <- data %>% mutate(state = case_when(rdiab_raw ==0 & robe ==0 & dead==0~ 1, ## (1) no obese or T2D
                                          rdiab_raw ==1 & robe==0  & dead==0  ~2,    ## (2) T2D, no obese
                                          rdiab_raw ==0 & robe==1  & dead==0 ~ 3,   ## (3) obese, no T2D
                                          rdiab_raw ==1 & robe==1  & dead==0 ~ 4, ## (4) T2D, obese    
                                          dead ==1 ~ 5))  %>%   
                 arrange(hhidpn, wave)


statetable.msm(state, hhidpn, data=data)

any(!is.finite(data$state))
any(!is.finite(data$year))

data <- data[is.finite(data$state) , ]


#################################################################################
##               baseline
#################################################################################

Q  <-  rbind ( c(0, 0.02,0.06,0,0.02),
               c(0.04, 0, 0,0.05,0.04),
               c(0.2,0, 0, 0.06,0.02),
               c(0,0.12,0.05, 0,0.05),
               c(0,0,0,0,0))

msm <- msm( state ~ age, subject=hhidpn, data = data, death = 5, qmatrix = Q, 
            control=list(fnscale=10000,maxit=5000),  center=FALSE)
   msm
 

#################################################################################
##                model with covariates and interactions                               
################################################################################# 
 
msm2 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
               covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D+ ec2_std + factor(edu) )
   msm2

     
  pmatrix <- data.frame(pmatrix.msm(msm2, covariates = list(age=0))) 
  
   pmatrix.msm(msm2, covariates = list(age=2))
   pmatrix.msm(msm2, covariates = list(age=4))
   pmatrix.msm(msm2, covariates = list(age=6))
   pmatrix.msm(msm2, covariates = list(age=8))
   pmatrix.msm(msm2, covariates = list(age=10))
   pmatrix.msm(msm2, covariates = list(age=12))
   pmatrix.msm(msm2, covariates = list(age=14))
   pmatrix.msm(msm2, covariates = list(age=16))
   pmatrix.msm(msm2, covariates = list(age=18))
   pmatrix.msm(msm2, covariates = list(age=20))
   pmatrix.msm(msm2, covariates = list(age=22))
   pmatrix.msm(msm2, covariates = list(age=24))
   pmatrix.msm(msm2, covariates = list(age=26))
   pmatrix.msm(msm2, covariates = list(age=28))
   pmatrix.msm(msm2, covariates = list(age=30))
   pmatrix.msm(msm2, covariates = list(age=32))
   pmatrix.msm(msm2, covariates = list(age=34))
   pmatrix.msm(msm2, covariates = list(age=36))
   pmatrix.msm(msm2, covariates = list(age=38))
   pmatrix.msm(msm2, covariates = list(age=40))
              
  
msm2_inter <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
                   covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D + ec2_std + ec2_std*PGS_BMI + factor(edu) )

msm2_inter2 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
                   covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D + ec2_std + ec2_std*PGS_T2D + factor(edu) )



#################################################################################
##               hidden Markov model                  
################################################################################# 

Q  <-  rbind ( c(0, 0.01,0.03,0,0.1),
               c(0, 0, 0,0.05,0.04),
               c(0.1,0, 0, 0.03,0.01),
               c(0,0.06,0, 0,0.02),
               c(0,0,0,0,0))

ematrix <- rbind(c(0, 0.001, 0.001, 0.001,0),
                 c(0.5, 0, 0.001, 0.001,0),
                 c(0.001, 0.5, 0, 0.001,0),
                 c(0.5, 0.001, 0.5, 0,0),
                 c(0,0,0,0,0))


msm3 <- msm( state ~ time, subject=hhidpn, data = data, death = 5, qmatrix = Q, ematrix = ematrix,covariates = ~age, 
            control=list(fnscale=10000,maxit=5000,reltol= 1e-08 ),  center=FALSE, method = "BFGS")
 msm3

msm4 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, ematrix = ematrix, death=5, control=list(fnscale=5000,maxit=500, reltol = 1e-08),method = "BFGS", 
             covariates= ~ age + factor(gender) + PGS_BMI + PGS_T2D + ec2_std + factor(edu))
 msm4
  