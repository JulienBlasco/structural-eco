
####Charger données
db=read.csv("/Users/paul-armand.veillon/Desktop/Struct_Eco/9f3c724eb6559842417329162d4acab0/g98_5years.csv",sep=",")

db=read.csv("/home/odran/Dropbox/ENSAE/Structural econometrics/g98_5years.csv",sep=",")

db=read.csv("./structural-eco/g98_5years.csv",sep=",")


###Créations des dummies schooling
db=db[!is.na(db$sch_lev),]
db$sch_lev= factor(db$sch_lev)
dummies = data.frame(model.matrix(~db$sch_lev+0))
db=cbind(db,dummies)


### Transform monthly earnings
db$wage_36=db$wage_36*12
db$wage_60=db$wage_60*12
db$wage_first=db$wage_first*12

### get the log of it

db$lnwage_36=log(db$wage_36)
db$lnwage_60=log(db$wage_60)
db$lnwage_first=log(db$wage_first)

### Create As

s = c(db$sch_lev)
A = 16+(s-1)

### Calculate the phi of the wage regression

data = data.frame(db[,c("numenq","wage_first","wage_36","wage_60")])
                     
wide_format = reshape(db,idvar="numenq",varying=c("lnwage_first","lnwage_36","lnwage_60"),
        v.names = "lnwage",
        timevar = "time",
        times = c("lnwage_first","lnwage_36","lnwage_60"),
        direction="long")

wide_format$t_A    = 0
wide_format$t_A    = ifelse(wide_format$time=="lnwage_36",3,wide_format$t_A )
wide_format$t_A    = ifelse(wide_format$time=="lnwage_60",5,wide_format$t_A )
wide_format$t_Asquare = wide_format$t_A^2   

coef_wage = lm(wide_format$lnwage~wide_format$t_A+wide_format$t_Asquare+wide_format$s+
     wide_format$late6)

phi = coef_wage$coefficients
phis = c(0,coef_wage$coefficients[c("wide_format$s2","wide_format$s3",
          "wide_format$s4","wide_format$s5","wide_format$s6",
          "wide_format$s7","wide_format$s8")])

### Selection de données

S = 10000
Selection_variable = c(3:20,22:23)

db=db[rowSums(is.na(db[,Selection_variable]))==0,]

db = db[sample(1:dim(db)[1],S,replace=FALSE),]

s = c(db$sch_lev)

### Calculate the phi of Vw(s,As)

A1 = 16
A2 = 17
A3 = 18
A4 = 19
A5 = 20
A6 = 21
A7 = 22
A8 = 23
T = 50


somme_geometrique = function(q,n){
  (1-q^(n+1))/(1-q)
}

somme_geometrique_derivee_premiere = function(q,n){
  (n*(q^(n+1))-(n+1)*(q^n)+1)/((1-q)^2)
}

somme_geometrique_derivee_seconde = function(q,n){
  (-n*(n-1)*(q^(n+1))+2*((n^2)-1)*(q^n)-n*(n+1)*q^(n-1)+2)/((1-q)^3)
  }

Vw_s_As = function(beta,s,As){
  somme_geometrique(beta,T-As)*(phi["(Intercept)"]+phis[s]+
                                  phi["wide_format$late6"]*db$late6)+
    phi["wide_format$t_A"]*beta*(somme_geometrique_derivee_premiere(beta,T-As))+
    phi["wide_format$t_Asquare"]*(beta^2)*(somme_geometrique_derivee_seconde(beta,T-As))+
    phi["wide_format$t_Asquare"]*beta*(somme_geometrique_derivee_premiere(beta,T-As))
}

Vw_s_As(beta,1,A1)

### Calculate the probabilities

### It is said that epsilon should that epsilon should follow
### a standard normal. Then it means looking at the utility
### that heterogeneities in school taste is fixed. It seems
### strange and unnecessary to me. We don't have the problem
### of identifying the scale parameter here since we have got
### real units for the log wages.

###############################################################################
#######################    Vraisemblance   ####################################
###############################################################################

#### Covariates

X = as.matrix(cbind(db[,Selection_variable]))

vraisemblance=function(parametres){

beta     = parametres[1]
delta_s  = parametres[2:9]
delta_l  = parametres[10]
delta_x  = cbind(parametres[11:(11+dim(X)[2]-1)])
  

## Shouldn't we add a beta here in front of Vw_s_As(beta,8,A8)?

c8 = Vw_s_As(beta,7,A7)-Vw_s_As(beta,8,A8)-delta_s[8]-delta_l*db$late6-c(X%*%delta_x)
  
Pd8 = 1 - pnorm(c8)

EV8condV7 = delta_s[8]+delta_l*db$late6+Vw_s_As(beta,8,A8)+dnorm(c8)/(1-pnorm(c8))

EV8 = Pd8 * EV8condV7 + (1-Pd8)*Vw_s_As(beta,7,A7)

## Shouldn't we add a beta here in front of the expectation?

c7 = Vw_s_As(beta,6,A6)-EV8-delta_s[7]-delta_l*db$late6-c(X%*%delta_x)

Pd7 = 1 - pnorm(c7)

EV7condV6 = delta_s[7]+delta_l*db$late6+EV8+dnorm(c7)/(1-pnorm(c7))

EV7 = Pd7 * EV7condV6 + (1-Pd7)*Vw_s_As(beta,6,A6)

## Shouldn't we add a beta here in front of the expectation?

c6 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd6 = 1 - pnorm(c6)

EV6condV5 = delta_s[6]+delta_l*db$late6+EV7+dnorm(c6)/(1-pnorm(c6))

EV6 = Pd6 * EV6condV5 + (1-Pd6)*Vw_s_As(beta,5,A5)

## Shouldn't we add a beta here in front of the expectation?

c5 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd5 = 1 - pnorm(c5)

EV5condV4 = delta_s[5]+delta_l*db$late6+EV6+dnorm(c5)/(1-pnorm(c5))

EV5 = Pd5 * EV5condV4 + (1-Pd5)*Vw_s_As(beta,4,A4)

## Shouldn't we add a beta here in front of the expectation?

c4 = Vw_s_As(beta,4,A4)-EV6-delta_s[5]-delta_l*db$late6-c(X%*%delta_x)

Pd4 = 1 - pnorm(c4)

EV4condV3 = delta_s[4]+delta_l*db$late6+EV5+dnorm(c4)/(1-pnorm(c4))

EV4 = Pd4 * EV4condV3 + (1-Pd4)*Vw_s_As(beta,3,A3)

## Shouldn't we add a beta here in front of the expectation?

c3 = Vw_s_As(beta,3,A3)-EV7-delta_s[4]-delta_l*db$late6-c(X%*%delta_x)

Pd3 = 1 - pnorm(c3)

EV3condV2 = delta_s[3]+delta_l*db$late6+EV4+dnorm(c3)/(1-pnorm(c3))

EV3 = Pd3 * EV3condV2 + (1-Pd3)*Vw_s_As(beta,2,A2)

## Shouldn't we add a beta here in front of the expectation?

c2 = Vw_s_As(beta,2,A2)-EV7-delta_s[3]-delta_l*db$late6-c(X%*%delta_x)

Pd2 = 1 - pnorm(c2)

############
############
############

likelihood_ind =  (1-Pd2)*(s==1)+
    Pd2*(1-Pd3)*(s==2)+
    Pd2*Pd3*(1-Pd4)*(s==3)+
    Pd2*Pd3*Pd4*(1-Pd5)*(s==4)+
    Pd2*Pd3*Pd4*Pd5*(1-Pd6)*(s==5)+
    Pd2*Pd3*Pd4*Pd5*Pd6*(1-Pd7)*(s==6)+
    Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*(1-Pd8)*(s==7)+
    Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*Pd8*(s==8)

likelihood = sum(log(likelihood_ind))

return(-likelihood)

}

start= c(0.9,rep(-2,8),0,rep(0,dim(X)[2]))

resultat=optim(start,vraisemblance,control=list(maxit=20000))

coefficient = resultat$par

names(coefficient) = c(rep(0,10),names(db[,Selection_variable]))

####################################################################
####################################################################  
####################################################################  
##### Question 1 : prédiction du niveau de schooling ###############
####################################################################  
####################################################################  
####################################################################

### On calcule les probabilités conditionnelles, Pdi où i ={1...8}

parametres = coefficient

beta     = parametres[1]
delta_s  = parametres[2:9]
delta_l  = parametres[10]
delta_x  = cbind(parametres[11:(11+dim(X)[2]-1)])


## Shouldn't we add a beta here in front of Vw_s_As(beta,8,A8)?

c8 = Vw_s_As(beta,7,A7)-Vw_s_As(beta,8,A8)-delta_s[8]-delta_l*db$late6-c(X%*%delta_x)

Pd8 = 1 - pnorm(c8)

EV8condV7 = delta_s[8]+delta_l*db$late6+Vw_s_As(beta,8,A8)+dnorm(c8)/(1-pnorm(c8))

EV8 = Pd8 * EV8condV7 + (1-Pd8)*Vw_s_As(beta,7,A7)

## Shouldn't we add a beta here in front of the expectation?

c7 = Vw_s_As(beta,6,A6)-EV8-delta_s[7]-delta_l*db$late6-c(X%*%delta_x)

Pd7 = 1 - pnorm(c7)

EV7condV6 = delta_s[7]+delta_l*db$late6+EV8+dnorm(c7)/(1-pnorm(c7))

EV7 = Pd7 * EV7condV6 + (1-Pd7)*Vw_s_As(beta,6,A6)

## Shouldn't we add a beta here in front of the expectation?

c6 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd6 = 1 - pnorm(c6)

EV6condV5 = delta_s[6]+delta_l*db$late6+EV7+dnorm(c6)/(1-pnorm(c6))

EV6 = Pd6 * EV6condV5 + (1-Pd6)*Vw_s_As(beta,5,A5)

## Shouldn't we add a beta here in front of the expectation?

c5 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd5 = 1 - pnorm(c5)

EV5condV4 = delta_s[5]+delta_l*db$late6+EV6+dnorm(c5)/(1-pnorm(c5))

EV5 = Pd5 * EV5condV4 + (1-Pd5)*Vw_s_As(beta,4,A4)

## Shouldn't we add a beta here in front of the expectation?

c4 = Vw_s_As(beta,4,A4)-EV6-delta_s[5]-delta_l*db$late6-c(X%*%delta_x)

Pd4 = 1 - pnorm(c4)

EV4condV3 = delta_s[4]+delta_l*db$late6+EV5+dnorm(c4)/(1-pnorm(c4))

EV4 = Pd4 * EV4condV3 + (1-Pd4)*Vw_s_As(beta,3,A3)

## Shouldn't we add a beta here in front of the expectation?

c3 = Vw_s_As(beta,3,A3)-EV7-delta_s[4]-delta_l*db$late6-c(X%*%delta_x)

Pd3 = 1 - pnorm(c3)

EV3condV2 = delta_s[3]+delta_l*db$late6+EV4+dnorm(c3)/(1-pnorm(c3))

EV3 = Pd3 * EV3condV2 + (1-Pd3)*Vw_s_As(beta,2,A2)

## Shouldn't we add a beta here in front of the expectation?

c2 = Vw_s_As(beta,2,A2)-EV7-delta_s[3]-delta_l*db$late6-c(X%*%delta_x)

Pd2 = 1 - pnorm(c2)

### On a les probabilités conditionnelles on calcul les
### Pi proba d'être à chaque niveau

P1 = 1-Pd2
P2 = Pd2*(1-Pd3)
P3 = Pd2*Pd3*(1-Pd4)
P4 = Pd2*Pd3*Pd4*(1-Pd5)
P5 = Pd2*Pd3*Pd4*Pd5*(1-Pd6)
P6 = Pd2*Pd3*Pd4*Pd5*Pd6*(1-Pd7)
P7 = Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*(1-Pd8)
P8 = Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*Pd8

P_predit = cbind(P1,P2,P3,P4,P5,P6,P7,P8)

## Nos prédictions comparées aux bons résultats

colMeans(P_predit)*10000
tabulate(s)

sum(colMeans(P_predit)*c(1:8))
mean(s)

####################################################################
####################################################################  
####################################################################  
##### Question 2 : effect of technological change ##################
####################################################################  
####################################################################  
####################################################################  

### Effect of education on log wages multiplied by 2
### cela veut dire que les phi_s sont multipliés par 2
### au vu de la question je pense qu'il ne fallait qu'un seul
### phi_s (faire un effet linéaire de l'éducation et non
### un truc par tranche mais bon, on pourrait le changer)
P_predit_vrai = P_predit
phis = 2*phis

### On calcule les probabilités conditionnelles, Pdi où i ={1...8}

parametres = coefficient

beta     = parametres[1]
delta_s  = parametres[2:9]
delta_l  = parametres[10]
delta_x  = cbind(parametres[11:(11+dim(X)[2]-1)])


## Shouldn't we add a beta here in front of Vw_s_As(beta,8,A8)?

c8 = Vw_s_As(beta,7,A7)-Vw_s_As(beta,8,A8)-delta_s[8]-delta_l*db$late6-c(X%*%delta_x)

Pd8 = 1 - pnorm(c8)

EV8condV7 = delta_s[8]+delta_l*db$late6+Vw_s_As(beta,8,A8)+dnorm(c8)/(1-pnorm(c8))

EV8 = Pd8 * EV8condV7 + (1-Pd8)*Vw_s_As(beta,7,A7)

## Shouldn't we add a beta here in front of the expectation?

c7 = Vw_s_As(beta,6,A6)-EV8-delta_s[7]-delta_l*db$late6-c(X%*%delta_x)

Pd7 = 1 - pnorm(c7)

EV7condV6 = delta_s[7]+delta_l*db$late6+EV8+dnorm(c7)/(1-pnorm(c7))

EV7 = Pd7 * EV7condV6 + (1-Pd7)*Vw_s_As(beta,6,A6)

## Shouldn't we add a beta here in front of the expectation?

c6 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd6 = 1 - pnorm(c6)

EV6condV5 = delta_s[6]+delta_l*db$late6+EV7+dnorm(c6)/(1-pnorm(c6))

EV6 = Pd6 * EV6condV5 + (1-Pd6)*Vw_s_As(beta,5,A5)

## Shouldn't we add a beta here in front of the expectation?

c5 = Vw_s_As(beta,5,A5)-EV7-delta_s[6]-delta_l*db$late6-c(X%*%delta_x)

Pd5 = 1 - pnorm(c5)

EV5condV4 = delta_s[5]+delta_l*db$late6+EV6+dnorm(c5)/(1-pnorm(c5))

EV5 = Pd5 * EV5condV4 + (1-Pd5)*Vw_s_As(beta,4,A4)

## Shouldn't we add a beta here in front of the expectation?

c4 = Vw_s_As(beta,4,A4)-EV6-delta_s[5]-delta_l*db$late6-c(X%*%delta_x)

Pd4 = 1 - pnorm(c4)

EV4condV3 = delta_s[4]+delta_l*db$late6+EV5+dnorm(c4)/(1-pnorm(c4))

EV4 = Pd4 * EV4condV3 + (1-Pd4)*Vw_s_As(beta,3,A3)

## Shouldn't we add a beta here in front of the expectation?

c3 = Vw_s_As(beta,3,A3)-EV7-delta_s[4]-delta_l*db$late6-c(X%*%delta_x)

Pd3 = 1 - pnorm(c3)

EV3condV2 = delta_s[3]+delta_l*db$late6+EV4+dnorm(c3)/(1-pnorm(c3))

EV3 = Pd3 * EV3condV2 + (1-Pd3)*Vw_s_As(beta,2,A2)

## Shouldn't we add a beta here in front of the expectation?

c2 = Vw_s_As(beta,2,A2)-EV7-delta_s[3]-delta_l*db$late6-c(X%*%delta_x)

Pd2 = 1 - pnorm(c2)

### On a les probabilités conditionnelles on calcul les
### Pi proba d'être à chaque niveau

P1 = 1-Pd2
P2 = Pd2*(1-Pd3)
P3 = Pd2*Pd3*(1-Pd4)
P4 = Pd2*Pd3*Pd4*(1-Pd5)
P5 = Pd2*Pd3*Pd4*Pd5*(1-Pd6)
P6 = Pd2*Pd3*Pd4*Pd5*Pd6*(1-Pd7)
P7 = Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*(1-Pd8)
P8 = Pd2*Pd3*Pd4*Pd5*Pd6*Pd7*Pd8

P_predit = cbind(P1,P2,P3,P4,P5,P6,P7,P8)

## Nos prédictions comparées aux bons résultats

colMeans(P_predit)*10000
tabulate(s)

sum(colMeans(P_predit)*c(1:8))
sum(colMeans(P_predit_vrai)*c(1:8))
mean(s)

#### Effet du changement technologique (everybody gets a Phd...)

sum(colMeans(P_predit)*c(1:8))-sum(colMeans(P_predit_vrai)*c(1:8))