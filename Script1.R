####Charger données
db=read.csv("./g98_5years.csv",sep=",")



###Créations des dummies schooling
db=db[!is.na(db$sch_lev),]
db$sch_lev= factor(db$sch_lev)
dummies = data.frame(model.matrix(~db$sch_lev+0))
db=cbind(db,dummies)


### Transform monthly earnings
db$wage_36=db$wage_36*12
db$wage_60=db$wage_60*12
db$wage_first=db$wage_first*12
