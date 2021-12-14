#########------------------------- GLOBAL DATA AND PARAMETERS -------------------------#############

nSIMUL=30 #100
npop=15
EXPE<-c(103)
#EXPE<-c(103,303,503,703,903,1103)
nYears=40
nInit=10
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
viridis <- rev(hcl.colors(n=6))
gap <- c(-.2,-.1,0, .1, .2, .3)
dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
dat<-dat[-1,]
Type=dat$Type
Type[7]<-"source"

#########--------------------------- PACKAGES AND FUNCTIONS ---------------------------#############
library(visNetwork)
library(reshape2)

pow <- function(x,y) {x^y}

coefSURV=function(pGres, kappaRIV, sigRIV) {(exp(-kappaRIV*pow(pGres,sigRIV))-exp(-kappaRIV))/(1-exp(-kappaRIV))}

#########---------------------------- CODE FOR FIGURES ----------------------------#############

##############################################################
### Fig.2 A) Theoretical optimal value of growth potential ###
##############################################################

#Populations parameters
# Survival in River
#with river trade off
SP0=0.9841606*1.0025 
SP1=0.9914398*1.0025 
SP1S=0.9967923*1.002
SP1M=0.9863295*1.002
SPnM=0.9911798*1.002
SPn=0.99775*1.002

#without river trade off
SP02=0.9841606
SP12=0.9914398
SP1S2=0.9967923
SP1M2=0.9863295
SPnM2=0.9911798
SPn2=0.99775

b=0.31
dr=0.5 #0.427 #edit alamarins new environmental conditions
Tlr=6
gr=0.208
Tur=24.5
bdens=142.7
lwa=3.804
lwb=0.32
ppercfm=0.12

#with river trade off
maxRIV=5 # maxRIV maximum growth rate, aded for the incorporation of trade-off growth/survival
sigRIV=3.7 #10 # sigRIV shape of the trade-off function
kappaRIV=0.001 #921.6 # kappaRIV shape of the trade-off function

#without river trade off
maxRIV2=50
sigRIV2=100
kappaRIV2=0.001

kg=0.003057
wmax=8500
lwa_sea=3.82568
lwb_sea=0.333779
alphaS=2.533333
betaS=-0.524
aNegg = 0.86
bNegg = 1.63

#Model fitness for a range of pG values (pG = pGsea) - trade-off in river only

#individuals born 1st april
pG=sort(rnorm(100000,0,0.2)) ##careful here gG #if range of pG:sort(rnorm(100000,0,(0.2/sqrt(0.14))))
pGres=exp(pG)/maxRIV;
pGres=ifelse(pGres>1,1,pGres)

pGres2=exp(pG)/maxRIV2;
pGres2=ifelse(pGres2>1,1,pGres2)

winit=rep(0.1,100000) #initial weight
linit=rep(0.1,100000) #initial length
fatinit=rep(0.1,100000) #initial fat reserves

#growth april-october
#omega=dr*(15-Tlr)*(1-exp(gr*(15-Tur)))
omega=sum(rep(dr*(15-Tlr)*(1-exp(gr*(15-Tur))),183)) #temperature effect #15 #183
denseffect=.6 #density effect (adjusted here)
activity=exp(pG) #activity effect - summer:1
w_parr= (winit^b + b*omega*activity*denseffect/100)^(1/b)
fat_parr = fatinit + (w_parr-winit)*0.2 
l_parr = exp(lwa+lwb*log((w_parr-fat_parr)/(1-ppercfm)))
# plot(w_parr~pG,type='l', xlab="pG",ylab="weight after summer")
# plot(fat_parr~pG,type='l', xlab="pG",ylab="fat after summer")
# plot(l_parr~pG,type='l', xlab="pG",ylab="size after summer")

#survival april-october
nb_days=183#183
P_survival_parr=pow(SP0*coefSURV(pGres, kappaRIV, sigRIV),nb_days);
P_survival_parr2=pow(SP02*coefSURV(pGres2, kappaRIV2, sigRIV2),nb_days);

plot(pG,P_survival_parr,type='l',ylab="Survival",xlab="Growth potential",main="Survival summer")
lines(pG,P_survival_parr2, col="red")
legend("bottomright", cex=.8, legend=c("With TO","Without TO"), border=NA,fill=c("black","red"), bty="n")

#growth november-april with decision to migrate to sea
#omega=dr*(15-Tlr)*(1-exp(gr*(15-Tur)))
omega=sum(rep(dr*(11-Tlr)*(1-exp(gr*(11-Tur))),182)) #182
denseffect=.6
activity=exp(pG)*0.725 #smolt
w_smolt= (w_parr^b + b*omega*activity*denseffect/100)^(1/b)
fat_smolt = fat_parr + (w_smolt-w_parr)*0
l_smolt = exp(lwa+lwb*log((w_smolt-fat_smolt)/(1-ppercfm)))
# plot(w_smolt~pG,type='l', xlab="pG",ylab="weight after winter")
# plot(fat_smolt~pG,type='l', xlab="pG",ylab="fat after winter")
# plot(l_smolt~pG,type='l', xlab="pG",ylab="size after winter")

nb_days=182#182
P_survival_smolt=pow(SP1S*coefSURV(pGres, kappaRIV, sigRIV),nb_days);
P_survival_smolt2=pow(SP1S2*coefSURV(pGres2, kappaRIV2, sigRIV2),nb_days);

plot(pG,P_survival_smolt,type='l',ylab="Survival",xlab="Growth potential",main="Survival winter")
lines(pG,P_survival_smolt2,col="red")
legend("bottomright", cex=.8, legend=c("With TO","Without TO"), border=NA,fill=c("black","red"), bty="n")

P_survival_river = P_survival_parr*P_survival_smolt
P_survival_river2 = P_survival_parr2*P_survival_smolt2

plot(pG,P_survival_river,type='l',ylab="Survival rate",xlab="Growth potential",main="Survival during 1 year in river")
lines(pG,P_survival_river2,col="red")
legend("bottomright", cex=.8, legend=c("With TO","Without TO"), border=NA,fill=c("black","red"), bty="n")

#sea growth
w_previous= w_smolt
fat_previous= fat_smolt
SP_sea = 1
for (t in 1:(365)) {
  w_sea= w_previous + (kg * w_previous * log(exp(pG) * wmax / w_previous)*1) 
  fat_sea = fat_previous + (w_sea-w_previous)*0.15
  l_sea = exp(lwa_sea+lwb_sea*log((w_sea-fat_sea)/(1-ppercfm)))
  w_previous = w_sea
  fat_previous = fat_sea
  SP_sea = SP_sea * (1 - alphaS*(l_sea/exp(lwa_sea))^(betaS/lwb_sea))^(1/30) #no tradeoff in sea
  
}
# plot(w_sea~pG,type='l', xlab="pG",ylab="weight after winter")
# plot(fat_sea~pG,type='l', xlab="pG",ylab="fat after winter")
# plot(l_sea~pG,type='l', xlab="pG",ylab="size after winter")

#sea survival
SP_sea[SP_sea=="NaN"]<-0
plot(pG,SP_sea,type='l',ylab="Survival rate",xlab="Growth potential",main='Survival during 1 year in sea')

#Neggs female all reproduce
neggs = exp(aNegg * log(w_sea) + bNegg)
plot(pG,neggs,type='l',ylab="Neggs",xlab="Growth potential",main="Eggs number by female")

#total fitness
fitness= P_survival_river * SP_sea * neggs
fitness2= P_survival_river2 * SP_sea * neggs

plot(NULL,ylab="LRS (Lifetime Recruitment Success)",xlab="Growth potential", xlim=c(-1,1), ylim=c(0,5))
lines(pG[1:60000],fitness2[1:60000], lty=2, lwd=3)
lines(pG,fitness,lwd=3, lty=1)

legend("topright", cex=.8, legend=c("With Trade-off","Without Trade-off"), border=NA,lty=c(1,2), bty="n")


########################################################
### Fig.2 B) Dispersal kernel for 20% dispersal rate ###
########################################################

#dispersal laplace distribution
mu=0
beta=29.5

distance=matrix(data=NA, ncol=npop, nrow=npop)
colnames(distance)=rownames(distance)=pops

dist <- dat$Distance  #edit al - 20/03/21

for (i in 1:npop){
  for (j in 1:npop){
    distance[i,j]=abs(dist[j]-dist[i])
  }}

area_log=log10(dat$Area)

ratio_area=matrix(, nrow=1, ncol= npop); colnames(ratio_area)=pops
for (i in 1:npop){
  ratio_area[1,i]=area_log[i]/(sum(area_log))} 

lap = function(mu, beta, distance) {(1/(2*beta))*exp(-(distance-mu)/beta)} # Laplace


connect_kernel=matrix(, nrow=npop, ncol=npop)
rownames(connect_kernel)=colnames(connect_kernel)=pops
connect=list()

h=0.8

#eg Leff towards gradient of distance
distance<-seq(from=0, to=135, by=15)
attrac_d=lap(mu, beta, distance)
ratio_area1=rep(0.1,10) # small #big=0.08
attrac_ad=attrac_d*ratio_area1 # influence of area on attractivity
proba_attract_ad=attrac_ad/sum(attrac_ad) #proba of attractivity du site fonction distance at aire
sum(proba_attract_ad)
proba_dispersion=proba_attract_ad*(1-h)

plot(proba_dispersion~distance, ylim=c(0,0.1),pch=19,cex=1.3, xlab="Distance from donnor population (km)",ylab="Dispersal probability")
lines(proba_dispersion~distance)

distance<-seq(from=0, to=135, by=15)
attrac_d=lap(mu, beta, distance)
ratio_area2=c(rep(0.05,5),rep(0.15,5)) # small #big=0.08
attrac_ad2=attrac_d*ratio_area2 # influence of area on attractivity
proba_attract_ad2=attrac_ad2/sum(attrac_ad2) #proba of attractivity du site fonction distance at aire
sum(proba_attract_ad2)
proba_dispersion2=proba_attract_ad2*(1-h)

points(proba_dispersion2~distance,pch=17,cex=1.3)
lines(proba_dispersion2~distance, lty=2)
legend("topright",cex=1,title="Recipient populations",legend=c("Equal size","Variable size"),col="black",pch=c(19,17), border=NA, bty='n',lty=c(1,2))


########################################################################
### Fig.3 A) Ratio Immigrants/Emigrants and Proportion of Immigrants ###
########################################################################

nReturns_scn <- list()
Mig_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/DEMOGRAPHY",EXPE[iEXPE],".RData"))
  nReturns_scn[[iEXPE]]<-nReturns
  Mig_scn[[iEXPE]]<-Mig
}

I <- P.im <- array(,dim=c(nSIMUL,npop, length(EXPE)))
for (scn in 1:length(EXPE)){
  for (simul in 1:nSIMUL){
    ratio<-Mig_scn[[scn]][[simul]]$NIm/Mig_scn[[scn]][[simul]]$NEm
    ratio[which(ratio=="Inf")]<-NA
    I[simul, ,scn] <- apply(ratio[45:50,],2,median, na.rm=T) #5last years #mean , na.rm=T
    P.im[simul, ,scn] <- apply((Mig_scn[[scn]][[simul]]$NIm/nReturns_scn[[scn]][[simul]])[45:50,]*100,2,median, na.rm=T)
  }
}

par(mar = c(5, 4, 4, 4) + 0.3)
plot(NULL, xlim=c(0,16), ylim=c(0,4), xlab="Population",ylab="Median of Ratio Immigrants/Emigrants", xaxt='n')
mtext(pops, side = 1, line = 1, outer = FALSE, at = 1:15, cex=.9)

for (scn in 1:length(EXPE)){
  for (pop in 1:npop) {
    points(pop+gap[scn], median(I[,pop,scn],na.rm=TRUE), pch=20, cex=1.5,col=viridis[scn]);
  }
}
abline(h=1,col="black",lty=2)
#add Pim on same graph
par(new = TRUE)                             # Add new plot
plot(NULL, xlim=c(0,16), ylim=c(0,50), xlab="",ylab="", axes=FALSE)
axis(side = 4, at = c(0,10,20,30,40,50), col="red",col.axis="red")      # Add second axis
mtext("Median of Immigrants Proportion (%)", side = 4, line = 3, col="red") 
for (pop in 1:npop) {
  points(pop+-.1, median(P.im[,pop,2],na.rm=TRUE), pch=17, cex=1,col="red");
}

###########################################################
### Fig.3 B) Populations network for 10% dispersal rate ###
###########################################################

#Migrants flows for arrows
a <- array(,dim=c(npop,npop,nSIMUL, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      a[pop,,simul,scn]<-apply(Mig_scn[[scn]][[simul]]$Im[[pop]][46:50,],2,mean)
      for (pop2 in 1:npop) {
        if (a[pop,pop2,simul,scn]=="NaN" || is.na(a[pop,pop2,simul,scn])) {
          a[pop,pop2,simul,scn]=0
        }
      }
    }
  }
}
arr<-mflow<-list()
for (scn in 1:length(EXPE)) {
  arr[[scn]] <- array( unlist(a[,,,scn]) , c(15,15,100) )
  mflow[[scn]]<-apply( arr[[scn]] , 1:2 , mean )#mean simulations
  colnames(mflow[[scn]])<-pops
  rownames(mflow[[scn]])<-pops
}


scenario<-2 #SCENARIO CHOICE
mflow2<-melt(mflow[[scenario]]) #for absolute values 
mflow3<-mflow2[c(2,1,3)] #for absolute values

# mflow2<-t(mflow[[2]]) #for proportion
# for (l in 1:15) {
#  mflow2[l,]<-mflow2[l,]/sum(mflow2[l,]) #by proportion
# }
# mflow3<-melt(mflow2)

mflow3[,1]<-as.factor(mflow3[,1])
mflow3[,2]<-as.factor(mflow3[,2])

#populations size for nodes
Parpop<-array(dim=c(nSIMUL, npop, length(EXPE)))
Parpop2<-array(dim=c(npop, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      Parpop[simul, ,scn] <- colMeans(nReturns_scn[[scn]][[simul]][46:50,], na.rm=T) #5last years
    }
  }
  Parpop2[,scn]<-colMeans(Parpop[,,scn])
}

#network
nodes<-data.frame(id=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15"),
                  pop=pops,
                  type=c("1","2","1","2","1","3","2","3","2","1","2","1","2","1","3"), #control
                  #type=c("2","1","1","1","1","3","2","2","1","1","2","1","2","1","3"),
                  type.label=Type,
                  size=Parpop2[,scenario]#choice of scenarios
)
mflow3$Var2=as.character(mflow3$Var2)
mflow3$Var1=as.character(mflow3$Var1)

mflow3$Var2[which(mflow3$Var2=="Leff")]<-"p1" ;mflow3$Var1[which(mflow3$Var1=="Leff")]<-"p1"
mflow3$Var2[which(mflow3$Var2=="Trieux")]<-"p2" ;mflow3$Var1[which(mflow3$Var1=="Trieux")]<-"p2"
mflow3$Var2[which(mflow3$Var2=="Jaudy")]<-"p3" ;mflow3$Var1[which(mflow3$Var1=="Jaudy")]<-"p3"
mflow3$Var2[which(mflow3$Var2=="Leguer")]<-"p4" ;mflow3$Var1[which(mflow3$Var1=="Leguer")]<-"p4"
mflow3$Var2[which(mflow3$Var2=="Yar")]<-"p5" ;mflow3$Var1[which(mflow3$Var1=="Yar")]<-"p5"
mflow3$Var2[which(mflow3$Var2=="Douron")]<-"p6" ;mflow3$Var1[which(mflow3$Var1=="Douron")]<-"p6"
mflow3$Var2[which(mflow3$Var2=="Penze")]<-"p7" ;mflow3$Var1[which(mflow3$Var1=="Penze")]<-"p7"
mflow3$Var2[which(mflow3$Var2=="Elorn")]<-"p8" ;mflow3$Var1[which(mflow3$Var1=="Elorn")]<-"p8"
mflow3$Var2[which(mflow3$Var2=="Aulne")]<-"p9" ;mflow3$Var1[which(mflow3$Var1=="Aulne")]<-"p9"
mflow3$Var2[which(mflow3$Var2=="Goyen")]<-"p10" ;mflow3$Var1[which(mflow3$Var1=="Goyen")]<-"p10"
mflow3$Var2[which(mflow3$Var2=="Odet")]<-"p11" ;mflow3$Var1[which(mflow3$Var1=="Odet")]<-"p11"
mflow3$Var2[which(mflow3$Var2=="Aven")]<-"p12" ;mflow3$Var1[which(mflow3$Var1=="Aven")]<-"p12"
mflow3$Var2[which(mflow3$Var2=="Laita")]<-"p13" ;mflow3$Var1[which(mflow3$Var1=="Laita")]<-"p13"
mflow3$Var2[which(mflow3$Var2=="Scorff")]<-"p14" ;mflow3$Var1[which(mflow3$Var1=="Scorff")]<-"p14"
mflow3$Var2[which(mflow3$Var2=="Blavet")]<-"p15" ;mflow3$Var1[which(mflow3$Var1=="Blavet")]<-"p15"

links<-mflow3
colnames(links)<-c("from","to","weigth")
links<-links[-which(links$weigth==0),]
vis.nodes<-nodes
vis.links<-links

vis.nodes$shape  <- "dot"  
vis.nodes$shadow <- TRUE # Nodes will drop shadow
#vis.nodes$title  <- nodes$pop # Text on click
vis.nodes$label  <- nodes$pop # Node label
vis.nodes$size   <- nodes$size/5 # Node size
vis.nodes$borderWidth <- 2 # Node border width
vis.nodes$color.background <- c("#EE5C42", "#4F94CD", "#000000")[nodes$type]
vis.nodes$color.border <- "orange"
vis.nodes$color.highlight.background <- "red"
vis.nodes$color.highlight.border <- "gold"
tryagain<-links$weigth*2 #+1 for absolute values #*10 for proportion
vis.links$width <- tryagain # line width
vis.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
vis.links$smooth <- TRUE    # should the edges be curved?
vis.links$shadow <- FALSE    # edge shadow

#nodes (pop) geographical coordinates
lat<-c(48.701791,48.677231, 48.714110, 48.651108, 48.646443, 48.636448, 48.600963, 48.479118, 48.204225, 48.040924, 48.002719, 47.869874, 47.870947, 47.859305, 47.834972)
lon<-c(-3.056085, -3.157408 , -3.262510, -3.420931, -3.577629, -3.659678, -3.937004, -4.191071, -4.050434,-4.478585, -4.111829, -3.726740,  -3.545196, -3.400563, -3.207688)
#plot(lon, lat)
vis.nodes$x<- lon*1000
vis.nodes$y <- -lat*1000

a<-visNetwork(vis.nodes, vis.links)
a<-visEdges(a, arrows=list(to=list(enable=T, scaleFactor=1.5)),color = list(color = "black", highlight = "red"), smooth = list(enabled = TRUE, type = "diagonalCross"))
a<-visNodes(a,fixed = TRUE,physics=T, font=list(color="black", size=0))
a<-visOptions(a, highlightNearest = TRUE, selectedBy = "type.label")
a<-visLegend(a,addNodes = data.frame(label = c("Sink","Source","Neutral"), shape = c("dot"), 
                                     size = c(10,10,10), 
                                     color = c("red","blue","black"),font.size=20),
             addEdges = data.frame(label = "Emigration",color="black"), useGroups = FALSE)
a

###########################################################
############### Fig. 4 A) Portfolio effect ################
###########################################################

PE_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/DEMOGRAPHY",EXPE[iEXPE],".RData"))
  PE_scn[[iEXPE]]<-PE
}

plot(NULL, xlim=c(0,length(EXPE)+1), ylim=c(0,4), xlab="Dispersal rate (%)",ylab="Portfolio effect (PE)", xaxt='n')
for (scn in 1:length(EXPE)){
  points(scn-.2, median(PE_scn[[scn]][,"pe"]), pch=20, cex=2, col=viridis[scn])
  segments(scn-.2, quantile(PE_scn[[scn]][,"pe"],probs=.025),scn-.2,quantile(PE_scn[[scn]][,"pe"],probs=.975),lwd=1,col=viridis[scn])
  segments(scn-.2, quantile(PE_scn[[scn]][,"pe"],probs=.25),scn-.2,quantile(PE_scn[[scn]][,"pe"],probs=.75),lwd=2,col=viridis[scn])
}
#loess trend
nodiv<-array(,dim=c(nSIMUL*6,2))
for (scn in 1:6) {
  nodiv[(nSIMUL*(scn-1)+1):(nSIMUL*scn),1] <- PE_scn[[scn]][,"pe"]
  nodiv[(nSIMUL*(scn-1)+1):(nSIMUL*scn),2] <- rep(scn,nSIMUL)  
}
nodiv<-as.data.frame(nodiv)
colnames(nodiv) <- c("PE","Scn")
nodiv_loess <- loess(PE ~ Scn, nodiv)
lines(rep((1:6)-.2,each=100),nodiv_loess$fitted, col="black",lwd=3)
mtext(c("0","10","20","30","40","50"), side = 1, line = 1, outer = FALSE, at = 1:6,las=1)
abline(h=1, lty=2)

###########################################################
################### Fig. 4 B) Synchrony ###################
###########################################################

plot(NULL, xlim=c(0,length(EXPE)+1), ylim=c(0,.3), xlab="Dispersal rate (%)",ylab="Synchrony", xaxt='n')
for (scn in 1:length(EXPE)){
  points(scn-.2, median(PE_scn[[scn]][,"sync"]), pch=20, cex=2, col=viridis[scn])
  segments(scn-.2, quantile(PE_scn[[scn]][,"sync"],probs=.025),scn-.2,quantile(PE_scn[[scn]][,"sync"],probs=.975),lwd=1,col=viridis[scn])
  segments(scn-.2, quantile(PE_scn[[scn]][,"sync"],probs=.25),scn-.2,quantile(PE_scn[[scn]][,"sync"],probs=.75),lwd=2,col=viridis[scn])
}
#loess trend
nodiv<-array(,dim=c(nSIMUL*6,2))
for (scn in 1:6) {
  nodiv[(nSIMUL*(scn-1)+1):(nSIMUL*scn),1] <- PE_scn[[scn]][,"sync"]
  nodiv[(nSIMUL*(scn-1)+1):(nSIMUL*scn),2] <- rep(scn,nSIMUL)  
}
nodiv<-as.data.frame(nodiv)
colnames(nodiv) <- c("PE","Scn")
nodiv_loess <- loess(PE ~ Scn, nodiv)
lines(rep((1:6)-.2,each=nSIMUL),nodiv_loess$fitted, col="black",lwd=3)
mtext(c("0","10","20","30","40","50"), side = 1, line = 1, outer = FALSE, at = 1:6,las=1)

###########################################################
################# Fig. 4 C) CV Abundance ##################
###########################################################

Npops <- array(,dim=c(npop,nSIMUL, length(EXPE)))
for (scn in 1:length(EXPE)){
  for (simul in 1:nSIMUL){
    for (pop in 1:npop) {
      Npops[pop,simul,scn] <- sd(residuals(lm(nReturns_scn[[scn]][[simul]][10:50,pop] ~ c(10:50))))/mean(nReturns_scn[[scn]][[simul]][10:50,pop]) #nb returns metapop per year, simu and scenario
    }
  }
}

Npops_median <- array(,dim=c(npop, length(EXPE)))
for (scn in 1:length(EXPE)){
  for (pop in 1:npop) {
    Npops_median[pop,scn] <- median(Npops[pop,,scn],na.rm=T)
  }
}

colnames(Npops_median) <- c("0","10","20","30","40","50")
Npops_median_df<-as.data.frame(Npops_median)
Npops_median_df$Type<-dat$Type
Npops_median_df[7,7]<-"source"

plot(NULL, xlim=c(0.5,3), ylim=c(0,.5), xlab="",ylab="CV Abundance", xaxt='n') #15
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
#sink
sink<-Npops_median_df[Npops_median_df$Type=="sink",]
colnames(sink) <- c(0,10,20,30,40,50)
sink_melted<-melt(sink)
sink_melted$variable <- as.numeric(sink_melted$variable)
lw_sink<-loess(value ~ variable, sink_melted, span=1)
for (pop in which(Npops_median_df$Type=="sink")) {
  points(1+gap[1:6],Npops_median_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5)
#neutral
neutral<-Npops_median_df[Npops_median_df$Type=="neutral",]
colnames(neutral) <- c(0,10,20,30,40,50)
neutral_melted<-melt(neutral)
neutral_melted$variable <- as.numeric(neutral_melted$variable)
lw_neutral<-loess(value ~ variable, neutral_melted, span=1)
for (pop in which(Npops_median_df$Type=="neutral")) {
  points(1.7+gap[1:6],Npops_median_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1.7+gap[1:6],each=3),lw_neutral$fitted, col="black",lwd=5)
#source
source<-Npops_median_df[Npops_median_df$Type=="source",]
colnames(source) <- c(0,10,20,30,40,50)
source_melted<-melt(source)
source_melted$variable <- as.numeric(source_melted$variable)
lw_source<-loess(value ~ variable, source_melted, span=1)
for (pop in which(Npops_median_df$Type=="source")) {
  points(2.4+gap[1:6],Npops_median_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(2.4+gap[1:6],each=6),lw_source$fitted, col="black",lwd=5)


###########################################################
######## Fig. 4 D) Quasi-Extinction probability ###########
###########################################################

##Pextinction
Rmax=rep(10,npop)
Rmax5=Rmax*5/100
theta=Rmax5/100
Area=datapop$Area
rPROP=.25

nParr_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/DEMOGRAPHY",EXPE[iEXPE],".RData"))
  nParr_scn[[iEXPE]]<-nParr
}

density_parr<-array(dim=c(nYears+nInit, npop,nSIMUL, length(EXPE)))
density_parr_median<-array(dim=c(npop,nSIMUL, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      for (year in 1:(nYears+nInit)) {
        density_parr[year, pop,simul,scn] <- nParr_scn[[scn]][[simul]][year,pop]/(Area[pop]*rPROP)
      }
      density_parr_median[pop,simul,scn] <- median(density_parr[45:50,pop,simul,scn], na.rm=T)
    }
  }
}
density_parr_simul<-array(dim=c(npop, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (pop in 1:npop) {
    density_parr_simul[pop,scn] <- median(density_parr_median[pop,,scn], na.rm=T)
  }
}

density_parr_simul<-as.data.frame(density_parr_simul)
density_parr_simul$Type<-datapop$Type
density_parr_simul[7,7]<-"source"

## Supp Mat Density of parr
plot(NULL, xlim=c(0.5,3), ylim=c(0.05,0.3), xlab="",ylab="Density of parr (by m2)", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
#sink
sink<-density_parr_simul[density_parr_simul$Type=="sink",]
colnames(sink) <- c(0,10,20,30,40,50)
sink_melted<-melt(sink)
sink_melted$variable <- as.numeric(sink_melted$variable)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
for (pop in which(density_parr_simul$Type=="sink")) {
  points(1+gap[1:6],density_parr_simul[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5)
#neutral
neutral<-density_parr_simul[density_parr_simul$Type=="neutral",]
colnames(neutral) <- c(0,10,20,30,40,50)
neutral_melted<-melt(neutral)
neutral_melted$variable <- as.numeric(neutral_melted$variable)
lw_neutral<-loess(value ~ variable, neutral_melted,span=1)
for (pop in which(density_parr_simul$Type=="neutral")) {
  points(1.7+gap[1:6],density_parr_simul[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1.7+gap[1:6],each=3),lw_neutral$fitted, col="black",lwd=5)
#source
source<-density_parr_simul[density_parr_simul$Type=="source",]
colnames(source) <- c(0,10,20,30,40,50)
source_melted<-melt(source)
source_melted$variable <- as.numeric(source_melted$variable)
lw_source<-loess(value ~ variable, source_melted,span=1)
for (pop in which(density_parr_simul$Type=="source")) {
  points(2.4+gap[1:6],density_parr_simul[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(2.4+gap[1:6],each=6),lw_source$fitted, col="black",lwd=5)

### Pextinction considering 2 consecutive years
Ext <- array(,dim=c(nSIMUL,length(pops), length(EXPE)))
PQext <- array(,dim=c(length(pops), length(EXPE)));colnames(PQext) <- EXPE;rownames(PQext)<-pops
j=0
for (scn in which(names(nReturns) %in% EXPE)){
  j=j+1
  for (simul in 1:nSIMUL){
    for (pop in 1:npop) {
      tmp <- density_parr[,pop,simul,scn]
      c=0
      for (t in 1:((nYears+nInit)-1)) {
        if (tmp[t]<theta[pop] && tmp[t+1]<theta[pop]) {
          c=c+1
        }
      }
      if (c>=1) {
        Ext[simul,pop,j]<-1
      }
      else {
        Ext[simul,pop,j]<-0
      }
    }
  }
  PQext[,j] <- colSums(Ext[,,j])/nSIMUL
}  

PQext_df<-as.data.frame(PQext)
PQext_df$Type<-datapop$Type
PQext_df[7,7]<-"source"

plot(NULL, xlim=c(0.5,3), ylim=c(0,20), xlab="",ylab="Quasi-Extinction risk (%, threshold 5% Rmax)", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
#sink
sink<-PQext_df[PQext_df$Type=="sink",]
colnames(sink) <- c(0,10,20,30,40,50)
sink_melted<-melt(sink)
sink_melted$variable <- as.numeric(sink_melted$variable)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
for (pop in which(PQext_df$Type=="sink")) {
  points(1+gap[1:6],PQext_df[pop,1:6]*100, col=viridis[1:6], pch=20)
}
lines(rep(1+gap[1:6],each=6),lw_sink$fitted*100, col="black",lwd=5)
#neutral
neutral<-PQext_df[PQext_df$Type=="neutral",]
colnames(neutral) <- c(0,10,20,30,40,50)
neutral_melted<-melt(neutral)
neutral_melted$variable <- as.numeric(neutral_melted$variable)
lw_neutral<-loess(value ~ variable, neutral_melted,span=1)
for (pop in which(PQext_df$Type=="neutral")) {
  points(1.7+gap[1:6],PQext_df[pop,1:6]*100, col=viridis[1:6], pch=20)
}
lines(rep(1.7+gap[1:6],each=3),lw_neutral$fitted*100, col="black",lwd=5)
#source
source<-PQext_df[PQext_df$Type=="source",]
colnames(source) <- c(0,10,20,30,40,50)
source_melted<-melt(source)
source_melted$variable <- as.numeric(source_melted$variable)
lw_source<-loess(value ~ variable, source_melted,span=1)
for (pop in which(PQext_df$Type=="source")) {
  points(2.4+gap[1:6],PQext_df[pop,1:6]*100, col=viridis[1:6], pch=20)
}
lines(rep(2.4+gap[1:6],each=6),lw_source$fitted*100, col="black",lwd=5)


###########################################################
################ Fig. 5 A) Smolt size #####################
###########################################################

res_smolt_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/PHENOGENOTYPE",EXPE[iEXPE],".RData"))
  res_smolt_scn[[iEXPE]]<-res_smolt
}

size_smolt_scn<-list()
for (scn in 1:length(EXPE)){
  size_smolt<-list()
  for (pop in 1:npop) {
    size_smolt[[pop]] <- res_smolt_scn[[scn]][[1]][[pop]]$Lf[res_smolt_scn[[scn]][[1]][[pop]]$years>45]
    for (simul in 2:nSIMUL) {
      size_smolt[[pop]] <- c(size_smolt[[pop]], res_smolt_scn[[scn]][[simul]][[pop]]$Lf[res_smolt_scn[[scn]][[simul]][[pop]]$years>45])
    }
  }
  size_smolt_scn[[scn]] <- size_smolt
}

pop_group1<-c(1,3,5,10,12,14) #sink
pop_group2 <- c(6,8,15) #neutral
pop_group3 <- c(2,4,7,9,11,13) #source

group1<-array(,dim=c(length(pop_group1),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group1) {
    i=i+1
    group1[i,scn] <- median(size_smolt_scn[[scn]][[pop]])
  }
}
group2<-array(,dim=c(length(pop_group2),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group2) {
    i=i+1
    group2[i,scn] <- median(size_smolt_scn[[scn]][[pop]])
  }
}
group3<-array(,dim=c(length(pop_group3),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group3) {
    i=i+1
    group3[i,scn] <- median(size_smolt_scn[[scn]][[pop]])
  }
}

plot(NULL, xlim=c(0.5,3), ylim=c(130,145), xlab="Population",ylab="Smolt size (mm)", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
for (scn in 1:length(EXPE)) {
  points(rep(1+gap[scn],each=nrow(group1)), group1[,scn], col=viridis[scn],pch=20)
  points(rep(1.7+gap[scn],each=nrow(group2)), group2[,scn], col=viridis[scn],pch=20)
  points(rep(2.4+gap[scn],each=nrow(group3)), group3[,scn], col=viridis[scn],pch=20)
}

colnames(group1) <- c(0,10,20,30,40,50)
sink_melted<-melt(group1)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

colnames(group2) <- c(0,10,20,30,40,50)
sink_melted<-melt(group2)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1.7+gap[1:6],each=3),lw_sink$fitted, col="black",lwd=5) #3 #5

colnames(group3) <- c(0,10,20,30,40,50)
sink_melted<-melt(group3)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(2.4+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

###########################################################
################# Fig. 5 B) Adult 1SW size ################
###########################################################

res_1SW_scn <- list()
res_MSW_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/PHENOGENOTYPE",EXPE[iEXPE],".RData"))
  res_1SW_scn[[iEXPE]]<-res_1SW
  res_MSW_scn[[iEXPE]]<-res_MSW
}


size_1SW_scn<-list()
size_1SW_scn_residents<-list()
for (scn in 1:length(EXPE)) {
  size_1SW<-list()
  size_1SW_residents<-list()
  for (pop in 1:npop) {
    size_1SW[[pop]] <- res_1SW_scn[[scn]][[1]][[pop]]$Lf[res_1SW_scn[[scn]][[1]][[pop]]$years > 45]
    size_1SW_residents[[pop]] <- res_1SW_scn[[scn]][[1]][[pop]]$Lf[res_1SW_scn[[scn]][[1]][[pop]]$years > 45 & res_1SW_scn[[scn]][[1]][[pop]]$CollecID == pop]
    for (simul in 2:nSIMUL) {
      size_1SW[[pop]] <- c(size_1SW[[pop]], res_1SW_scn[[scn]][[simul]][[pop]]$Lf[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45])
      size_1SW_residents[[pop]] <- c(size_1SW_residents[[pop]], res_1SW_scn[[scn]][[simul]][[pop]]$Lf[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45 & res_1SW_scn[[scn]][[simul]][[pop]]$CollecID == pop])
    }
  }
  size_1SW_scn[[scn]] <- size_1SW
  size_1SW_scn_residents[[scn]] <- size_1SW_residents
}

pop_group1<-c(1,3,5,10,12,14) #sink
pop_group2 <- c(6,8,15) #neutral
pop_group3 <- c(2,4,7,9,11,13) #source

group1<-array(,dim=c(length(pop_group1),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group1) {
    i=i+1
    group1[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group1[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}
group2<-array(,dim=c(length(pop_group2),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group2) {
    i=i+1
    group2[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group2[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}
group3<-array(,dim=c(length(pop_group3),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group3) {
    i=i+1
    group3[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group3[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}

plot(NULL, xlim=c(0.5,3), ylim=c(595,615), xlab="Population",ylab="Adult 1SW size (mm)", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
for (scn in 1:length(EXPE)) {
  points(rep(1+gap[scn],each=nrow(group1)), group1[,scn], col=viridis[scn],pch=20)#color_transparent(col)
  points(rep(1.7+gap[scn],each=nrow(group2)), group2[,scn], col=viridis[scn],pch=20)
  points(rep(2.4+gap[scn],each=nrow(group3)), group3[,scn], col=viridis[scn],pch=20)
}

colnames(group1) <- c(0,10,20,30,40,50)
sink_melted<-melt(group1)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

colnames(group2) <- c(0,10,20,30,40,50)
sink_melted<-melt(group2)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1.7+gap[1:6],each=3),lw_sink$fitted, col="black",lwd=5) #3 #5

colnames(group3) <- c(0,10,20,30,40,50)
sink_melted<-melt(group3)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(2.4+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5


####################################################################
#### Fig. 5 C) Adult 1SW genetic male sea maturation threshold #####
####################################################################

size_1SW_scn<-list()
size_1SW_scn_residents<-list()
for (scn in 1:length(EXPE)) {
  size_1SW<-list()
  size_1SW_residents<-list()
  for (pop in 1:npop) {
    size_1SW[[pop]] <- res_1SW_scn[[scn]][[1]][[pop]]$gFmid3[res_1SW_scn[[scn]][[1]][[pop]]$years > 45]
    size_1SW_residents[[pop]] <- res_1SW_scn[[scn]][[1]][[pop]]$gFmid3[res_1SW_scn[[scn]][[1]][[pop]]$years > 45 & res_1SW_scn[[scn]][[1]][[pop]]$CollecID == pop]
    for (simul in 2:nSIMUL) {
      size_1SW[[pop]] <- c(size_1SW[[pop]], res_1SW_scn[[scn]][[simul]][[pop]]$gFmid3[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45])
      size_1SW_residents[[pop]] <- c(size_1SW_residents[[pop]], res_1SW_scn[[scn]][[simul]][[pop]]$gFmid3[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45 & res_1SW_scn[[scn]][[simul]][[pop]]$CollecID == pop])
    }
  }
  size_1SW_scn[[scn]] <- size_1SW
  size_1SW_scn_residents[[scn]] <- size_1SW_residents
}

pop_group1<-c(1,3,5,10,12,14) #sink
pop_group2 <- c(6,8,15) #neutral
pop_group3 <- c(2,4,7,9,11,13) #source

group1<-array(,dim=c(length(pop_group1),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group1) {
    i=i+1
    group1[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group1[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}
group2<-array(,dim=c(length(pop_group2),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group2) {
    i=i+1
    group2[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group2[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}
group3<-array(,dim=c(length(pop_group3),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group3) {
    i=i+1
    group3[i,scn] <- median(size_1SW_scn_residents[[scn]][[pop]])
    #group3[i,scn] <- median(size_1SW_scn[[scn]][[pop]])
  }
}

plot(NULL, xlim=c(0.5,3), ylim=c(30,50), xlab="Population",ylab="Adult 1SW genotypic male sea maturation threshold", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
for (scn in 1:length(EXPE)) {
  points(rep(1+gap[scn],each=nrow(group1)), group1[,scn], col=viridis[scn],pch=20)#color_transparent(col)
  points(rep(1.7+gap[scn],each=nrow(group2)), group2[,scn], col=viridis[scn],pch=20)
  points(rep(2.4+gap[scn],each=nrow(group3)), group3[,scn], col=viridis[scn],pch=20)
}

colnames(group1) <- c(0,10,20,30,40,50)
sink_melted<-melt(group1)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

colnames(group2) <- c(0,10,20,30,40,50)
sink_melted<-melt(group2)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1.7+gap[1:6],each=3),lw_sink$fitted, col="black",lwd=5) #3 #5

colnames(group3) <- c(0,10,20,30,40,50)
sink_melted<-melt(group3)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(2.4+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

####################################################################
############ Fig.5 D) Proportion of 1SW in returns #################
####################################################################

ratio <- array(,dim=c(nYears+nInit,npop,nSIMUL, length(EXPE)))
Npops <- array(,dim=c(npop,nSIMUL, length(EXPE)))
for (scn in 1:length(EXPE)){
  for (simul in 1:nSIMUL){
    for (pop in 1:npop) {
      for (year in 1:(nYears+nInit)) {
        res_1SW_year <- subset(res_1SW_scn[[scn]][[simul]][[pop]], years == year)
        res_MSW_year <- subset(res_MSW_scn[[scn]][[simul]][[pop]], years == year)
        #only residents
        ratio[year,pop,simul,scn] <- nrow(res_1SW_year[res_1SW_year$CollecID==pop,]) / (nrow(res_1SW_year[res_1SW_year$CollecID==pop,]) + nrow(res_MSW_year[res_MSW_year$CollecID==pop,]))
        # with immigrants
        #ratio[year,pop,simul,scn] <- nrow(res_1SW_year) / (nrow(res_1SW_year) + nrow(res_MSW_year))
      }
      Npops[pop,simul,scn] <- median(ratio[45:50,pop,simul,scn], na.rm=T)
    }
  }
}

pop_group1<-c(1,3,5,10,12,14) #sink
pop_group2 <- c(6,8,15) #neutral
pop_group3 <- c(2,4,7,9,11,13) #source

group1<-array(,dim=c(length(pop_group1),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group1) {
    i=i+1
    group1[i,scn] <- median(Npops[pop,,scn], na.rm=T)
  }
}
group2<-array(,dim=c(length(pop_group2),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group2) {
    i=i+1
    group2[i,scn] <- median(Npops[pop,,scn], na.rm=T)
  }
}
group3<-array(,dim=c(length(pop_group3),length(EXPE)))
for (scn in 1:length(EXPE)) {
  i=0
  for (pop in pop_group3) {
    i=i+1
    group3[i,scn] <- median(Npops[pop,,scn], na.rm=T)
  }
}

plot(NULL, xlim=c(0.5,3), ylim=c(0.78,.88), xlab="Population",ylab="Proportion of 1SW in returns", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
for (scn in 1:length(EXPE)) {
  points(rep(1+gap[scn],each=nrow(group1)), group1[,scn], col=viridis[scn],pch=20)
  points(rep(1.7+gap[scn],each=nrow(group2)), group2[,scn], col=viridis[scn],pch=20)
  points(rep(2.4+gap[scn],each=nrow(group3)), group3[,scn], col=viridis[scn],pch=20)
}

colnames(group1) <- c(0,10,20,30,40,50)
sink_melted<-melt(group1)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

colnames(group2) <- c(0,10,20,30,40,50)
sink_melted<-melt(group2)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(1.7+gap[1:6],each=3),lw_sink$fitted, col="black",lwd=5) #3 #5

colnames(group3) <- c(0,10,20,30,40,50)
sink_melted<-melt(group3)
sink_melted$variable <- as.numeric(sink_melted$Var2)
lw_sink<-loess(value ~ variable, sink_melted,span=1)
lines(rep(2.4+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=5) #6 #5

####################################################################
####### Fig.6 Genetic variance of adult 1SW growth potential #######
####################################################################

genVar_scn<-array(,dim=c(nSIMUL,npop,length(EXPE)))
genVar_scn_residents<-array(,dim=c(nSIMUL,npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (pop in 1:npop) {
    for (simul in 1:nSIMUL) {
      individuals <- res_1SW_scn[[scn]][[simul]][[pop]]$gG[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45] #gFmid4
      individuals_residents <- res_1SW_scn[[scn]][[simul]][[pop]]$gG[res_1SW_scn[[scn]][[simul]][[pop]]$years > 45 & res_1SW_scn[[scn]][[simul]][[pop]]$CollecID == pop]
      
      genVar_scn[simul,pop,scn] <- var(individuals)
      genVar_scn_residents[simul,pop,scn] <- var(individuals_residents)
    }
  }
}
genVar_scn_residents_median <- array(,dim=c(npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (pop in 1:npop) {
    genVar_scn_residents_median[pop,scn] <- median(genVar_scn_residents[,pop,scn], na.rm = T)
  }
}

genVar_df<-as.data.frame(genVar_scn_residents_median)
genVar_df$Type<-datapop$Type

genVar_df[7,7]<-"source"

plot(NULL, xlim=c(0.5,3), ylim=c(0.003,.006), xlab="",ylab="Genetic variance of gG", xaxt='n')
mtext(c("Sink","Neutral","Source"), side = 1, line = 2, outer = FALSE, at = c(1,1.7,2.4))
mtext(c("0","10","20","30","40","50","0","10","20","30","40","50","0","10","20","30","40","50"), side = 1, line = .8, 
      outer = FALSE, at = c(1+gap,1.7+gap,2.4+gap),cex = .9)
#sink
sink<-genVar_df[genVar_df$Type=="sink",]
colnames(sink) <- c(0,10,20,30,40,50)
sink_melted<-melt(sink)
sink_melted$variable <- as.numeric(sink_melted$variable)
lw_sink<-loess(value ~ variable, sink_melted)
for (pop in which(genVar_df$Type=="sink")) {
  points(1+gap[1:6],genVar_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1+gap[1:6],each=6),lw_sink$fitted, col="black",lwd=4)

#neutral
neutral<-genVar_df[genVar_df$Type=="neutral",]
colnames(neutral) <- c(0,10,20,30,40,50)
neutral_melted<-melt(neutral)
neutral_melted$variable <- as.numeric(neutral_melted$variable)
lw_neutral<-loess(value ~ variable, neutral_melted)
for (pop in which(genVar_df$Type=="neutral")) {
  points(1.7+gap[1:6],genVar_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(1.7+gap[1:6],each=3),lw_neutral$fitted, col="black",lwd=4)

#source
source<-genVar_df[genVar_df$Type=="source",]
colnames(source) <- c(0,10,20,30,40,50)
source_melted<-melt(source)
source_melted$variable <- as.numeric(source_melted$variable)
lw_source<-loess(value ~ variable, source_melted)
for (pop in which(genVar_df$Type=="source")) {
  points(2.4+gap[1:6],genVar_df[pop,1:6], col=viridis[1:6], pch=20)
}
lines(rep(2.4+gap[1:6],each=6),lw_source$fitted, col="black",lwd=4)



