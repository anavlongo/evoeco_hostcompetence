#####Code to reproduce analyses from Longo, Lips, and Zamudio Philosophical Transactions Royal Society B
#Title: Evolutionary ecology of host competence after a chytrid outbreak in a naïve amphibian community
#Please contact corresponding author, Ana V. Longo, for questions regarding the code and data (ana.longo@ufl.edu)

library(ggplot2)
library(dplyr)
library(sjPlot)
library(MASS)
library(robustbase)
library(gvlma)
library(ape)
library(geiger)
library(picante)
library(vegan)
library(phytools)
library(coin)
library(pgirmess)
library(indicspecies)
library(meantables)
library(mosaic)

##Summarized information for Table S1
Individualdata <- read.table("Longo_et_al_PTRSB_individualdata.txt", header=T)
names(Individualdata)

##Getting Bd Prevalence and mean infection intensity per species from individual data
Individualdata %>% group_by(Stage2,Family, Species_Name) %>%
    summarise(Bdpos=sum(Bd_presabs), n = n(), avgzoosp = mean(Load), sdzp = sd(Load), maxzp = max(Load)) %>%
    mutate(Prevalence = Bdpos/n) %>% tab_df(title="Table S1")

##Comparing prevalence between extirpated species versus species that persisted after the epizootic
Individualdata %>%  filter(Stage2=="1_Epizootic") %>%
  group_by(Persistence_afterEpizootic) %>%
  summarise(Bdpos=sum(Bd_presabs), n = n()) %>%
  mutate(Prevalence = Bdpos/n) 
mymatrix <- matrix(c(282, 109, 245, 262),nrow=2,byrow=TRUE)
chisq.test(mymatrix,correct=F)

###Linear Models with summarized data
summarizeddata <- read.table("Longo_et_al_PTRSB_summarizeddata.txt",h=TRUE)
names(summarizeddata)
attach(summarizeddata)

lm1<- lm(prev~Survival)
gvmodel1 <- gvlma(lm1) 
summary(gvmodel1)

lm2<- lm(log10(zoosp+1)~Survival)
gvmodel2 <- gvlma(lm2) 
summary(gvmodel2)

# Figure 1A
plot(prev ~ Survival, ylim=c(0, 100), xlim=c(0, 1), las=1, type = "n", xlab = "Trait-based Persistence probability", ylab="Bd Prevalence (%)")
points(Survival[develop=="AL"], prev[develop=="AL"], cex=1.5)
points(Survival[develop=="DD"], prev[develop=="DD"], cex=1.5, pch=16)
abline(lm1, lwd=2)

# Figure 1B
plot(Survival, ylim=c(-1.5,5), xlim=c(0,1), las=1, type = "n", ylab = "Infection Intensity (log10)", xlab="Trait-based Persistence Probability")
points(Survival[develop=="AL"], log10(zoosp[develop=="AL"]), cex=1.5)
points(Survival[develop=="DD"], log10(zoosp[develop=="DD"]), cex=1.5, pch=16)
abline(lm2, lwd=2)

######Comparative methods following phytools
amphibiantree<-read.nexus("Longo_et_al_PTRSB_TreeCope.nex")
dichotomousphylogeny <- multi2di(amphibiantree, random = TRUE)

#tree with branches equal to 1
dichotomousphylogeny <- compute.brlen(dichotomousphylogeny, 1)
name.check(dichotomousphylogeny, summarizeddata)
is.binary(dichotomousphylogeny)
plot(dichotomousphylogeny)

##setting up the variables for the phylogenetic analyses 
prev <- summarizeddata$prev
zoosp <- summarizeddata$zoosp
Survival <- summarizeddata$Survival

names(prev) <- row.names(summarizeddata)
names(zoosp) <- row.names(summarizeddata)
names(Survival) <- row.names(summarizeddata)

dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Craugastor_azueroensis") 
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Craugastor_megacephalus")
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Craugastor_tabasarae")
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Ecnomiohyla_miliaria" )
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Pristimantis_museosus")
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Smilisca_sila" )
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Craugastor_punctariolus" )
dichotomousphylogeny <- drop.tip(dichotomousphylogeny, "Craugastor_bransfordii" )

plot(dichotomousphylogeny)

#calculate contrasts with expected variance
Contrastprev <- pic(prev, dichotomousphylogeny)
Contrastload <- pic(log(zoosp+1), dichotomousphylogeny)
ContrastSurvival <- pic(Survival, dichotomousphylogeny)

###Regressions accounting for evolutionary history
RegressPICS1 <- lm(Contrastprev ~ ContrastSurvival-1)
summary.lm(RegressPICS1)

# Figure 1C
plot(Contrastprev~ContrastSurvival, las=1, ylab = "Contrasts Bd Prevalence", xlab="Contrasts Trait-based Persistence Probability", pch=18, cex=1.5)
abline(RegressPICS1, lwd=2, col="black")

RegressPICS2 <- lm(Contrastload ~ ContrastSurvival -1)
summary.lm(RegressPICS2)

# Figure 1D
plot(Contrastload ~ ContrastSurvival, ylab = "Contrasts Infection intensity (log10)", xlab="Contrasts Trait-based Persistence Probability", pch=18, cex=1.5)
abline(RegressPICS2, lwd=2, col="black")

###Phylogenetic signal tests following phytools
amphibiantree<-read.nexus("Longo_et_al_PTRSB_TreeCope.nex")
dichotomousphylogeny <- multi2di(amphibiantree, random = TRUE)
dichotomousphylogeny <- compute.brlen(dichotomousphylogeny, 1)
logload <- log10(Individualdata$Load+1)
names(logload) <- row.names(Individualdata$Species_Name)
#Log10(inf + 1) 
xbar1 <- as.matrix(tapply(logload, Individualdata$Species_Name, mean))
xvar1 <- as.matrix(tapply(logload, Individualdata$Species_Name, var))
n1 <- as.matrix(tapply(logload, Individualdata$Species_Name, length))
names(xbar1) <- row.names(xbar1)
names(xvar1) <- row.names(xvar1)
names(n1) <- row.names(n1)
# replace NA with mean (m) or pooled (p) variance
xvarm1<-xvarp1<-xvar1
xvarm1[is.na(xvar1)]<-mean(xvar1,na.rm=TRUE)
xvarp1[is.na(xvar1)]<-0
xvarp1[is.na(xvar1)]<-  sum((n1-1)*xvarp1/(sum(n1[n1>1])-length(n1[n1>1])))
#BlombergsK for infection parameters
# compute K ignoring sampling error
phylosig(dichotomousphylogeny,xbar1, test=TRUE, nsim=1000)
# compute K with sampling error, pooled variance for n=1 for infection intensity
phylosig(dichotomousphylogeny,xbar1,se=sqrt(xvarp1/n1), test=TRUE, nsim=1000)

#separating trees
DDtree <- drop.tip(dichotomousphylogeny, c("Bolitoglossa_schizodactyla", 
                                           "Bolitoglossa_colonnea", "Lithobates_warszewitschii", "Nelsonophryne_aterrima",
                                           "Leptodactylus_pentadactylus", "Atelopus_zeteki", "Cochranella_euknemos", "Ecnomiohyla_miliaria", "Espadarana_prosoblepon",         
                                           "Hyalinobatrachium_colymbiphyllum", "Hylomantis_lemur",  "Hyloscirtus_colymba",       
                                           "Hyloscirtus_palmeri", "Incilius_coniferus", "Rhaebo_haematiticus", "Rhinella_marina",                 
                                           "Sachatamia_albomaculata", "Sachatamia_ilex", "Smilisca_phaeota", "Smilisca_sila",
                                           "Gastrotheca_cornuta", "Dendrobates_auratus", "Silverstoneia_nubicola", "Silverstoneia_flotator",
                                           "Colostethus_panamansis", "Colostethus_pratti"))
ALtree <- drop.tip(dichotomousphylogeny, c("Bolitoglossa_schizodactyla", 
                                           "Bolitoglossa_colonnea", "Lithobates_warszewitschii", "Nelsonophryne_aterrima",
                                           "Leptodactylus_pentadactylus", "Gastrotheca_cornuta", "Dendrobates_auratus", 
                                           "Silverstoneia_nubicola", "Silverstoneia_flotator", "Colostethus_panamansis", 
                                           "Colostethus_pratti", "Craugastor_azueroensis", "Craugastor_bransfordii",        
                                           "Craugastor_crassidigitus", "Craugastor_gollmeri", "Craugastor_megacephalus", "Craugastor_noblei" ,              
                                           "Craugastor_punctariolus", "Craugastor_tabasarae", "Craugastor_talamancae",                 
                                           "Diasporus_diastema", "Diasporus_vocator", "Pristimantis_caryophyllaceus",    
                                           "Pristimantis_cerasinus", "Pristimantis_cruentus","Pristimantis_gaigeae", "Pristimantis_museosus",          
                                           "Pristimantis_pardalis" , "Pristimantis_ridens", "Strabomantis_bufoniformis"))
phylosig(DDtree,xbar1, test=TRUE)
phylosig(DDtree,xbar1,se=sqrt(xvarp1/n1), test=TRUE)
phylosig(ALtree,xbar1, test=TRUE)
phylosig(ALtree,xbar1,se=sqrt(xvarp1/n1), test=TRUE)

#K for prevalence
amphibiantree<-read.nexus("Longo_et_al_PTRSB_TreeCope.nex")
dichotomousphylogeny <- multi2di(amphibiantree, random = TRUE)
dichotomousphylogeny <- compute.brlen(dichotomousphylogeny, 1)
prev <- summarizeddata$prev
names(prev) <- row.names(summarizeddata)
phylosig(dichotomousphylogeny,prev, test=TRUE)
phylosig(DDtree,prev, test=TRUE)
phylosig(ALtree,prev, test=TRUE)

#K for Persistence
surv <- summarizeddata$Survival
names(surv) <- row.names(summarizeddata)
phylosig(dichotomousphylogeny,surv, test=TRUE)
phylosig(DDtree,surv, test=TRUE)
phylosig(ALtree,surv, test=TRUE)

###Kruskal wallis Test by families
Individualdata %>%
  group_by(Family)%>% tally() %>% filter(n>10) 

KWT<-Individualdata %>%
  filter(Family %in% c("Bufonidae", "Centrolenidae","Craugastoridae","Dendrobatidae", "Eleutherodactylidae", "Hylidae", "Plethodontidae", "Strabomantidae")) %>% 
  as.data.frame()

KWT$Family <- factor(KWT$Family)
kruskal_test(KWT$Load ~ KWT$Family, distribution = approximate(nresample = 10000))
kruskalmc(KWT$Load ~ KWT$Family)

library(indicspecies)
Cope <- read.table("Longo_et_al_PTRSB_CopeAbundanceprop.txt", header=T)
Copeenv <-read.table("Longo_et_al_PTRSB_CopeEnv.txt", header=T)
copept <- multipatt(Cope, func="r.g", Copeenv$DiseaseStatus, control = how(nperm=999))
summary(copept, indvalcomp=TRUE)

##Figure 2
Individualdata %>%  filter(Stage2=="1_Epizootic") %>%
  group_by(Development) %>%
  summarise(Bdpos=sum(Bd_presabs), n = n()) %>%
  mutate(Prevalence = Bdpos/n) 
mymatrix <- matrix(c(282, 109, 245, 262),nrow=2,byrow=TRUE)
chisq.test(mymatrix,correct=F)

ALvsDD <- Individualdata %>% filter(Stage2=="1_Epizootic") %>%
  group_by(Development) %>%
  summarise( 
    n=n(),
    Bdpos=sum(Bd_presabs)) %>%
    mutate(Prevalence = Bdpos/n,
      lower = lapply(n, prop.test, n = sum(n)), 
      upper = sapply(lower, function(x) x$conf.int[2]), 
      lower = sapply(lower, function(x) x$conf.int[1]))

ALvsDD2 <- Individualdata %>% filter(Stage2=="1_Epizootic") %>%
  group_by(Development) %>%
  mean_table(Bd_presabs)

##Figure 2A
ggplot(ALvsDD2) +
  geom_bar( aes(x=group_cat, y=mean, fill=group_cat), stat="identity", color="black") +
  geom_errorbar( aes(x=group_cat, ymin=lcl, ymax=ucl), width=0.4, colour="black", alpha=0.9, size=1.5) + 
  scale_fill_manual(values=c("white", "gray47"))+ scale_y_continuous(breaks=seq(0,1, 0.2), limits=c(0, 1))+
   ylab("Probability of Infection")+theme_classic()+ theme(axis.text=element_text(size=12), axis.title.y = element_text(face="bold", size=16))

##Barplot Figure 2B
Fig2b <- Individualdata %>% filter(Stage2=="1_Epizootic") %>%
      filter(Family %in% c("Bufonidae", "Centrolenidae","Craugastoridae","Dendrobatidae", "Eleutherodactylidae", "Hylidae", "Plethodontidae", "Strabomantidae")) 
Figure2 <- ggplot(Fig2b, aes(as.factor(Family), fill = factor(Level))) +
  geom_bar(colour="black", position = position_fill(reverse = TRUE)) + scale_fill_manual(values=c("#FFFFFF", "#b7e3e7", "#5cbfc9","#0fa1af", "#4e857b", "#a65d33","#e14303"))+
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size = 0.5)
  )+coord_flip()

positions <- c("Centrolenidae", "Dendrobatidae", "Bufonidae", "Hylidae", "Strabomantidae", "Craugastoridae", "Eleutherodactylidae", "Plethodontidae")
Figure2 + scale_x_discrete(limits = positions)

#get sample sizes for each bar
Fig2b %>%
  group_by(Family)%>% tally()

###plotting asymmetries###
Asymmetries_file <- read.table("Longo_et_al_PTRSB_TableS5_asymmetries.txt", header=T)
head(Asymmetries_file)

#filtered final table for observations with N larger than 10 individuals and only occuring on year 2004
Final_asymmetries <- Asymmetries_file[which(Asymmetries_file$N > 10 & Asymmetries_file$Year == 2004),]

#Figure 3
Figure4 <- ggplot(Final_asymmetries, aes(Infection, pi)) 
Figure4 + geom_point(aes(colour=Shedding, size=Abundance)) + scale_size(range = c(5, 10))+
  scale_colour_gradient(low="#0FA1AF", high="#E14303") + theme_bw(base_size = 20)+theme(strip.background = element_rect(fill="white"))

#Figure 4 - analyses for community versus 3 focal species
Figure4A <- ggplot(Individualdata, aes(as.factor(Time), fill = factor(Level))) +
  geom_bar(colour="black", position = position_fill(reverse = TRUE)) + scale_fill_manual(values=c("#FFFFFF", "#b7e3e7", "#5cbfc9","#0fa1af", "#4e857b", "#a65d33","#e14303"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size = 0.5)
  )
Figure4A

Fig4_Atelopus <- Individualdata %>% filter(Species_Name=="Atelopus_zeteki")
Figure4B <- ggplot(Fig4_Atelopus, aes(as.factor(Time), fill = factor(Level))) +
  geom_bar(colour="black", position = position_fill(reverse = TRUE)) + scale_fill_manual(values=c("#FFFFFF", "#b7e3e7", "#5cbfc9","#0fa1af", "#4e857b", "#a65d33","#e14303"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size = 0.5)
  )
Figure4B

Fig4_Espadarana <- Individualdata %>% filter(Species_Name=="Espadarana_prosoblepon")
Figure4C <- ggplot(Fig4_Espadarana, aes(as.factor(Time), fill = factor(Level))) +
  geom_bar(colour="black", position = position_fill(reverse = TRUE)) + scale_fill_manual(values=c("#FFFFFF", "#b7e3e7", "#5cbfc9","#0fa1af", "#4e857b", "#a65d33","#e14303"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size = 0.5)
  )
Figure4C

Fig4_Pristi <- Individualdata %>% filter(Species_Name=="Pristimantis_cruentus")
Figure4D <- ggplot(Fig4_Pristi, aes(as.factor(Time), fill = factor(Level))) +
  geom_bar(colour="black", position = position_fill(reverse = TRUE)) + scale_fill_manual(values=c("#FFFFFF", "#b7e3e7", "#5cbfc9","#0fa1af", "#4e857b", "#a65d33","#e14303"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size = 0.5)
  )
Figure4D

##Density
Densitydata <- read.table("Longo_et_al_PTRSB_density.txt", header=T)
names(Densitydata)

#Figure4E
boxplot(Densitydata$Density~Densitydata$mo, las=2, col="white")
means1 <- tapply(Densitydata$Density, Densitydata$mo, mean)
points(means1, pch=20, cex=1.5)
#Figure4F
boxplot(Densitydata$Atelopus_zeteki~Densitydata$mo, las=2, col="white")
means2 <- tapply(Densitydata$Atelopus_zeteki, Densitydata$mo, mean)
points(means2, pch=20, cex=1.5)
#Figure4G
boxplot(Densitydata$Espadarana_prosoblepon~Densitydata$mo, las=2, col="white")
means3 <- tapply(Densitydata$Espadarana_prosoblepon, Densitydata$mo, mean)
points(means3, pch=20, cex=1.5)
#Figure4H
boxplot(Densitydata$Pristimantis_cruentus~Densitydata$mo, las=2, col="white")
means4 <- tapply(Densitydata$Pristimantis_cruentus, Densitydata$mo, mean)
points(means4, pch=20, cex=1.5)

#Fin! 
