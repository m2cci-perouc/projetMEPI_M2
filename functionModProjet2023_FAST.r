### Modele dont la sensibilite doit etre analysee dans le cadre du projet MODE-MPI 2023-2024

### Le modele est ici defini sous forme de fonction pour faciliter vos analyses de sensibilite (AS)
### La fonction renvoie les sorties ponctuelles qui sont a analyser dans l'AS

library(ggplot2)





# Fonction du modele ------------------------------------------------------


modAppli <- function(parametre){  

  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (ou simulations), en jours
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de parametres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)
  
  
  
  # boucle des scenarios de l'echantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 

    
    # STRUCTURE & PARAMETRES DU MODELE

    # Parametres demographiques
    K = parametre[i,1];		# nombre maximal d'individus que le milieu peut supporter
    sr = parametre[i,2];	# sex-ratio
    m1 = parametre[i,3];	# mortalité naturelle des nouveaux-nés
    m2 = parametre[i,4];	# mortalité naturelle des jeunes
    m3 = parametre[i,5];	# mortalité naturelle des adultes
    f2 = parametre[i,6];	# taux de fécondité des jeunes
    f3 = parametre[i,7];	# taux de fécondité des adultes
    portee = parametre[i,8];	# effectif maximal d'une portée
    t1 = parametre[i,9];	# probabilité de passage de la classe d'âge "nouveau-né" à "jeune"
    t2 = parametre[i,10];	# probabilité de passage de la classe d'âge "jeune" à "adulte"

    # Parametres lies a l'AP
    trans = parametre[i,11]; # force d'infection
    lat = parametre[i,12];	# taux de la latence
    rec = parametre[i,13];	# taux de passage à un état d'immmunité
    loss = parametre[i,14];	# taux de passage d'un état d'immunité à un état sensible (sain)
    madd = parametre[i,15];	# mortalité lié à l'infection par l'AP

    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'age en ligne (derniere ligne = pop tot), 
    #etat de sante en colonne, pas de temps (dimension 3)
    eff<- matrix(0,16,temps) ;#objet de stockage des effectifs 
    nvinf <- array(0, dim=c(temps));
    # conditions initiales (la population est a sa structure d'equilibre, calculee par ailleurs)
    #MAT[classe_age, etat_sante, jour]
    MAT[1,1,1] <- 27; # Effectif des individus nouveaux-nes sains le premier jour
    MAT[2,1,1] <- 23; # Effectif des individus jeunes sains le premier jour
    MAT[3,1,1] <- 36; # Effectif des individus adultes sains le premier jour 
      #Tous les individus sains
    
    MAT[3,3,1] <- 1;  # Effectif des individus adultes infectes infectieux le premier jour

    # effectifs par etat de sante
    MAT[4,1,1] <- sum(MAT[1:3,1,1]);#Effectif total des individus de toutes les classes d'age sains le premier jour
    MAT[4,2,1] <- sum(MAT[1:3,2,1]);#Pour le deuxieme etat de sante 
    MAT[4,3,1] <- sum(MAT[1:3,3,1]);#Pour le troisieme etat de sante 
    MAT[4,4,1] <- sum(MAT[1:3,4,1]);#Pour le quatrieme etat de sante 
    
    #Objet eff a remplir pour avoir les etats initiaux
    eff[1,1]<-MAT[1,1,1]
    eff[5,1]<-MAT[2,1,1]
    eff[9,1]<-MAT[3,1,1]
    eff[11,1]<-MAT[3,3,1]
    
    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
     # Les nouveaux-nes
      # RQ : les naissances sont XX, les nouveaux nes etant dans l'etat XX
      N <- sum(MAT[4,,t]);	# taille de la pop en t
	MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t] 
	                + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
	MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
	MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
	MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
     # Les jeunes
	MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
	MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
	MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
	MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
     # Les adultes
	MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
	MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
	MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
	MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
     # calcul des effectifs  N par etat de sante
	MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); 
	MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); 
	MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); 
	MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
	nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N

# 	#stockage des donnees effectif dans un objet eff -----------------------

	
	  #Nouveaux-nes
	eff[1,t+1]<-MAT[1,1,t+1];
	eff[2,t+1]<-MAT[1,2,t+1];
	eff[3,t+1]<-MAT[1,3,t+1];
	eff[4,t+1]<-MAT[1,4,t+1];
	  #Jeunes
	eff[5,t+1]<-MAT[2,1,t+1];
	eff[6,t+1]<-MAT[2,2,t+1];
	eff[7,t+1]<-MAT[2,3,t+1];
	eff[8,t+1]<-MAT[2,4,t+1];
	  #Adultes
	eff[9,t+1]<-MAT[3,1,t+1];
	eff[10,t+1]<-MAT[3,2,t+1];
	eff[11,t+1]<-MAT[3,3,t+1];
	eff[12,t+1]<-MAT[3,4,t+1];
	
	  #Total
	eff[13,t+1]<-MAT[4,1,t+1];
	eff[14,t+1]<-MAT[4,2,t+1];
	eff[15,t+1]<-MAT[4,3,t+1];
	eff[16,t+1]<-MAT[4,4,t+1];
	
	
	
    }# fin boucle temps

    # sorties ponctuelles a analyser
    # Taux des infectes au dernier jour
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # Incidence le dernier jour
    sortie2 <- nvinf[temps]
    # Maximum d'infectes infectieux au cours du temps
    sortie3 <- max(MAT[4,3,1:temps])
    # Incidence (total de nouveaux infectes par an)
    sortie4 <- sum(nvinf[1:365])
    
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
    
  }# fin boucle scenarios AS
  #return(sorties)
  return(sorties)
} # fin fonction du modele

# END



# Simulations du modele ---------------------------------------------------


#Valeur des parametres : 
ValNominale <- c(100,0.5,0.0014,0.00029,0.0019,0.0019,0.0082,5,1/365,1/365, #10 parametres
                 0.3,1/5,1/20,1/100,0.001)# 5 parametres maladie
PAR<-matrix(ValNominale, nrow=1)
PAR


# Tests -------------------------------------------------------------------

# Simulation 1
sim1<-modAppli(parametre = PAR)
# length(sim1) # 11680
sim1<-t(sim1)
colnames(sim1)<- 1:16 

# reshape dataframe to have only one value column
library(reshape2)
sim1<-melt(sim1) # Mettre toutes les données dans 1 seule colonne + colonne identifiant (1:16)
# create a x vector from 1 to the number of value for each variable (here 730)
# sim1$x<-rep(1:length(sim1$value[sim1$variable==1]))
npoints<- 730*16
npoints
sim1$x<-seq_len(npoints)

# Representations graphiques ----------------------------------------------
#https://stackoverflow.com/questions/43850567/how-can-one-plot-the-rows-of-a-two-dimensional-array-in-one-plot-using-ggplot2

# Nom des axes et des differents etats
library(tidyverse)

sim1<-sim1 %>% rename( Etat = Var2,
                 "Temps (j)" = Var1,
                 Effectif = value)
# sim1$Etat<-as.factor(sim1$Etat)
# sim1 %>% mutate(Etat = replace(Etat, Etat == 1, c("Snn")))
# sim1$Etat<-ifelse(sim1$Etat==1,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==2,"INnn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==3,"IInn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==4,"nn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==5,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==6,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==7,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==8,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==9,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==10,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==11,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==12,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==13,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==14,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==15,"Snn",sim1$Etat)
# sim1$Etat<-ifelse(sim1$Etat==16,"Snn",sim1$Etat)


# Faire le graphique avec ggplot2
library(RColorBrewer)
library(ggnewscale)

# ggplot(data=sim1, aes(x="Temps (j)", y=Effectif, color=Etat, linetype=Etat)) +
#   geom_line(data = sim1, mapping = aes(linetype = Etat)) +
#   scale_linetype_manual(values =c("Snn"= 2, "INnn"=3))+
#   new_scale("linetype")+
#   geom_line(
#     data= sim1 %>% mutate() 
#   )
#   scale_color_manual(values= c(brewer.pal(10, "Set3"), brewer.pal(2, "Set3")))+
#   ggtitle("Evolution des effectifs par classe d'age et état de santé")+
#   theme_bw()



# Analyse de sensibilité du modèle ----------------------------------------


# Methode OAT -------------------------------------------------------------



# Méthode Morris ----------------------------------------------------------


  

# Méthode FAST ------------------------------------------------------------

#  Lien pour la méthode sur R : https://rdrr.io/cran/fast/man/sensitivity.html 

# sensitivity(x, numberf, order = 4, make.plot = FALSE, show.legend
#                  = TRUE, plot.max = max(ff[-1]), include.total.variance
#                  = FALSE, cukier = TRUE, names = paste(sep = "", "P",
#                  1:numberf), main = "", xlab = "frequency", ylab =
#                  "Fourier Coef", pch = rep(0, numberf), col =
#                  (1:numberf) + 1, reorder = 1:numberf, ...)

# Exemple 
  
# Utilisation de la fonction fast99 du package sensitivity 
library(sensitivity)
  
# 2 manieres de lancer l'analyse FAST :
  # Soit en mettant directement le modele dans la fonction
  # Soit en appelant la fonction avec le modele nul, puis en utilisant la fonction tell()
# Utile dans le cas d'un mmodele assez lourd et gourmand en analyse

# Définir les paramètres et les plages de valeurs
parameters <- c("K","sr", "m1", "m2","m3","f2","f3","portee","t1","t2",
                "trans","lat","rec","loss","madd")
fact<-0.25  # On cree des bornes a +/- 25% de la valeur entree

# Version "lourde"
param_ranges<- list(
  list(min = PAR[1]-PAR[1]*fact, max = PAR[1]+PAR[1]*fact), # pour K
  list(min = PAR[2]-PAR[2]*fact, max = PAR[2]+PAR[2]*fact), # pour sr
  list(min = PAR[3]-PAR[3]*fact, max = PAR[3]+PAR[3]*fact), # pour m1
  list(min = PAR[4]-PAR[4]*fact, max = PAR[4]+PAR[4]*fact), # pour m2
  list(min = PAR[5]-PAR[5]*fact, max = PAR[5]+PAR[5]*fact), # pour m3
  list(min = PAR[6]-PAR[6]*fact, max = PAR[6]+PAR[6]*fact), # pour f2
  list(min = PAR[7]-PAR[7]*fact, max = PAR[7]+PAR[7]*fact), # pour f3
  list(min = PAR[8]-PAR[8]*fact, max = PAR[8]+PAR[8]*fact), # pour portee
  list(min = PAR[9]-PAR[9]*fact, max = PAR[9]+PAR[9]*fact), # pour t1
  list(min = PAR[10]-PAR[10]*fact, max = PAR[10]+PAR[10]*fact), # pour t2
  list(min = PAR[11]-PAR[11]*fact, max = PAR[11]+PAR[11]*fact), # pour trans
  list(min = PAR[12]-PAR[12]*fact, max = PAR[12]+PAR[12]*fact), # pour lat
  list(min = PAR[13]-PAR[13]*fact, max = PAR[13]+PAR[13]*fact), # pour rec
  list(min = PAR[14]-PAR[14]*fact, max = PAR[14]+PAR[14]*fact), # pour loss
  list(min = PAR[15]-PAR[15]*fact, max = PAR[15]+PAR[15]*fact) # pour madd
)# Creation de la liste avec les bornes par parametre

# Version boucle
# param_ranges<-list(min = PAR[1]-PAR[1]*fact, max = PAR[1]+PAR[1]*fact) # On cree l'objet contenant les bornes min et max de chaque parametre
# for(i in 2:length(PAR)-1){
#   param_ranges[i]<-list(param_ranges[i-1],list(min = PAR[i]-PAR[i]*fact, max = PAR[i]-PAR[i]*fact))
# }  ##### NE FONCTIONNE PAS, a perfectionner parce que la version lourde est vraiment pas opti 

View(param_ranges)

# Utiliser fast99 pour créer l'objet sa (avec model = NULL)
parameters
param_ranges
sa <- fast99(model = NULL, factors = parameters, n=15, q = "qunif", q.arg = param_ranges)
sa$X

# Utiliser tell pour générer le design expérimental
# experiments <- tell(sa, n = 10, distribution = "latin.hypercube", # pas sure que ce soit necessaire
#                     factor_levels = 5, params = parameters)

# Analyser la sensibilité avec FAST
# fast_results <- fast99(sa, x = experiments$X1, y = experiments$y)

# at this stage, only the design of experiment (sa$x) was generated
# the response is computed "manually":
sa$x
n <- nrow(sa$x)
y <- numeric(n)
tell(sa, sim1) # tell(x= sensitivity analysis object, y = the response)
print(sa)
plot(x)
# Afficher les résultats
print(fast_results)


# Analyse Sobol
y <- sobol.fun(sa$x) # at this place could be a
# call to an external code

# then, the sensitivity analysis:

tell(sa, y)
print(sa)