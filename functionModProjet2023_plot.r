### Modele dont la sensibilite doit etre analysee dans le cadre du projet MODE-MPI 2023-2024

### Le modele est ici defini sous forme de fonction pour faciliter vos analyses de sensibilite (AS)
### La fonction renvoie les sorties ponctuelles qui sont a analyser dans l'AS

modAppli <- function(parametre){  

  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de parametres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)
  #initialisation pour la sauvegarde des dynamique de la population
  sortiepop=matrix(0,16,2*365)

  # boucle des scenarios de l'echantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 

    # STRUCTURE & PARAMETRES DU MODELE

    # XX
    K = parametre[i,1];		# xx
    sr = parametre[i,2];	# xx
    m1 = parametre[i,3];	# mort naturelle classe 1
    m2 = parametre[i,4];	# mort naturelle classe 2
    m3 = parametre[i,5];	# mort naturelle classe 3
    f2 = parametre[i,6];	# xx
    f3 = parametre[i,7];	# xx
    portee = parametre[i,8];	# xx
    t1 = parametre[i,9];	# xx
    t2 = parametre[i,10];	# xx

    # XX
    trans = parametre[i,11]; # tx transmission
    lat = parametre[i,12];	# latence
    rec = parametre[i,13];	# xx
    loss = parametre[i,14];	# xx
    madd = parametre[i,15];	# xx

    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'age en ligne (derniere ligne = pop tot), etat de sante en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    
    #Plot
    plot(x = NA, y = NA, xlim = c(0, 730), ylim = c(0, 40),
         type = "n", xlab = "Temps (jours)", ylab = "Nombre d'individus")

    # conditions initiales (la population est a sa structure d'equilibre, calculee par ailleurs)
    MAT[1,1,1] <- 27; # xx
    MAT[2,1,1] <- 23; # xx
    MAT[3,1,1] <- 36; # xx
    MAT[3,3,1] <- 1;  # xx
    # effectifs par etat de sante
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); MAT[4,2,1] <- sum(MAT[1:3,2,1]); MAT[4,3,1] <- sum(MAT[1:3,3,1]); MAT[4,4,1] <- sum(MAT[1:3,4,1]);

    #condition initiale
    sortiepop[1,1]=27
    sortiepop[5,1]=23
    sortiepop[9,1]=36
    sortiepop[11,1]=1
    sortiepop[13,1]=sortiepop[1,1]+sortiepop[5,1]+sortiepop[9,1]
    sortiepop[14,1]=sortiepop[2,1]+sortiepop[6,1]+sortiepop[10,1]
    sortiepop[15,1]=sortiepop[3,1]+sortiepop[7,1]+sortiepop[11,1]
    sortiepop[16,1]=sortiepop[4,1]+sortiepop[8,1]+sortiepop[12,1]
    
    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
     # classe d'age 1 (nouveau-ne)
      # RQ : les naissances sont XX, les nouveaux nes etant dans l'etat XX
      N <- sum(MAT[4,,t]);	# taille de la pop en t
	MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
	MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
	MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
	MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
     # classe d'age 2 (jeune)
	MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
	MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
	MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
	MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
     # classe d'age 3 (adulte)
	MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
	MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
	MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
	MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
     # calcul des effectifs par etat de sante
	MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); 
	MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); 
	MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); 
	MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
	nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N

	#NN
	sortiepop[1,t+1]=MAT[1,1,t+1]
	sortiepop[2,t+1]=MAT[1,2,t+1]
	sortiepop[3,t+1]=MAT[1,3,t+1]
	sortiepop[4,t+1]=MAT[1,4,t+1]
  #J
	sortiepop[5,t+1]=MAT[2,1,t+1]
	sortiepop[6,t+1]=MAT[2,2,t+1]
	sortiepop[7,t+1]=MAT[2,3,t+1]
	sortiepop[8,t+1]=MAT[2,4,t+1]
  #A
	sortiepop[9,t+1]=MAT[3,1,t+1]
	sortiepop[10,t+1]=MAT[3,2,t+1]
	sortiepop[11,t+1]=MAT[3,3,t+1]
	sortiepop[12,t+1]=MAT[3,4,t+1]
  #pop totale
	sortiepop[13,t+1]=MAT[4,1,t+1]#sain
	sortiepop[14,t+1]=MAT[4,2,t+1]#IN
	sortiepop[15,t+1]=MAT[4,3,t+1]#II
	sortiepop[16,t+1]=MAT[4,4,t+1]#R
	
    }# fin boucle temps

    # sorties ponctuelles a analyser
    # taux des infectes au dernier jours
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # incidence du dernier jour
    sortie2 <- nvinf[temps]
    # nbr max d'infectes infectieux sur les 2 ans
    sortie3 <- max(MAT[4,3,1:temps])
    # incidence sur l'année
    sortie4 <- sum(nvinf[1:365])
    
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
  }# fin boucle scenarios AS
  lines(sortiepop[1,],col="green")
  lines(sortiepop[2,],col="purple")
  lines(sortiepop[3,],col="red")
  lines(sortiepop[4,],col="blue")
  lines(sortiepop[5,],lty = 2,col="green")
  lines(sortiepop[6,],lty = 2,col="purple")
  lines(sortiepop[7,],lty = 2,col="red")
  lines(sortiepop[8,],lty = 2,col="blue")
  lines(sortiepop[9,],lty = 3,col="green")
  lines(sortiepop[10,],lty = 3,col="purple")
  lines(sortiepop[11,],lty = 3,col="red")
  lines(sortiepop[12,],lty = 3,col="blue")
  #lines(sortiepop[13,],col="green")
  #lines(sortiepop[14,],col="purple")
  #lines(sortiepop[15,],col="red")
  #lines(sortiepop[16,],col="blue")
  
  return(sorties)
} # fin fonction du modele

# END
ValNominale = c(100, 0.5, 0.0014, 0.00029, 0.0019, 0.0019, 0.0082, 5, 1/365, 
                1/365, 0.3, 1/5, 1/20, 1/100, 0.001)

Par=matrix(ValNominale,nrow=1,ncol=15)
modAppli(Par)

ggplot(modAppli(Par))
