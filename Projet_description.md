# Projet MEPI  : Analyse de sensibilité d'un modèle épidémiologique


## Détermination du modèle 

On considère dans ce modèle :

**3 classes d'âge :** 
- *A* : Nouveaux-nés
- *B* : Jeunes
- *C* : Adultes
- *N* : Total

**4 états de santé :** 
- *1* : Sain
- *2* : Infecté neutre
- *3* : Infecté infectieux
- *4* : Remis

Sur une période de 2 ans soit 2x365 jours. 

En sortie on a une matrice avec en ligne les paramètres et 4 colonnes. 

*Boucle des scénarios de l'échantillonnage de l'AS*

15 paramètres au total avec





# Questions

1. 

a. Le modèle est un modèle déterministe (sans facteur aléatoire) à compartiments (états de santé et classe d'âge) à temps discret (par jour).

b. Les processus biologiques en jeu sont : 

* La mort naturelle
* La mort liée à la maladie
* La naissance de nouveaux-nés
* La croissance (changement de classe d'âge) linéaire unidirectionnelle (sans regression)
* transmission de la maladie 
* période de latence entre état infecté neutre et infecté infectieux (incubation)
* période d'immunité 


c. Equations associées : 
```

```
d. 
Schéma des transitions entre états

e. Les hypothèses principales du modèle sont : 

Paramètres démographiques : 

* La population se renouvelle (naissances et morts)
* La mortalité naturelle est la même par classe d'âge
* Les jeunes et les adultes se reproduisent dans des proportions différentes 
* Le passage d'une classe d'âge à une autre est identique pour tous les états de santé d'une classe d'âge 

Paramètres liés à l'AP : 

* L'AP peut uniquement infecter les individus sains, et de toutes les classes d'âge
* Une période de latence précède l'état infectieux de l'individu infecté
* Après la période d'immunité, l'individu redevient sensible à l'AP

f. 
Initialement il y a : 
* 27 nouveaux-nés sains
* 23 jeunes sains
* 36 adultes sains
* 1 adulte infecté infectieux 

g. Les sorties possibles du modèle sont : 

* Taux d'infectés (neutres et infectieux) le dernier jour 
* Incidence de l'AP au dernier 
* Journée avec le maximum d'infectés infectieux sur la période de temps étudiée
* Incidence de l'AP pour 1 an

h. 
Paramètres liés à la population étudiée :


- *K* : **100**, nombre maximal d'individus que le milieu peut supporter
- *sr* : **0,5**, sex-ratio
- *m1** : **0,0014**, mortalité naturelle des nouveaux-nés
- *m2* : **0,00029**, mortalité naturelle des jeunes
- *m3* : **0,0019**, mortalité naturelle des adultes
- *f2* : **0,0019**, taux de fécondité des jeunes
- *f3* : **0,0082**, taux de fécondité des adultes
- *portee* : **5**, effectif maximal d'une portée 
- *t1* : **1/365**, probabilité de passage de la classe d'âge "nouveau-né" à "jeune"
- *t2* : **1/365**, probabilité de passage de la classe d'âge "jeune" à "adulte"

**Paramètres liés à la maladie :**
- *trans* : **0,3**, force d'infection
- *lat* : **1/5**, taux de la latence
- *rec* : **1/20**, taux de passage à un état d'immmunité
- *loss* : **1/100**, taux de passage d'un état d'immunité à un état sensible (sain)
- *madd* : **0,001**, mortalité lié à l'infection par l'AP 

i. 


2. 

