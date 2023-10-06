# Projet MEPI  : Analyse de sensibilité d'un modèle épidémiologique


## Détermination du modèle 

On considère dans ce modèle :

**3 classes d'âge :** 
- *NN* : Nouveaux-nés
- *J* : Jeunes
- *A* : Adultes
- *N* : Population totale

**4 états de santé :** 
- *S* : Sain
- *IN* : Infecté neutre
- *II* : Infecté infectieux
- *R* : Remis

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

**Classe d'âge NN (nouveaux nés)**

$$S_{NN_{t+1}} = S_{NN_t} (1 - m_{1} - t_{1} - trans . N_{II_t}) + loss . R_{NN_t} + max( 0,\ s_{R} . portee . \sum_{J} . f_{2}) + \sum_{A} . f_{3} (1 - \frac{N}{K}) $$ 

$$IN_{NN_{t+1}} = IN_{NN_t}(1 - m_{1} - t{1} - \text{lat}) + \frac{\text{trans} . S_{NN_t} . II_{N_t}}{N}$$ 

$$II_{NN{t+1}} = II_{NN_t}(1 - m_{1} - madd - t_{1} - rec) + lat . IN_{NN_t}$$ 

$$R_{NN{t+1}} = R_{NN_t}(1 - m_{1} - t_{1} - loss) + rec . II_{NN_t}$$

**Classe d'âge J (jeunes)**

$$S_{J_{t+1}} = S_{NN_t}.t_{1} + S_{J_t}(1 - m_{2} - t_{2} - \frac{trans . II_{N_t}}{N}) + loss.R_{J_t}$$

$$IN_{J_{t+1}} = IN_{NN_t}.t_{1} + IN_{J_t}(1 - m_{2} - t_{2} - lat) + trans.S_{J_t}.II_{N_t}$$

$$II_{J_{t+1}} = II_{NN_t}.t_{1} + II_{J_t}(1 - m_{2} - madd - t_{2} - rec) + lat.IN_{J_t}$$

$$R_{J_{t+1}} = R_{NN_t}.t_{1} + R{J_t}(1 - m_{2} - t_{2} - loss) + rec.II_{J_t}$$

**Classe d'âge A (adultes)**

$$S_{A_{t+1}} = S_{J_t}.t_{2} + S_{A_t}(1 - m_{3} - \frac{trans.II_{N_t}}{N}) + loss.R_{A_t}$$

$$IN_{A_{t+1}} = IN_{J_t}.t_{2} + IN_{A_t}(1 - m_{3} - lat) + \frac{trans.S_{A_t}.II_{N_t}}{N}$$

$$II_{A_{t+1}}= II_{J_t}.t_{2} + II{A_t}(1 - m_{3} - madd - rec) + lat.IN_{A_t}$$

$$R_{A_{t+1}}=R_{J_t}.t_{2}+R_{A_t}(1-m_{3}-loss)+rec.II_{A_t}$$

**Total**
$$S_{N_{t+1}}=S_{NN_{t+1}}+S_{J_{t+1}}+S_{A_{t+1}}$$

$$IN_{N_{t+1}}=IN_{NN_{t+1}}+IN_{J_{t+1}}+IN_{A_{t+1}}$$

$$II_{N_{t+1}}=II_{NN_{t+1}}+II_{J_{t+1}}+II_{A_{t+1}}$$

$$R_{N_{t+1}}=R_{NN_{t+1}}+R_{J_{t+1}}+R_{A_{t+1}}$$

**Niveau d'"infectiosité" de l'agent pathogène par jour**

$$Nv_{inf_{t+1}}=trans.S_{N_t}(\frac{II_{N_t}}{N})$$

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


2. Analyse de sensibilité : Méthode OAT

3. Analyse de sensibilité : Méthode Morris

4. Analyse de sensibilité : Méthode FAST

Lien pour la méthode sur R : https://rdrr.io/cran/fast/man/sensitivity.html 


