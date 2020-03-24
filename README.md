# QuadsPython

## Initiation (7 étapes) :

1) Charger un .obj
2) Afficher un .obj
3) Afficher/souligner le bord
4) Afficher le 1-ring d'un sommet
5) Afficher des croix avec un angles aléatoires sur les sommets
6) Afficher les contraintes sur le bord et 0 à l'interieur
7) Afficher un champ lissé à l'interieur

### Etape 1 et 2 : Charger et afficher un .obj

![Alt text](/screenshots/ObjCreator.png?raw=true "Optional Title")

### Etape 3 : Afficher/souligner le bord

![Alt text](/screenshots/surface_entouree.png?raw=true "Optional Title")

### Etape 4a : Afficher le 1-ring d'un sommet de l'intérieur de la surface

![Alt text](/screenshots/capturewith1ring.png?raw=true "Optional Title")

### Etape 4b : Afficher le 1-ring d'un sommet du contour de la surface

![Alt text](/screenshots/1ringSurunbord.png?raw=true "Optional Title")


### Etape 5 : Croix sur les sommets avec un angle aléatoire

![Alt text](/screenshots/croixaleatoires.png?raw=true "Optional Title")

### Etape 5b : Vecteurs normaux sur les milieux d'arêtes

![Alt text](/screenshots/vecteursnormauxsurcanard.png?raw=true "Optional Title")


### Etape 5b : Vecteurs normaux sur les sommets

![Alt text](/screenshots/deuxiemeessai.png?raw=true "Optional Title")


### Etape 6 : Afficher les contraintes sur le bord et 0 à l'interieur

![Alt text](/screenshots/contraintes.png?raw=true "Optional Title")


### Etape 7 : Afficher un champ lissé à l'interieur

![Alt text](/screenshots/deuxiemeessai.png?raw=true "Optional Title")
![Alt text](/screenshots/duck_example.png?raw=true "Optional Title")

### Etape 8 et 9 : Afficher un champ lissé grace à la résolution d'un système linéaire

A partir de l'article https://hal.inria.fr/hal-01245657/file/framefield.pdf :

En partant des équations 

√π(X[2i]−X[2j]) = 0

√π(X[2i+1]−X[2j+1]) = 0

On les regroupe en une seule en posant Z[i] = X[2i] + j * X[2i + 1] avec j le nombre imaginaire pur unitaire: 

√π(Z[i]−Z[j]) = 0

Et on remplace les contraintes de bords 

CX[2i] = C * cos(4θ[i])
CX[2i+1] = C * sin (4θ[i]) 

-> Z[i] = e(j * 4θ[i]) avec j le nombre imaginaire pur unitaire

Algorithme pour déterminer le système linéaire:

Pour tout i,

  si i est une contrainte, on ajoute "l'équation": Z[i] = e(j * 4θ[i])
  
  sinon on ajoute l'équation : nbVoisins * Z[i] - somme(pour k les voisins de i n'appartenant pas au contour){ Z[k] } = somme(pour k les voisins de i appartenant au contour) { e(j * 4θ[k]) }

Ce système linéaire est équivalent à la matrice AZ = b avec A une matrice carré inversible.
La résolution est donc simplement Z = A^(-1) * b

On retrouve [cos(phase(Z[i]) / 4), sin(phase(Z[i])/ 4)] les coordonnées du vecteur représentant le champs sur le sommet n°i

On peut donc obtenir une solution exacte plutôt qu'une solution approchée : 

![Alt text](/screenshots/interpolationlineairecomplexe.png?raw=true "Optional Title")
![Alt text](/screenshots/duckcomplexe.png?raw=true "Optional Title")

### Etape 10: Déterminer les singularités 

### Etape 11 : QuadCover 

Decomposition intuitive : on part des sommets du bord, et on suit les lignes de champs (sans se soucier des singularités):

![Alt text](/screenshots/exempledecomposition.png?raw=true "Optional Title")

