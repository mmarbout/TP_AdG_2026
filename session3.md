# Session 3

## analyse d'une matrice d'interaction sous R

### les fichiers matrices

on peut commencer à travailler sur nos fichiers matrices (i.e.: les fichiers cool).

```sh
library(HiCExperiment)
library(HiContacts)
library(GenomicRanges)
library(ggplot2)
library(dplyr)


coolf1 <-("cool_files/exemple.mcool")
cf1 <- CoolFile(coolf1)
```

Plusieurs « emplacements » (c'est-à-dire des éléments d'information) sont associés à un objet `ContactFile` :

* Le chemin d'accès à la matrice de contacts stockée sur disque
* La résolution active (par défaut, la résolution la plus fine disponible dans une matrice de contacts multi-résolution)
* Optionnellement, le chemin d'accès à un fichier de paires correspondantes
* Certaines métadonnées.

```sh
cf1
resolution(cf1)
pairsFile(cf1)
metadata(cf1)
availableResolutions(cf1)
availableChromosomes(cf1)
```

NB: Les objets ContactFile ne sont que des connexions à un fichier HiC stocké sur disque. Bien que des métadonnées soient disponibles, ils ne contiennent pas les données elles-mêmes !


on peut ensuite créer un object HiCExperiment a partir de cette connexion:

```sh
hic1 <- import(cf1, resolution=5000)
```

```sh
interactions(hic1)
```

il est ensuite possible de mettre ces données sous forme de data frame.

```sh
data <- as.data.frame(hic1)
```

on peut également importer les données uniquement pour une région donnée:

```sh
hic1_zoom <- import(cf1, resolution=1000, focus="E_coli:1-50000")
interactions(hic1_zoom)
```

### plot d'une matrice d'interaction

il existe ensuite une fonction pour visualiser directement la matrice:

```sh
plotMatrix(hic1)
```

```sh
plotMatrix(hic1_zoom)
```

et oui c'est aussi simple que cela !!! 
mais il existe pleins d'arguments à la fonction plotMatrix qui permettent de modifier l'image (voir l'aide).

je vous laisse jouer un peu avec la fonction et realiser différentes représentation de cette matrice d'interaction.

### comparaisons de deux matrices d'interactions

on peut également comparer deux matrices (à condition bien sur qu'elles aient été faites a partir du même génome)

```sh
coolf2 <-("cool_files/XX.mcool")
cf2 <- CoolFile(coolf2)
hic2 <- import(cf2, resolution=5000)
```

```sh
div_contacts <- divide(hic2, by = hic1)
plotMatrix(div_contacts,
    use.scores = 'balanced.fc', 
        scale = 'log2', 
        limits = c(-1, 1),
        cmap = bwrColors()
    )
```

### filtres des matrices d'interactions

sortez de R studio et retournez sur le terminal.
activez l'environnement hicstuff et refaites tourner le pipeline en activant l'option --filter

* Q: Quels filtres le pipeline a appliqué sur les données ?
* Q: Quels est le taux de reads conservées après le filtre de vos données ?

maintenant, retournez sous R et comparez vos deux matrices (unfilter vs. filter)





