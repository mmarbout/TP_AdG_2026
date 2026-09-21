# Session 6


## étude des intéractions entre molécules d'ADN

pour cette session, vous allez commencer par me récupérer un fichier cool sur GAIA.


```sh
mkdir -p cool_file/
scp sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/session5/PAO_pJN105.mcool cool_file/
```

vous allez ensuite me répondre aux questions suivantes:

* Combien de molécules d'ADN contient votre fichier cool ?
* quelle est leur taille ?
* combien de contact intra-molécules contient chaque molécule d'ADN ?
* pouvez vous en déduire une abondance relative de chaque molécule ?


nous allons maintenant faire différentes analyses sur ce fichier ... ouvrez R studio et faites moi un plot de la matrice d'interaction.


```sh
library(HiCExperiment)
library(HiContacts)
library(GenomicRanges)
library(ggplot2)
library(dplyr)


coolf1 <-("HiC/exemple/exemple.mcool")
cf1 <- CoolFile(coolf1)

```
