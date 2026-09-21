# Session 6


## étude des intéractions entre molécules d'ADN

pour cette session, vous allez commencer par me récupérer un fichier cool sur GAIA.


```sh
mkdir -p cool_file/
scp sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/session5/PAO_pJN105.mcool cool_file/
```

vous allez ensuite me répondre aux questions suivantes:

* quelles résolutions sont disponibles dans votre fichiers mcool ?
* Combien de molécules d'ADN contient votre fichier mcool ?
* quelle est leur taille ?
* combien de contact intra-molécules contient chaque molécule d'ADN ?
* pouvez vous en déduire une abondance relative de chaque molécule ?


nous allons maintenant faire différentes analyses sur ce fichier ... ouvrez R studio et faites moi un plot de la matrice d'interaction dans son ensemble ainsi que de chaque molécule d'ADN.

<details><summary>Solution</summary>
<p>

```sh
library(HiCExperiment)
library(HiContacts)
library(GenomicRanges)
library(ggplot2)
library(dplyr)


coolf <-("cool_file/PAO_pJN105.mcool")
cf <- CoolFile(coolf)
hic <- import(cf, resolution=5000)
plotMatrix(hic)

hic_chr <- import(cf, resolution=5000, focus="LR657304.1")
plotMatrix(hic_chr)

hic_pls <- import(cf, resolution=1000, focus="pJN105")
plotMatrix(hic_pls)

```
</p>
</details>


nous allons maintenant voir comment visualiser les interactions d'une molécule d'ADN avec l'autre (4C plot).


```sh
v4C <- virtual4C(hic, viewpoint = GRanges("pJN105:1-6055"))
v4C


df <- as_tibble(v4C)
ggplot(df, aes(x = center, y = score)) + 
    geom_area(position = "identity", alpha = 0.5) + 
    theme_bw() + 
    labs(x = "Position", y = "Contacts with viewpoint") +
    scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-06)) + 
    facet_wrap(~seqnames, scales = 'free_y')
```
