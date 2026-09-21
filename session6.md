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


nous allons maintenant voir comment étudier les interactions entre deux molécules d'ADN.

Le rapport entre les interactions *cis* et les interactions *trans* est souvent utilisé pour évaluer les interactions entre molécules d'ADN. Il peut être calculé pour chaque molécule à l'aide de la fonction `cisTransRatio()`. Vous devrez fournir un objet `HiCExperiment` couvrant l'ensemble du génome pour estimer les rapports *cis*/*trans* !

```sh
hic <- import(cf, resolution = 1000)
ct <- cisTransRatio(hic) 
ct
```

Il peut être représenté graphiquement à l'aide de fonctions de visualisation basées sur ggplot2.


```sh
ggplot(ct, aes(x = chr, y = cis_pct)) + 
    geom_col(position = position_stack()) + 
    theme_bw() + 
    guides(x=guide_axis(angle = 90)) + 
    scale_y_continuous(labels = scales::percent) + 
    labs(x = 'Chromosomes', y = '% of cis contacts')
```

On peut également analyser le profil d'interaction d'un locus génomique d'intérêt ou d'une molécule d'ADN avec son environnement immédiat ou avec le reste du génome. Dans certains cas, cela peut aider à identifier et/ou à comparer des interactions régulatrices ou structurelles. Cela peut également servir a analyser l'interaction d'une petite molécule d'ADN avec le génome de l'hôte.

Par exemple, il est possible de calculer le profil d'interaction « 4C virtuel » à l'échelle du génome, ancré au niveau du plasmide pJN105.

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
