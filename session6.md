# Session 6


## étude des intéractions entre molécules d'ADN

pour cette session, vous allez commencer par me récupérer un fichier cool sur GAIA.


```sh
mkdir -p cool_file/
scp sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/session5/V_cholerae.mcool cool_file/
```

vous allez ensuite me répondre aux questions suivantes:

* quelles résolutions sont disponibles dans votre fichiers mcool ?
* Combien de molécules d'ADN contient votre fichier mcool ?
* quelle est leur taille ?
* combien de contact intra-molécules contient chaque molécule d'ADN ?
* pouvez vous en déduire une abondance relative de chaque molécule ?

et faire également les exercices suivants:

* faites moi un plot de la matrice dans son ensemble.
* faites moi une matrice de chaque molécule d'ADN individuellement.


<details><summary>Solution</summary>
<p>

```sh
library(HiCExperiment)
library(HiContacts)
library(GenomicRanges)
library(ggplot2)
library(dplyr)


coolf <-("cool_file/V_cholerae.mcool")
cf <- CoolFile(coolf)
hic <- import(cf, resolution=5000)
plotMatrix(hic)

hic_chr1 <- import(cf, resolution=5000, focus="V_cholerae_chr1")
plotMatrix(hic_chr1)

hic_chr2 <- import(cf, resolution=5000, focus="V_cholerae_chr2")
plotMatrix(hic_chr2)

```
</p>
</details>


nous allons maintenant voir comment étudier les interactions entre deux molécules d'ADN.

on peut déjà commencer par faire une matrice des interactions inter-chromosomiques 

```sh
hic_interchr<- import(cf, resolution=5000, focus="V_cholerae_chr1|V_cholerae_chr2")
plotMatrix(hic_interchr)
```

On peut également analyser le profil d'interaction d'un locus génomique d'intérêt ou d'une molécule d'ADN avec son environnement immédiat ou avec le reste du génome. Dans certains cas, cela peut aider à identifier et/ou à comparer des interactions régulatrices ou structurelles. Cela peut également servir a analyser l'interaction d'une petite molécule d'ADN avec le génome de l'hôte.

Par exemple, il est possible de calculer le profil d'interaction « 4C virtuel » à l'échelle du génome, ancré au niveau du chromosome 2.

```sh
v4C <- virtual4C(hic, viewpoint = GRanges("V_cholerae_chr2:1-1072315"))
df <- as_tibble(v4C)
ggplot(df, aes(x = center, y = score)) + 
    geom_area(position = "identity", alpha = 0.5) + 
    theme_bw() + 
    labs(x = "Position", y = "Contacts with viewpoint") +
    scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-06)) + 
    facet_wrap(~seqnames, scales = 'free_y')
```


on peut également faire cela uniquement pour la région Origine (disons 1-20000) ou Terminus (disons 525000-545000) du chromosome 2

a vous de jouer


<details><summary>Solution</summary>
<p>

```sh

v4C_ORI <- virtual4C(hic, viewpoint = GRanges("V_cholerae_chr2:1-20000"))
df_ORI <- as_tibble(v4C_ORI)

ggplot(df_ORI, aes(x = center, y = score)) + 
    geom_area(position = "identity", alpha = 0.5) + 
    theme_bw() + 
    labs(x = "Position", y = "Contacts with viewpoint") +
    scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-06)) + 
    facet_wrap(~seqnames, scales = 'free_y')

v4C_TER <- virtual4C(hic, viewpoint = GRanges("V_cholerae_chr2:525000-545000"))
df_TER <- as_tibble(v4C_TER)

ggplot(df_TER, aes(x = center, y = score)) + 
    geom_area(position = "identity", alpha = 0.5) + 
    theme_bw() + 
    labs(x = "Position", y = "Contacts with viewpoint") +
    scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-06)) + 
    facet_wrap(~seqnames, scales = 'free_y')
```
</p>
</details>

