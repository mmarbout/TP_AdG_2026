# Session 2

## génération du fichier matrice

le pipeline hicstuff permet de générer, à partir d'un génome (fasta) et de données de séuqnçage (HiC), un fihcier mcool (multicool) qui est le format standard des données HiC. 
[cooler package](https://github.com/open2c/cooler)

voici un schéma du pipeline:

![hicstuff](docs/images/pipeline.svg)


il peut se réaliser par étape ou d'un seul coup.

```sh
hicstuff --help
```

nous allons lancer la commande pipeline pour lancer l'ensemble du pipeline en une seule fois.

```sh
hicstuff pipeline --help
```

les arguments à donner obligatoirement sont le génome (ou l'index), les fichiers fastq. Il y a également toute une série d'argument optionnel que l'on peut fournir au pipeline.

![hicstuff_pipeline](docs/images/hicstuff_pipeline.png)

voici la commande à lancer (c'est un exemple pour le jeu de donnée Binome_1_1):

```sh
hicstuff pipeline --genome ref/PAO.fa \
				--binning 1000 \
				--distance_law \
				--duplicates \
				--enzyme DpnII,HinfI \
				--filter \
				--outdir HiC/binome_1_1/ \
				--plot \
				--prefix binome_X_1 \
				--threads 4 \
				--zoomify \
				--skip-count \
				fastq/Binome_1_1_R1.fq.gz \
				fastq/Binome_1_1_R2.fq.gz
```

cela risque de prendre un peu de temps ... 

donc pendant ce temps la nous allons installer les packages R nécessaire à l'analyse de nos données et faire quelques exercices sur des données tests.


## analyse d'une matrice d'interaction

commencer par copier le fichier suivant sur gaia

```sh
mkdir -p cool_files/
scp votrelogin@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/cool_files/exemple1.mcool cool_files/
scp votrelogin@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/cool_files/exemple2.mcool cool_files/
```
lancer R studio et mettez vous dans le répertoire adéquat (TP_HiC) puis installez les packages suivant:

```sh
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
```

```sh
BiocManager::install("HiCExperiment", ask = FALSE)
BiocManager::install("HiCool", ask = FALSE)
BiocManager::install("HiContacts", ask = FALSE)
BiocManager::install("GenomicRanges", ask = FALSE)
```

on peut ensuite commencer à travailler sur nos fichiers cool.

```sh
library(HiCExperiment)
library(HiContacts)
coolf <-("cool_files/exemple1.mcool)
cf <- CoolFile(coolf)
```

Plusieurs « emplacements » (c'est-à-dire des éléments d'information) sont associés à un objet `ContactFile` :

* Le chemin d'accès à la matrice de contacts stockée sur disque
* La résolution active (par défaut, la résolution la plus fine disponible dans une matrice de contacts multi-résolution)
* Optionnellement, le chemin d'accès à un fichier de paires correspondantes
* Certaines métadonnées.


```sh
cf
resolution(cf)
pairsFile(cf)
metadata(cf)
availableResolutions(cf)
availableChromosomes(cf)
interactions(hic)
```

on peut également importer les données uniquement pour une région donnée:

```sh
plotMatrix(hic)
```

```sh
plotMatrix(hic)
```

et oui c'est aussi simple que cela !!! 
mais il existe pleins d'arguments à la fonction plotMatrix qui permettent de modifier l'image (voir l'aide).







## analyse de vos matrice d'interaction

maintenant que c'est fait , on peut regarder ou en est notre pipeline hicstuff et explorer les fichiers de sorties.



```sh
ls HiC/binome_1_1/
```

je vous laisse explorer tout ca et répondre aux questions suivantes:

* Q: Combien de reads aviez vous dans votre librairies ?
* Q: Quel est le taux de mapping de vos données sur le génome de référence ?
* Q: Quel est le taux de duplicats de PCR ?
* Q: Quels est le taux de reads conservées après le filtre de vos données ?
* Q: combien votre matrice contient de contacts ?
