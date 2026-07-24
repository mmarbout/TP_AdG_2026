# Session 2

## génération du fichier matrice

le pipeline hicstuff permet de générer, à partir d'un génome (fasta) et de données de séuqnçage (HiC), un fihcier mcool (multicool) qui est le format standard des données HiC. 
[cooler package](https://github.com/open2c/cooler)

voici un schéma du pipeline:

![hicstuff](docs/images/pipeline.svg)

on n'oublie pas de se mettre dans le bon répertoire et d'activer l'environnement

```sh
cd ~/Bureau/TP_HiC/
micromamba activate hicstuff
```

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
hicstuff pipeline --genome ref/PAO.fa --binning 1000 --distance-law --duplicates --enzyme DpnII,HinfI --filter --outdir hic/binome1_1/ --plot --prefix binome1_1 --threads 4 --skip-count fastq/Binome1_1_R1.fq.gz fastq/Binome1_1_R2.fq.gz
```

cela ne devrait pas prendre plus de 10 minutes ...  

maintenant que c'est fait , vous pouvez regarder où en est votre pipeline hicstuff et explorer les fichiers de sorties.


```sh
ls HiC/binome1_1/
```

je vous laisse explorer tout ca et répondre aux questions suivantes:

* Q: Combien de reads aviez vous dans votre librairies ?
* Q: Quel est le taux de mapping de vos données sur le génome de référence ?
* Q: Quel est le taux de duplicats de PCR ?
* Q: Quels est le taux de reads conservées après le filtre de vos données ?
* Q: combien votre matrice contient de contacts ?




