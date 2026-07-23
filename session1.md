# Session 1

## mise en place de l'environnment et récupération des données

se placer sur le bureau de la Machine virtuelle

```sh
cd ~/Bureau/
```

créer un répertoire (l'option -p permet de créer des répertoires de manière recursive et dans des répertoires n'existant pas et évite les messages d'erreurs ... c'est parfois utile et c'est un réflexe chez moi)

```sh
mkdir -p TP_HiC/
```
rentrer dans le repertoire

```sh
cd  TP_HiC/
```

toutes les lignes de commande que vous verrez s'exécuteront depuis cet emplacement désormais !!!!


créer un répertoire pour y déposer les fichiers fastq et les fichiers log des différents logiciles

```sh
mkdir -p fastq/
mkdir -p log_files/
```

Vous allez travailler avec les fichiers de sorties de séquençage correspondants aux deux banques construites par votre binome: les reads en sens (forward) et en anti-sens (reverse) pour chaque banque construites. Vos fichiers sont nommés ainsi et se trouve sur l'espace GAIA:

* Binome_X_1_R1.fq.gz
* Binome_X_1_R2.fq.gz
* Binome_X_2_R1.fq.gz
* Binome_X_2_R2.fq.gz


récupérer les données correspondants à vos librairies en copiant les fichiers fastq (n'oubliez pas de changer le X !!!)

```sh
scp votrelogin@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/fastq/Binome_X_* fastq/
```

il faut également récupérer le fichier FastA correspondant à notre génome de référence.

```sh
mkdir -p ref/
scp votrelogin@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/ref/PAO.fa ref/
```



il nous reste une dernière chose à faire ... installer le pipeline hicstuff qui permet de traiter les données de HiC et construire les fichiers de matrice !!!

lien vers le github du programme:

[hicstuff](https://github.com/koszullab/hicstuff)

Pour cela nous allons utilisons micromamba qui est un système de gestion de paquets et d'environnement open-source qui fonctionne sous Windows, macOS et Linux. Micromamba installe, exécute et met à jour rapidement les paquets et leurs dépendances. Micromamba crée, enregistre, charge et bascule facilement entre les environnements sur votre ordinateur local.

lien vers le github du programme:

[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

```sh
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

une fois installé , nous pouvons vérifier que tout se passe bien en faisant un update de micromaba

```sh
micromamba self-update
```
il est alors possible de créer un environnement dédié à notre logiciel hicstuff.

```sh
micromamba create -n hicstuff bioconda::metator
```

il faut ensuite activer l'envrionnement pour pouvoir utiliser hicstuff

```sh
micromamba activate hicstuff
```

on peut vérifier facilement que tout est ok en appelant la doc de hicstuff.

```sh
hicstuff --help
```

pour désactiver l'environnement, il suffit de taper la commande suivant 

```sh
micromamba deactivate
```

gardez bien en tête qu'il faudra activer l'environnement à chaque fois que l'on veut utiliser hicstuff.


si tout est ok, alors on est prêt à commencer !!!





