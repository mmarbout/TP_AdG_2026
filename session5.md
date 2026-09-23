# Session 5

## Traitement et analyse de données RNAseq

on va commencer par récupérer des données issues d'une manip de RNAseq

```sh
scp sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2026-2027/HiC/session4/PAO_RNA_R1.fq.gz fastq/
```

vous pouvez vérifier que vous les avez bien récupérées 

```sh
ls fastq/
```

on va ensuite utiliser le package tinymapper qui contient tous les software nécessaire pour traiter ces données.

```sh
mkdir -p RNA_track/
```


```sh
micromamba activate tinymapper
```

on va utiliser le programme STAR pour l'alignement 

la première étape consiste à créer l'index du génome


```sh
STAR --runMode genomeGenerate --genomeFastaFiles ref/PAO.fa --genomeDir ref/PAO/STAR
```


ensuite on peut faire l'alignement (6-8 minutes)

```sh
STAR --genomeDir ref/PAO/STAR/ --readFilesCommand zcat --runThreadN 4 --readFilesIn fastq/MM301_R1.fq.gz --outFileNamePrefix RNA_track/MM301.bam. --outSAMtype BAM Unsorted --outSAMunmapped None --outSAMattributes Standard
```


il faut ensuite trier les reads, créer un index de ces reads (programme=samtools) (3-4 minutes)

```sh
samtools sort -@ 4 --output-fmt bam -l 9 -T RNA_track/bam/genome/MM301/MM301.bam_sorting -o RNA_track/bam/genome/MM301/MM301_filtered.bam RNA_track/MM301.bam.
samtools index -@ 4 RNA_track/MM301_filtered.bam
```

on peut ensuite créer nos profils de couverture de RNA avec le programme bamCoverage.

il y a énormément d'option à prendre en compte:

* --outFileFormat --> type de format de sortie. je vous conseille de mettre bedgraph car le format bigwig est un peu plus compliqué à prendre en main.
* --normalizeUsing --> CPM = Counts Per Million mapped reads ... le standard en RNAseq
* --binSize --> taille des bins pour lesquels on fait le calcul
* --filterRNAstrand --> doit on filtrer les reads en fonction de leur orientation ?


on va commencer doucement ... 


```sh
bamCoverage --bam RNA_track/MM301_filtered.bam --outFileFormat bedgraph --outFileName RNA_track/MM301_unstranded_bin1000.CPM.bed --binSize 1000 --numberOfProcessors 4  --normalizeUsing CPM --skipNonCoveredRegions  --ignoreDuplicates
```


si c'est bon pour vous , on va allez jeter un oeil à cela.


```sh
head RNA_track/MM301_unstranded_bin1000.CPM.bed
```


on va maintenant pouvoir traiter ces données sous R

ouvrez R studio.


```sh
library(ggplot2)

# import des données
track_rna=read.table("MM301_unstranded_bin1000.CPM.bed")

#creation du plot vaec ggplot
ggplot(track_rna,aes(x=track_rna$V2,y=track_rna$V4))+
+ geom_line()

```

évidemment, c'est mieux si on fait un beau plot avec des légendes ... 
mais ca on le garde pour juste après ;)


## exercice

refaites les fichiers bedgraph pour les gènes en forward et en reverse.
puis installer le package patchwork ou cowplot et essayer de me combiner ces deux plots ensemble avec la matrice d'interaction de la souche WT en phase expo.

--> https://r-graph-gallery.com/package/patchwork.html
--> https://r-graph-gallery.com/package/cowplot.html



<details><summary>Solution</summary>
<p>

```sh
install.packages("cowplot")
library(cowplot)

rna_for=read.table("MM301_for_bin1000.CPM.bed")
rna_rev=read.table("MM301_rev_bin1000.CPM.bed")

p0 <- plotMatrix(
	hic1,
	maxDistance = 200000,
	limits = c(0, 0.04),
	scale="linear",
	use.scores = "balanced",
	caption = FALSE
	) + coord_cartesian(xlim = c(1, 6271000))


p1<-ggplot(rna_for,aes(x=rna_for$V2,y=rna_for$V4))+
	geom_area(col = NA, fill = '#ee766f') +
	coord_cartesian(xlim = c(1, 6271000), expand = FALSE)+
	labs(x = NULL) +
	ylim(0,1000) +
	labs(y="RNA signal (CPM)")

p2<-ggplot(rna_rev,aes(x=rna_rev$V2,y=rna_rev$V4))+
	geom_area(col = NA, fill = '#83a8f7') +
	coord_cartesian(xlim = c(1, 6271000), expand = FALSE)+
	labs(x = "Genomic Coordinates") +
	ylim(0,1000) +
	labs(y="RNA signal (CPM)")

final_plot <- cowplot::plot_grid(p0, p1, p2, ncol = 1, align = "v", axis = "tb")

ggplot2::ggsave("plot/WT_RNA_5kb.pdf", plot = final_plot, width = 10, height = 8)

```

</p>
</detail>




