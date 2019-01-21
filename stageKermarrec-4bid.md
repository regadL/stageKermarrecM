---
title: "Stage de Maxime Kermarrec (L3-BI)"
author: "Leslie REGAD et Maxime Kermarrec"
date: '2019-01-21'
output:
  html_document:
    code_folding: show #hide
    self_contained: yes
    fig_caption: yes
    highlight: pygments #pour les sorties R
    theme: spacelab
    toc: yes  #for add table of contents
    toc_depth: 3
    toc_float: yes
    df_print: paged
    keep_md: yes  #pour garder le .md after run
  slidy_presentation:
    smart: no
    slide_level: 2
    self_contained: yes
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    highlight: tango
    incremental: no
    keep_md: yes
    smaller: yes
    theme: cerulean
    toc: yes
    widescreen: yes
  ioslides_presentation:
    slide_level: 2
    self_contained: no
    colortheme: dolphin
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    fonttheme: structurebold
    highlight: tango
    smaller: yes
    toc: yes
    widescreen: yes
  pdf_document:
    fig_caption: yes
    highlight: zenburn
    toc: yes
    toc_depth: 3
  beamer_presentation:
    colortheme: dolphin
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    fonttheme: structurebold
    highlight: tango
    incremental: no
    keep_tex: no
    slide_level: 2
    theme: Montpellier
    toc: yes
font-import: http://fonts.googleapis.com/css?family=Risque
subtitle: DTest
font-family: Garamond
transition: linear
---




# Objectif

* Travail 1 :

1. Etudier la conservation du packing cristallin dans les 13 structures de PR2 et les 15 structures de PR1
2. Mettre en relation cette conservation avec les RMSD
3. Lien entre conservation de l’asymétrie et neqSL ?

* Travail 2 : Comparaison des interfaces PR1 et PR2
* Travail 3 : Etude de la variabilité structurale de PR1 et PR2 – SA-conf
* Travail 4 : Etude de la flexibilité de PR1 et PR2 – distance au SL


# Semaine 1 : Vendredi 18/01/2018

## Objectifs 
Etudier la conservation du packing cristallin dans les 13 structures de PR2 et les 15 structures de PR1.
Pour cela, on va identifier les résidus qui sont : 

* inclus dans le packing cristallin chez toutes (ou la majorité (=80%)) des structures
* inclus dans le packing cristallin chez toutes (ou la majorité (=80%)) des structures de PR1
* inclus dans le packing cristallin chez toutes (ou la majorité (=80%)) des structures de PR2


## Données

* 26 fichiers qui contient les atomes impliqués dans le packing cristallin pour chaque structure de PR1 et PR2
* le type de PR pour chaque code PDB

PDB code | type | remarque
---------|------|---------
3s45 | PR2 |
1hhp | PR1 |
1hih | PR1 |
1hii | PR2 |
1hiv | PR1 |
1hpv | PR1 |
1hsh | PR2 |
1hsi | PR2 | 
1ivp | PR2 |
1sdt | PR1 | (peu avoir des pbl de numérotation des résidus)
2hb3 | PR1 |
2hb4 | PR1 | (monomère)
2hpe | PR2 |
2hpf | PR2 |
2ien | PR1 |  (peu avoir des pbl de numérotation des résidus)
2mip | PR2 |
2nph | PR1 |  (peu avoir des pbl de numérotation des résidus)
3phv | PR1 | (monomère)
2z4o | PR1 |  (peu avoir des pbl de numérotation des résidus)
3ebz | PR2 |  (peu avoir des pbl de numérotation des résidus)
3ec0 | PR2 |  (peu avoir des pbl de numérotation des résidus)
3ecg | PR2 | (peu avoir des pbl de numérotation des résidus)
3ekv | PR1 |
3nu3 | PR1 |  (peu avoir des pbl de numérotation des résidus)
4hla | PR1 |
4ll3 | PR1 |
1hsh | PR2 | 
2pc0 | PR1 | 

## To do

### Création d'une matrice sous R qui contient en lignes les protéines et en colonnes les positions

La case de la matrice contient : 

* 1 quand le résidu a au moins 1 atomes impliqués dans le packing dans la structure, 
* 0 sinon.


Attention, 

* certaines structures sont des monomères 
* certaines structures ont des mauvaises numérotations des résidus. Normalement les résidus sont numérotés de 1 à 99 pour la chaîne A et la chaîne B. Certaines structures ont les résidus de la chaîne B numérotés de 101 à 109.
Ainsi le résidus 101B = 1B


1. Etape 1 : récupération de la liste des fichiers PDB

```r
listFile = dir("fileByProt3/")
```

2. Exemple sur le premier fichier


```r
filein = listFile[1]
M = read.table(paste("fileByProt3",filein,sep="/"))
listAtom = unique(paste(as.character(M[,6]), as.character(M[,5]), sep="_"))
listAtomSyn = listAtom
listAtomSyn
```

```
 [1] "1_A"  "2_A"  "3_A"  "4_A"  "5_A"  "6_A"  "7_A"  "8_A"  "9_A"  "23_A" "24_A" "25_A" "26_A" "27_A" "29_A" "35_A" "37_A" "42_A" "48_A" "49_A" "50_A" "51_A" "52_A" "53_A" "54_A" "55_A" "56_A" "57_A" "61_A" "67_A" "69_A" "72_A" "79_A" "87_A" "90_A" "91_A" "92_A" "93_A" "94_A" "95_A" "96_A" "97_A"
[43] "98_A" "99_A"
```


Récupération de la liste des résidus pour chaque protéine


```r
v1 =(1:26)
v2 =(1:160)

for (i in v1){
  for (j in v2) {
  filein = listFile[i]
  M = read.table(paste("fileByProt3",filein,sep="/"))
  listAtom = unique(paste(as.character(M[,6]), as.character(M[,5]), sep="_"))
  listAtomSyn = unique(c(listAtom,listAtomSyn))
  }
  
}

listAtomSyn[160] = ("Protease")
length(listAtomSyn)
```

```
[1] 160
```

```r
length(listAtom)
```

```
[1] 69
```
*Remarques* : 
Je ne comprends pas pourquoi dans votre code ci-dessus, il y a une double boucle sachant que l'indice `j` vous ne l'utilisez pas.

Faire la matrice

*Remarques*

Pour simplifier la matrice, il faudrait donner des acronymes aux trois dernières lignes pour simplifier les rownames : 
* NbStructinPC : Nb de structure dans lequel le résidu est implique",
* NbStructinPC.PR1 :  "Nb de structure PR1 dans lequel le résidu est implique", 
* NbStructinPC.PR2 : "Nb de structure PR2 dans lequel le résidu est implique.



```r
v3 = c("1hhp","1hih","1hii","1hiv","1hpv","1hsh","1hsi","1ivp","1sdt","2hb3","2hb4","2hpe","2hpf","2ien","2mip","2nph","2z4o","3ebz","3ec0","3ecg","3ekv","3nu3","3phv","3s45","4hla","4ll3","* NbStructinPC", "NbStructinPC.PR1", "NbStructinPC.PR2Nb")


matrice <- matrix(c(rep(0,(length(listAtomSyn)))*(length(v3))), nrow=length(v3), ncol=length(listAtomSyn))

#vous pouvez simplifier avec la commande suivante : 
matrice <- matrix(0, nrow=length(v3), ncol=length(listAtomSyn))


rownames(matrice) <- c(v3)
colnames(matrice) <- c(listAtomSyn)
#matrice
#dans les lignes ci-dessus les fonctions c ne sont pas nécessaires


filein = listFile[2]
M = read.table(paste("fileByProt3",filein,sep="/"))
Atome <- c(paste(as.character(M[,6]), as.character(M[,5]), sep="_"))
Atome
```

```
  [1] "2_A"  "2_A"  "2_A"  "2_A"  "2_A"  "2_A"  "4_A"  "6_A"  "6_A"  "6_A"  "6_A"  "6_A"  "12_A" "12_A" "14_A" "14_A" "14_A" "16_A" "16_A" "16_A" "18_A" "18_A" "18_A" "18_A" "19_A" "19_A" "34_A" "34_A" "34_A" "34_A" "34_A" "37_A" "37_A" "37_A" "37_A" "38_A" "38_A" "38_A" "38_A" "39_A" "39_A" "39_A"
 [43] "39_A" "39_A" "39_A" "40_A" "42_A" "42_A" "42_A" "43_A" "43_A" "43_A" "44_A" "44_A" "44_A" "46_A" "46_A" "46_A" "46_A" "53_A" "53_A" "53_A" "53_A" "53_A" "53_A" "53_A" "55_A" "55_A" "61_A" "61_A" "61_A" "63_A" "63_A" "63_A" "67_A" "67_A" "67_A" "67_A" "68_A" "68_A" "70_A" "72_A" "72_A" "72_A"
 [85] "72_A" "72_A" "72_A" "72_A" "79_A" "79_A" "79_A" "79_A" "80_A" "80_A" "81_A" "81_A" "81_A" "81_A" "81_A" "81_A" "81_A" "82_A" "91_A" "91_A" "91_A" "91_A" "92_A" "92_A" "92_A" "92_A" "92_A" "92_A" "92_A" "4_B"  "4_B"  "4_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B"  "6_B" 
[127] "7_B"  "7_B"  "7_B"  "7_B"  "12_B" "12_B" "14_B" "16_B" "16_B" "16_B" "17_B" "17_B" "17_B" "17_B" "18_B" "19_B" "19_B" "19_B" "19_B" "19_B" "19_B" "19_B" "20_B" "20_B" "20_B" "35_B" "35_B" "35_B" "35_B" "35_B" "36_B" "37_B" "37_B" "37_B" "37_B" "38_B" "38_B" "39_B" "39_B" "39_B" "39_B" "39_B"
[169] "39_B" "39_B" "40_B" "40_B" "41_B" "41_B" "41_B" "41_B" "41_B" "42_B" "42_B" "42_B" "42_B" "42_B" "42_B" "42_B" "42_B" "43_B" "43_B" "44_B" "44_B" "44_B" "44_B" "44_B" "45_B" "45_B" "45_B" "45_B" "45_B" "45_B" "46_B" "46_B" "46_B" "46_B" "49_B" "50_B" "52_B" "52_B" "53_B" "53_B" "53_B" "53_B"
[211] "53_B" "53_B" "55_B" "55_B" "57_B" "58_B" "63_B" "63_B" "63_B" "70_B" "70_B" "70_B" "70_B" "70_B" "70_B" "71_B" "72_B" "79_B" "79_B" "79_B" "79_B" "80_B" "80_B" "81_B" "81_B" "81_B" "81_B" "81_B" "92_B" "92_B" "94_B" "96_B"
```

```r
for (i in 1:length(v1)) {
  filein = listFile[i]
  M = read.table(paste("fileByProt3",filein,sep="/"))
  Atome <- unique(c(paste(as.character(M[,6]), as.character(M[,5]), sep="_")))
  for (j in 1:length(v2)) {
    for (y in 1:length(Atome)) {
      if ((listAtomSyn[j]) == (Atome[y])){
        matrice[i,j] = 1
      }
    }
  }
} 




#residu = function(n){

#matrice
```


Solution 2 : 


```r
v3 = c("1hhp","1hih","1hii","1hiv","1hpv","1hsh","1hsi","1ivp","1sdt","2hb3","2hb4","2hpe","2hpf","2ien","2mip","2nph","2z4o","3ebz","3ec0","3ecg","3ekv","3nu3","3phv","3s45","4hla","4ll3","* NbStructinPC", "NbStructinPC.PR1", "NbStructinPC.PR2Nb")

#creation de la matrice contenant que des 0 
matrice2 <- matrix(0, nrow=length(v3), ncol=length(listAtomSyn))
rownames(matrice2) <- v3
colnames(matrice2) <- listAtomSyn
#matrice
#dans les lignes ci-dessus les fonctions c ne sont pas nécessaires



for (i in 1:length(v1)) {
  filein = listFile[i]
  M = read.table(paste("fileByProt3",filein,sep="/"))
  residue <- unique(c(paste(as.character(M[,6]), as.character(M[,5]), sep="_")))
  pdbCode <- unlist(strsplit(filein,"_"))[1]
  matrice2[pdbCode, residue] = 1
} 
```


*Remarques*

* je pense qu'il serait préférable de stocker ces informations : "Nb de structure dans lequel le résidu est implique", "Nb de structure PR1 dans lequel le résidu est implique" et "Nb de structure PR2 dans lequel le résidu est implique", dans une autre matrice



### Coloration de la matrice
Utilisation la commande `pheatmap` du package `pheatmap`


```r
library(pheatmap)
?pheatmap
pheatmap::pheatmap(matrice)
```

<img src="figures/07_tests_multiplescolorMat-1.png" style="display: block; margin: auto;" />

```r
pheatmap(matrice[-27:-29,], cluster_rows = TRUE, cluster_cols = FALSE, br=-1:1, col=c("white", "red"))
```

<img src="figures/07_tests_multiplescolorMat-2.png" style="display: block; margin: auto;" />
*Remarques:*

* Pour cette première visuation, je ne mettrais pas les lignes trois dernières lignes de la matrice pour faciliter la lecture*

* Il faudrait aussi trier les colonnes de la matrice en donnant : 
- les résidus de la chaîne A
- puis les résidus de la chaîne B


### Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin 
Le  nombre de structure dans lequel un résidu est impliqué dans le packing cristallin = la somme des colonnes de la matrice

1. Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin

*Remarques* :
La je ne comprendre pas bien ce que vous faites...



```r
for (i in 1:length(v3)){
    if (matrice[i,] == "3s45"|| matrice[i,] == "1hii"|| matrice[i,] == "1hsh"|| matrice[i,] == "1hsi"|| matrice[i,] == "1ivp"|| matrice[i,] == "3s45"|| matrice[i,] == "2hpe"|| matrice[i,] == "2hpf"|| matrice[i,] == "2mip"|| matrice[i,] == "3ebz"|| matrice[i,] == "3ec0"|| matrice[i,] == "3ecg"){
      
      matrice[i,160] = 2
    }
    if (matrice[i,] == "1hhp" || matrice[i,] == "1hih" || matrice[i,] == "1hiv" || matrice[i,] == "1hpv" || matrice[i,] == "1sdt" || matrice[i,] == "2hb3" || matrice[i,] == "2hb4" || matrice[i,] == "2ien" || matrice[i,] == "2nph" || matrice[i,] == "3phv" || matrice[i,] == "2z4o" || matrice[i,] == "3ekv" || matrice[i,] == "3nu3" || matrice[i,] == "4hla" || matrice[i,] == "4ll3" || matrice[i,] == "2pc0"){
      
    matrice[i,160] = 1
    }
  
}
matrice
```

```
                   2_A 4_A 6_A 7_A 8_A 12_A 14_A 16_A 17_A 18_A 19_A 20_A 21_A 35_A 37_A 38_A 39_A 40_A 43_A 44_A 45_A 46_A 47_A 48_A 51_A 52_A 53_A 60_A 61_A 63_A 70_A 71_A 72_A 79_A 81_A 92_A 94_A 96_A 2_B 4_B 6_B 7_B 8_B 12_B 14_B 16_B 17_B 18_B 19_B 20_B 21_B 37_B 38_B 39_B 40_B 44_B 46_B 47_B
1hhp                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0    0    0    1    1    1    1    0    1    0    0    0    1    1    0    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1hih                 1   1   1   0   0    1    1    1    0    1    1    0    0    0    1    1    1    1    1    1    0    1    0    0    0    0    1    0    1    1    1    0    1    1    1    1    0    0   0   1   1   1   0    1    1    1    1    1    1    1    0    1    1    1    1    1    1    0
1hii                 1   1   1   1   0    0    1    1    1    1    1    0    0    0    0    0    1    1    1    0    1    1    0    0    0    0    1    1    1    0    1    0    0    0    0    1    1    0   1   1   1   1   0    0    0    0    0    0    0    0    0    1    1    1    1    1    1    0
1hiv                 1   1   1   0   0    1    1    1    0    1    1    0    1    0    1    1    1    0    0    1    1    1    0    0    0    1    1    0    1    1    1    1    1    1    1    1    0    0   0   1   1   1   0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    0
1hpv                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    1    0    0    0    1    1    0    1    1    1    1    1    0    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    0    1    1
1hsh                 1   0   1   1   1    1    0    1    1    0    0    0    0    1    1    0    0    1    1    1    1    1    0    0    1    0    0    1    1    1    0    0    1    1    1    0    0    0   0   0   1   1   1    0    0    1    1    1    1    0    0    0    1    1    0    1    1    0
1hsi                 1   1   1   1   0    0    0    1    1    1    0    0    0    0    0    0    1    1    0    0    0    1    0    0    0    0    1    0    0    1    0    0    1    0    0    0    0    0   1   0   0   1   1    1    0    1    1    1    1    0    1    1    1    1    0    1    1    0
1ivp                 1   1   1   1   0    1    1    0    0    0    1    1    0    0    1    0    1    0    0    1    1    1    0    0    0    0    1    1    0    0    0    0    0    1    1    1    1    0   0   0   1   1   0    0    1    0    1    0    1    0    0    0    0    1    1    1    0    0
1sdt                 0   1   1   1   0    0    0    0    0    1    0    0    0    1    1    0    0    0    1    1    1    1    0    0    0    0    1    0    0    0    0    0    1    0    0    1    1    0   1   0   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
2hb3                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    1    0    0    1    0    1    0    1    0    1    0    0    1    1    0   1   0   1   0   0    1    1    0    1    1    1    0    0    0    0    0    1    1    0    0
2hb4                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    1    1    1    0    1    0    0    0    1    1    1    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
2hpe                 1   1   1   1   0    0    0    0    0    0    1    0    0    1    1    0    1    1    0    1    0    0    0    0    1    1    0    1    1    0    0    0    0    1    1    0    0    0   1   0   1   1   0    1    0    0    1    1    1    1    1    1    0    1    1    0    1    0
2hpf                 0   0   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    1    1    0    1    1    0    0    0    0    0   0   0   1   1   0    1    1    1    1    0    1    1    1    1    0    1    0    0    1    0
2ien                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    0    0    1    1    0    1    0    1    0    1    0    0    1    1    0   0   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
2mip                 0   1   1   1   1    0    1    1    1    0    1    0    0    0    0    0    0    1    0    1    1    1    1    0    1    0    1    1    0    0    0    0    1    0    1    0    0    0   1   0   1   1   1    1    1    1    1    1    1    0    1    1    0    0    1    1    1    1
2nph                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    1    1    1    0    1    1    0    1    1    1    1    1    1    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0
2z4o                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    0    0    1    1    0    0    0    0    0    1    0    0    1    1    0   1   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
3ebz                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    1    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ec0                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    1    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    1    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ecg                 1   1   1   1   0    0    0    0    0    1    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    0    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ekv                 1   1   1   0   0    1    1    1    1    1    1    0    0    0    1    1    1    1    1    1    0    1    0    0    0    1    1    0    1    1    1    1    1    1    1    1    0    0   0   1   1   1   0    0    0    1    1    1    1    0    0    1    1    1    1    1    1    0
3nu3                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    1    0    0    1    0    0    0    0    0    1    0    0    1    1    0   1   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
3phv                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0    0    0    1    1    1    0    0    1    0    0    0    1    1    0    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3s45                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    0    1    0    0    1    0    0    0    0   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
4hla                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
4ll3                 1   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   1   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
* NbStructinPC       0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR1     0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR2Nb   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                   48_B 52_B 53_B 61_B 63_B 70_B 71_B 72_B 79_B 92_B 96_B 41_A 35_B 43_B 51_B 57_B 42_A 54_A 55_A 57_A 58_A 65_A 98_A 99_A 10_B 41_B 42_B 45_B 55_B 81_B 94_B 98_B 99_B 1_A 3_A 5_A 9_A 24_A 25_A 26_A 27_A 29_A 49_A 50_A 56_A 67_A 69_A 87_A 90_A 91_A 93_A 95_A 97_A 34_A 36_A 73_A 88_A
1hhp                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
1hih                  0    1    1    0    1    1    1    1    1    1    1    0    1    1    0    1    1    0    1    0    0    0    0    0    0    1    1    1    1    1    1    0    0   0   0   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    0    0    0    1    0    0    0
1hii                  0    0    1    1    0    0    0    1    1    0    0    1    0    1    0    1    1    0    1    0    0    0    0    0    0    1    1    1    0    1    0    0    1   1   1   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    0    0    0    0    0
1hiv                  0    1    1    0    1    1    0    1    1    1    1    1    1    0    0    1    1    0    1    0    0    0    0    0    0    1    1    1    0    1    1    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    0
1hpv                  1    1    1    1    1    1    1    1    0    1    0    0    1    1    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1hsh                  0    1    1    1    1    1    0    0    1    0    0    1    0    0    1    1    1    0    1    1    0    0    0    1    0    0    1    1    1    0    0    0    0   0   0   0   0    0    0    0    0    1    0    0    0    0    0    1    0    0    0    0    0    1    0    0    0
1hsi                  0    0    1    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1ivp                  0    0    1    1    0    0    0    1    1    0    0    1    0    0    0    1    1    0    1    0    0    1    0    0    0    0    1    0    0    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    1    1    0    0    0    0    0    0
1sdt                  0    0    0    1    1    0    1    1    1    1    0    0    0    1    0    0    0    0    0    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    1    1
2hb3                  0    0    0    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    1    0    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    0    0    0    0    0    0    1    0    0    0    1    1    1    1
2hb4                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
2hpe                  0    1    1    1    0    0    0    1    0    0    1    1    0    1    1    1    1    0    1    1    1    0    1    1    0    1    1    0    1    0    1    1    0   1   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0
2hpf                  0    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    1    0    0    0    0   0   0   1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
2ien                  0    0    0    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    1    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    1    0    0    0    0    0    1    0    0    0    1    1    1    1
2mip                  1    1    1    1    1    1    0    1    1    0    0    1    0    1    1    0    1    1    1    1    1    0    0    1    0    1    1    0    1    0    0    0    0   0   0   0   0    0    0    0    0    1    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0
2nph                  0    1    1    1    1    1    1    1    1    1    0    0    1    1    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0   0   1   0   0    0    0    0    0    0    1    0    0    1    0    0    0    0    0    0    0    0    0    0    0
2z4o                  0    0    0    1    1    1    1    1    1    1    0    0    0    1    0    0    0    0    0    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    1    0    0    0    0    0    1    0    0    0    1    1    1    1
3ebz                  1    0    1    1    0    0    0    0    1    0    1    1    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3ec0                  0    0    1    1    1    1    0    0    1    0    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0
3ecg                  1    0    1    0    0    1    0    0    1    0    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3ekv                  1    1    1    0    1    1    1    1    1    1    1    0    1    0    0    1    1    0    1    0    0    1    0    0    0    0    1    1    0    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    0    0    0    1    1    0    0
3nu3                  0    0    0    1    1    1    1    1    1    1    0    0    0    1    0    0    0    0    1    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    0    0    0    0    0    0    1    0    0    0    1    1    1    1
3phv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
3s45                  1    0    1    0    0    0    0    0    1    0    0    1    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
4hla                  1    1    1    1    1    1    1    1    1    1    0    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
4ll3                  1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
* NbStructinPC        0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR1      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR2Nb    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                   89_A 34_B 65_B 68_B 80_B 91_B 15_A 68_A 80_A 36_B 49_B 58_B 62_B 78_B 1_B 54_B 95_B 60_B 11_B 82_B 3_B 30_A 82_A 29_B 30_B 76_B 13_B 74_B 83_B 87_B 88_B 59_A 74_A 50_B 56_B 59_B 23_A 66_A 11_A 69_B 73_B 67_B 10_A 78_A Protease
1hhp                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0        0
1hih                  0    0    0    0    1    0    0    1    1    1    1    1    0    0   0    0    0    0    0    0   0    0    1    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0        0
1hii                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   1    0    0    1    0    0   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    1    0    0    1    0    0        0
1hiv                  0    0    0    0    1    0    0    1    1    0    1    0    0    0   0    0    0    0    0    0   0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    0        0
1hpv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
1hsh                  0    0    1    1    0    0    0    1    0    0    0    0    0    0   0    1    0    0    0    0   0    0    1    1    1    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    0    0    1    1        0
1hsi                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   1    0    0    0    0    0   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0        0
1ivp                  0    0    0    1    0    0    0    1    1    0    0    1    0    0   0    0    0    1    0    0   0    0    0    0    0    0    0    1    0    0    0    1    1    0    0    1    0    0    1    1    1    0    0    0        0
1sdt                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2hb3                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2hb4                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0        0
2hpe                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   0    0    1    1    0    1   0    0    0    0    0    0    0    1    1    0    0    1    1    1    1    1    0    0    0    0    0    0    0    0        0
2hpf                  0    0    0    0    0    1    0    0    0    0    0    0    0    0   0    1    0    1    0    0   0    0    0    0    0    0    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2ien                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2mip                  0    0    1    1    0    0    0    0    0    0    0    0    0    0   1    1    0    0    0    0   1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2nph                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    1   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2z4o                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ebz                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   1    1    0    1    1    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ec0                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   1    1    0    1    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ecg                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   1    1    1    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ekv                  0    0    0    0    1    0    1    1    1    1    1    1    1    1   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3nu3                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3phv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3s45                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
4hla                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
4ll3                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
* NbStructinPC        0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
NbStructinPC.PR1      0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
NbStructinPC.PR2Nb    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
```

```r
for (i in length(listAtomSyn)){
  y = 0
  for (j in length(v3)) {
    ma = matrice[j,i]
    if (ma == 1){
        y = y+1
    }
  ma = 0
  }
  matrice[length(v3)-2,i] = y
}
```


```r
NbStructinPC = apply(matrice[1:26,],2,sum)
sort(NbStructinPC)
```

```
Protease     15_A     62_B     78_B     11_B     30_A     76_B     13_B     88_B     56_B     66_A     73_B     78_A     36_B     95_B     29_B     30_B     83_B     87_B     59_A     74_A     50_B     23_A     69_B     67_B     10_A     47_A      9_A     24_A     25_A     26_A     27_A     56_A 
       0        1        1        1        1        1        1        1        1        1        1        1        1        2        2        2        2        2        2        2        2        2        2        2        2        2        3        3        3        3        3        3        3 
    90_A     95_A     97_A     49_B     58_B     74_B     59_B     11_A      8_B     47_B     51_B     10_B      5_A     80_A      3_B     82_A     21_A     98_B     99_B      1_A      3_A     29_A     50_A     87_A     93_A     73_A     88_A     89_A     34_B     82_B      8_A     60_A     35_B 
       3        3        3        3        3        3        3        3        4        4        4        4        4        4        4        4        5        5        5        5        5        5        5        5        5        5        5        5        5        5        6        6        6 
    65_A     69_A     91_B      1_B     54_B     60_B     38_A     71_A     96_A     94_B     67_A     36_A     65_B     51_A     48_B     96_B     54_A     98_A     49_A     80_B     68_A     12_A     17_A     20_A     81_A     34_A     14_A     38_B     52_B     99_A     16_A     19_A     39_A 
       6        6        6        6        6        6        7        7        7        7        7        7        7        8        8        8        8        8        8        8        8        9        9        9        9        9       10       10       10       10       11       11       11 
    48_A     94_A      4_B     16_B     71_B     57_B     58_A     45_B     20_B     57_A     52_A     63_A     21_B     39_B     92_B     91_A     35_A     70_A     41_A     41_B     68_B     18_A     40_A     45_A     70_B     42_A     55_B     61_A     79_A      2_B     18_B     37_B     61_B 
      11       11       11       11       11       11       11       11       12       12       13       13       13       13       13       13       14       14       14       14       14       15       15       15       15       15       15       16       16       16       16       16       16 
    63_B     72_B     43_B     81_B     43_A     92_A     14_B     46_B     42_B     37_A      6_B      7_B     53_B      2_A     44_A     12_B     79_B     40_B     44_B     55_A     72_A     46_A     17_B     19_B      7_A     53_A      4_A      6_A 
      16       16       16       16       17       17       17       17       17       18       18       18       18       19       19       19       19       20       20       20       21       22       22       22       23       23       24       26 
```

D'après ces résultats on voit que le résidus 6_A est impliqués dans le packing cristallin dans toutes les structures.

D'autres résidus sont retrouvés dans la majorité des structures (>80%) : `names(which(NbStructinPC/26 > 0.8))`


2. Représenter ces valeurs graphiquement



3. Calculer la moyenne et écart type de ce nombre





### Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin pour les PR1

1. Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin pour les PR1



```r
for (i in length(listAtomSyn)){
  y = 0
  for (j in length(v3)) {
    ma = matrice[j,i]
    if (ma == 1){
        y = y+1
    }
  ma = 0
  }
  matrice[length(v3)-2,i] = y
}

for (i in length(listAtomSyn)){
  y = 0
  for (j in length(v3)) {
    if (matrice[j,length(listAtomSyn)] == 1){
      if (matrice[j,i] == 1){
        y = y+1
      }
    }
    
  }
  matrice[length(v3)-1,i] = y
}
```


2. Représenter ces valeurs graphiquement
3. Calculer la moyenne et écart type de ce nombre



### Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin pour les PR2

1. Calculer le nombre de structure dans lequel un résidu est impliqué dans le packing cristallin pour les PR2



```r
for (i in length(listAtomSyn)){
  y = 0
  for (j in length(v3)) {
    if (matrice[j,listAtomSyn[160]] == 2){
      if (matrice[j,i] == 1){
        y = y+1
      }
    }
    
  }
  matrice[length(v3),i] = y
}

matrice
```

```
                   2_A 4_A 6_A 7_A 8_A 12_A 14_A 16_A 17_A 18_A 19_A 20_A 21_A 35_A 37_A 38_A 39_A 40_A 43_A 44_A 45_A 46_A 47_A 48_A 51_A 52_A 53_A 60_A 61_A 63_A 70_A 71_A 72_A 79_A 81_A 92_A 94_A 96_A 2_B 4_B 6_B 7_B 8_B 12_B 14_B 16_B 17_B 18_B 19_B 20_B 21_B 37_B 38_B 39_B 40_B 44_B 46_B 47_B
1hhp                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0    0    0    1    1    1    1    0    1    0    0    0    1    1    0    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1hih                 1   1   1   0   0    1    1    1    0    1    1    0    0    0    1    1    1    1    1    1    0    1    0    0    0    0    1    0    1    1    1    0    1    1    1    1    0    0   0   1   1   1   0    1    1    1    1    1    1    1    0    1    1    1    1    1    1    0
1hii                 1   1   1   1   0    0    1    1    1    1    1    0    0    0    0    0    1    1    1    0    1    1    0    0    0    0    1    1    1    0    1    0    0    0    0    1    1    0   1   1   1   1   0    0    0    0    0    0    0    0    0    1    1    1    1    1    1    0
1hiv                 1   1   1   0   0    1    1    1    0    1    1    0    1    0    1    1    1    0    0    1    1    1    0    0    0    1    1    0    1    1    1    1    1    1    1    1    0    0   0   1   1   1   0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    0
1hpv                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    1    0    0    0    1    1    0    1    1    1    1    1    0    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    0    1    1
1hsh                 1   0   1   1   1    1    0    1    1    0    0    0    0    1    1    0    0    1    1    1    1    1    0    0    1    0    0    1    1    1    0    0    1    1    1    0    0    0   0   0   1   1   1    0    0    1    1    1    1    0    0    0    1    1    0    1    1    0
1hsi                 1   1   1   1   0    0    0    1    1    1    0    0    0    0    0    0    1    1    0    0    0    1    0    0    0    0    1    0    0    1    0    0    1    0    0    0    0    0   1   0   0   1   1    1    0    1    1    1    1    0    1    1    1    1    0    1    1    0
1ivp                 1   1   1   1   0    1    1    0    0    0    1    1    0    0    1    0    1    0    0    1    1    1    0    0    0    0    1    1    0    0    0    0    0    1    1    1    1    0   0   0   1   1   0    0    1    0    1    0    1    0    0    0    0    1    1    1    0    0
1sdt                 0   1   1   1   0    0    0    0    0    1    0    0    0    1    1    0    0    0    1    1    1    1    0    0    0    0    1    0    0    0    0    0    1    0    0    1    1    0   1   0   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
2hb3                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    1    0    0    1    0    1    0    1    0    1    0    0    1    1    0   1   0   1   0   0    1    1    0    1    1    1    0    0    0    0    0    1    1    0    0
2hb4                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    1    1    1    0    1    0    0    0    1    1    1    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
2hpe                 1   1   1   1   0    0    0    0    0    0    1    0    0    1    1    0    1    1    0    1    0    0    0    0    1    1    0    1    1    0    0    0    0    1    1    0    0    0   1   0   1   1   0    1    0    0    1    1    1    1    1    1    0    1    1    0    1    0
2hpf                 0   0   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    1    1    0    1    1    0    0    0    0    0   0   0   1   1   0    1    1    1    1    0    1    1    1    1    0    1    0    0    1    0
2ien                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    0    0    1    1    0    1    0    1    0    1    0    0    1    1    0   0   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
2mip                 0   1   1   1   1    0    1    1    1    0    1    0    0    0    0    0    0    1    0    1    1    1    1    0    1    0    1    1    0    0    0    0    1    0    1    0    0    0   1   0   1   1   1    1    1    1    1    1    1    0    1    1    0    0    1    1    1    1
2nph                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    1    1    1    0    1    1    0    1    1    1    1    1    1    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0
2z4o                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    0    0    1    1    0    0    0    0    0    1    0    0    1    1    0   1   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
3ebz                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    1    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ec0                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    1    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    1    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ecg                 1   1   1   1   0    0    0    0    0    1    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    1    1    0    0    1    0    0    0    1   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
3ekv                 1   1   1   0   0    1    1    1    1    1    1    0    0    0    1    1    1    1    1    1    0    1    0    0    0    1    1    0    1    1    1    1    1    1    1    1    0    0   0   1   1   1   0    0    0    1    1    1    1    0    0    1    1    1    1    1    1    0
3nu3                 0   1   1   1   0    0    0    0    0    1    0    1    0    1    1    0    0    0    1    1    1    1    0    1    0    0    1    0    0    0    0    0    1    0    0    1    1    0   1   1   1   0   0    1    1    0    1    0    1    0    0    0    0    0    1    1    0    0
3phv                 1   1   1   1   1    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0    0    0    1    1    1    0    0    1    0    0    0    1    1    0    1    1    1   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3s45                 1   1   1   1   0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    0    1    0    0    1    0    0    0    1    0    0    1    0    0    0    0   1   0   0   1   0    1    1    0    1    1    1    1    1    1    0    0    1    1    1    0
4hla                 1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    0    1    0    0   1   1   1   1   0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
4ll3                 1   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   1   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
* NbStructinPC       0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR1     0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR2Nb   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                   48_B 52_B 53_B 61_B 63_B 70_B 71_B 72_B 79_B 92_B 96_B 41_A 35_B 43_B 51_B 57_B 42_A 54_A 55_A 57_A 58_A 65_A 98_A 99_A 10_B 41_B 42_B 45_B 55_B 81_B 94_B 98_B 99_B 1_A 3_A 5_A 9_A 24_A 25_A 26_A 27_A 29_A 49_A 50_A 56_A 67_A 69_A 87_A 90_A 91_A 93_A 95_A 97_A 34_A 36_A 73_A 88_A
1hhp                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
1hih                  0    1    1    0    1    1    1    1    1    1    1    0    1    1    0    1    1    0    1    0    0    0    0    0    0    1    1    1    1    1    1    0    0   0   0   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    0    0    0    1    0    0    0
1hii                  0    0    1    1    0    0    0    1    1    0    0    1    0    1    0    1    1    0    1    0    0    0    0    0    0    1    1    1    0    1    0    0    1   1   1   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    0    0    0    0    0
1hiv                  0    1    1    0    1    1    0    1    1    1    1    1    1    0    0    1    1    0    1    0    0    0    0    0    0    1    1    1    0    1    1    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    0    0
1hpv                  1    1    1    1    1    1    1    1    0    1    0    0    1    1    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1hsh                  0    1    1    1    1    1    0    0    1    0    0    1    0    0    1    1    1    0    1    1    0    0    0    1    0    0    1    1    1    0    0    0    0   0   0   0   0    0    0    0    0    1    0    0    0    0    0    1    0    0    0    0    0    1    0    0    0
1hsi                  0    0    1    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
1ivp                  0    0    1    1    0    0    0    1    1    0    0    1    0    0    0    1    1    0    1    0    0    1    0    0    0    0    1    0    0    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    1    1    0    0    0    0    0    0
1sdt                  0    0    0    1    1    0    1    1    1    1    0    0    0    1    0    0    0    0    0    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    1    1    1
2hb3                  0    0    0    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    1    0    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    0    0    0    0    0    0    1    0    0    0    1    1    1    1
2hb4                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
2hpe                  0    1    1    1    0    0    0    1    0    0    1    1    0    1    1    1    1    0    1    1    1    0    1    1    0    1    1    0    1    0    1    1    0   1   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0
2hpf                  0    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    1    0    0    0    0   0   0   1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
2ien                  0    0    0    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    1    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    1    0    0    0    0    0    1    0    0    0    1    1    1    1
2mip                  1    1    1    1    1    1    0    1    1    0    0    1    0    1    1    0    1    1    1    1    1    0    0    1    0    1    1    0    1    0    0    0    0   0   0   0   0    0    0    0    0    1    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0
2nph                  0    1    1    1    1    1    1    1    1    1    0    0    1    1    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0   0   1   0   0    0    0    0    0    0    1    0    0    1    0    0    0    0    0    0    0    0    0    0    0
2z4o                  0    0    0    1    1    1    1    1    1    1    0    0    0    1    0    0    0    0    0    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    1    0    0    0    0    0    1    0    0    0    1    1    1    1
3ebz                  1    0    1    1    0    0    0    0    1    0    1    1    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3ec0                  0    0    1    1    1    1    0    0    1    0    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0
3ecg                  1    0    1    0    0    1    0    0    1    0    1    1    0    1    0    0    1    1    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
3ekv                  1    1    1    0    1    1    1    1    1    1    1    0    1    0    0    1    1    0    1    0    0    1    0    0    0    0    1    1    0    1    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    1    0    0    0    1    0    0    0    1    1    0    0
3nu3                  0    0    0    1    1    1    1    1    1    1    0    0    0    1    0    0    0    0    1    1    1    0    0    0    0    1    1    0    1    1    0    0    0   0   0   0   0    0    0    0    0    0    1    0    0    0    0    0    0    1    0    0    0    1    1    1    1
3phv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0   1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0
3s45                  1    0    1    0    0    0    0    0    1    0    0    1    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
4hla                  1    1    1    1    1    1    1    1    1    1    0    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
4ll3                  1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
* NbStructinPC        0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR1      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
NbStructinPC.PR2Nb    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                   89_A 34_B 65_B 68_B 80_B 91_B 15_A 68_A 80_A 36_B 49_B 58_B 62_B 78_B 1_B 54_B 95_B 60_B 11_B 82_B 3_B 30_A 82_A 29_B 30_B 76_B 13_B 74_B 83_B 87_B 88_B 59_A 74_A 50_B 56_B 59_B 23_A 66_A 11_A 69_B 73_B 67_B 10_A 78_A Protease
1hhp                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0        0
1hih                  0    0    0    0    1    0    0    1    1    1    1    1    0    0   0    0    0    0    0    0   0    0    1    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0        0
1hii                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   1    0    0    1    0    0   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    1    0    0    1    0    0        0
1hiv                  0    0    0    0    1    0    0    1    1    0    1    0    0    0   0    0    0    0    0    0   0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    1    0        0
1hpv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
1hsh                  0    0    1    1    0    0    0    1    0    0    0    0    0    0   0    1    0    0    0    0   0    0    1    1    1    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    0    0    1    1        0
1hsi                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   1    0    0    0    0    0   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0        0
1ivp                  0    0    0    1    0    0    0    1    1    0    0    1    0    0   0    0    0    1    0    0   0    0    0    0    0    0    0    1    0    0    0    1    1    0    0    1    0    0    1    1    1    0    0    0        0
1sdt                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2hb3                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2hb4                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0        0
2hpe                  0    0    0    1    0    0    0    1    0    0    0    0    0    0   0    0    1    1    0    1   0    0    0    0    0    0    0    1    1    0    0    1    1    1    1    1    0    0    0    0    0    0    0    0        0
2hpf                  0    0    0    0    0    1    0    0    0    0    0    0    0    0   0    1    0    1    0    0   0    0    0    0    0    0    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2ien                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2mip                  0    0    1    1    0    0    0    0    0    0    0    0    0    0   1    1    0    0    0    0   1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2nph                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    1   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
2z4o                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    1   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ebz                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   1    1    0    1    1    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ec0                  0    0    0    1    0    0    0    0    0    0    0    0    0    0   1    1    0    1    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ecg                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   1    1    1    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3ekv                  0    0    0    0    1    0    1    1    1    1    1    1    1    1   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3nu3                  1    1    1    1    1    1    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3phv                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
3s45                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
4hla                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
4ll3                  0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
* NbStructinPC        0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
NbStructinPC.PR1      0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
NbStructinPC.PR2Nb    0    0    0    0    0    0    0    0    0    0    0    0    0    0   0    0    0    0    0    0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0        0
```

2. Représenter ces valeurs graphiquement



3. Calculer la moyenne et écart type de ce nombre



### Déterminer les résidus impliqués dans le packing cristallin dans toutes les structures


### Déterminer les résidus impliqués dans le packing cristallin dans toutes les structures de PR1


### Déterminer les résidus impliqués dans le packing cristallin dans toutes les structures de PR2


###Etudier le lien entre la conservation des résidus impliqués dans le packing cristallin et l'espace cristallo des structures

1. Déterminer l'espace cristallographique de chaque structure en allant sur la site de la PDB (rcsb.org)

