---
title: "Stage de Maxime Kermarrec (L3-BI)"
author: "Leslie REGAD et Maxime KERMARREC"
date: '2019-02-27'
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


# Projet 1 : Etude du packing cristallin chez PR1 et PR2
du Vendredi 18/01/2018 au ..  

[Résultat](projet1.html)


## Planning pour le vendredi 11/02/2019 
A finir : 

* les représentations des positions impliquées dans le packing sur les structures 3D  
* Extraire les informations sur le packing cristallin et sur la résolution
* Classer les protéines suivant leur ligand. Pour cela, réaliser une classification hiérarchique des ligands [code](projet1.html#classement_des_protéines_suivant_leur_ligand)




Il faudrait que ces parties soient finies pour la semaine prochaine.
Que je puisse voir les résultats, car je vais devoir envoyer le papier en fin de semaine prochaine.



# Projet 2 : Comparaison de la flexibilité des atomes des ligands chez PR1 et PR2

début du projet : vendredi 01/03/2019  


[Résultat](projet2.html)


## Planning pour le vendredi 22/02/2019 

* A l'aide du site www.rcsb.org, déterminer pour chaque protéine :    
    + son groupe cristallin. Exemple pour 1hsi:  	P 1 21 1  
    + sa résolution. Exemple pour 1hsi: 2.5  
* Déterminer les groupes de la classification hiérarchique obtenue à partir du packing cristallin. Pour cela, à l'aide de la fonction  `cutree` couper l'arbre crée en 5 groupes.  
* Pour chaque groupe,   
    + Est-ce que les protéines appartenant aux mêmes groupes appartiennent au même groupe cristallin ?   
    + Est-ce que les protéines appartenant aux mêmes groupes ont une résolution proche    
    + Est-ce que les protéines appartenant aux mêmes groupes sont complexés au même ligand    
Pour répondre à ces questions, il faudrait faire la figure 5 du papier [Triki et al., 2018](biblio/Triki_scReports_2018_asym.pdf).  
Avec un logiciel de gestion de figure gimp, inkscape power point ou openoffice, il faudrait réalisé une figure comme la figure 5 où :    
   + la 1ere figure correspond à la classification hiérarchique   
   + la 1ere ligne de cases de couleur représente le groupe cristallin  
   + la 2ème ligne de cases de couleur représente la résolution  
   + la 3ème ligne de cases de couleur représente la classe de ligand (cf. la classification des ligands que l'on a réalisé)  
   + la 4ème ligne de cases de couleur représente PR1/PR2  
* Pour l'évaluation de votre stage, vous aurez une soutenance à faire. Commencer à résumer ce travail fait sur le packing cristallin dans une présentation. Pensez que les membres du jury seront mixtes, informaticiens, biologistes et bioinformaticiens, statisticiens.
    




## Planning pour le vendredi 01/03/2019 

A partir d'aujourd'hui, on change de projet.
Cette fois-ci vous allez travailler sur la comparaison de la flexibilité des atomes de ligands.

* Etape 1 : Classer les atomes de ligands des PR1  
* Etape 2 : Classer les atomes de ligands des PR2  

