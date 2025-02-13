 # Projet : Analyse du Séquençage NGS du SARS-CoV-2

## Description

Ce projet a pour objectif d'analyser et d'interpréter des données de séquençage haut débit (NGS) du SARS-CoV-2 à travers différentes étapes, allant du prétraitement des reads à l'identification des variants.

L'analyse comprend :

- La préparation des données (récupération et nettoyage des reads).
- L'assemblage de novo du génome du virus.
- Le mapping des reads contre une séquence de référence.
- L'identification des variants présents dans les séquences obtenues.

## Prérequis

Avant d’exécuter les scripts, il faut s'assurer que les outils suivants sont installés :

Python 3.8+
FastQC (évaluation de qualité des reads)
Trimmomatic (nettoyage des reads)
SPAdes (assemblage de novo)
QUAST (évaluation de l’assemblage)
BWA (alignement des reads)
Samtools (traitement des fichiers BAM/SAM)
GATK (détection des variants)
VCFtools (analyse des variants)
SnpEff (annotation des variants)

## Utilisation


### Auteur
Dinan BRISSON, ISEN-TMS
