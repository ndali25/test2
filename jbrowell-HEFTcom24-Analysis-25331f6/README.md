HEFTcom24-Analysis

Analyse de la compétition Hybrid Renewable Energy Forecasting and Tracing Competition 2024, reproduisant et étendant les résultats présentés dans cet article
.

1. Récupération des données

Les fichiers nécessaires doivent être téléchargés depuis Zenodo et placés dans le dossier data/ du projet.

Fichiers à télécharger depuis

https://doi.org/10.5281/zenodo.13950764

trades.csv

pinball.csv

forecasts.csv

Energy_Data_20200920_20240118.csv

Energy_Data_20240119_20240519.csv

overall_leaderboard.csv

HEFTcom Reports.csv

2. Exécution du code R

Le script principal analysis.R utilise ces fichiers pour réaliser l’analyse.

Pré-requis

R (version 4.2.1 recommandée)

RStudio (optionnel mais recommandé)

Packages R suivants (assurez-vous qu’ils sont installés) :
patchwork, latex2exp, xtable, ggridges, ggplot2, rstudioapi, data.table, dplyr

Exemple pour installer les packages (dans R) :
install.packages(c("patchwork", "latex2exp", "xtable", "ggridges", "ggplot2", "rstudioapi", "data.table", "dplyr"))

Lancer l’analyse :

Ouvrez RStudio ou R, puis exécutez :

source("analysis.R")

3. Informations supplémentaires
Citation du dépôt

Ce dépôt est archivé sur Zenodo : https://doi.org/10.5281/zenodo.14247209

DOI général : 10.5281/zenodo.14247209

@misc{Browell2024HEFTcomAnalysis,
    title = {{jbrowell/HEFTcom24-Analysis}},
    year = {2024},
    author = {Browell, Jethro},
    publisher = {Zenodo},
    doi = {10.5281/zenodo.14247209}
}

Citation de l’article HEFTcom
@misc{browell2025hybridrenewableenergyforecasting,
      title={The Hybrid Renewable Energy Forecasting and Trading Competition 2024}, 
      author={Jethro Browell and Dennis van der Meer and Henrik Kälvegren and Sebastian Haglund and Edoardo Simioni and Ricardo J. Bessa and Yi Wang},
      year={2025},
      eprint={2507.01579},
      archivePrefix={arXiv},
      primaryClass={stat.AP},
      url={https://arxiv.org/abs/2507.01579} 
}

4. Environnement R utilisé

L’analyse a été réalisée avec :

R version 4.2.1 (Windows 10 x64)

RStudio 2023.09.1+494 "Desert Sunflower"

Packages et versions :

> sessionInfo()
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 26100)

attached base packages:
[1] stats graphics grDevices utils datasets methods base

other attached packages:
[1] patchwork_1.2.0 latex2exp_0.9.6 xtable_1.8-4 ggridges_0.5.4 ggplot2_3.5.1 rstudioapi_
