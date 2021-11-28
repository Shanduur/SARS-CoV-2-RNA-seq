# Master Thesis

*Analysis of smoking status impact on molecular mechanisms of SARS-CoV-2 viral entry through single-cell sequencing experiments.*

Repository containg Master's Thesis source and data.

# Dataset

## Source

The data that is located in `Data/` directory, is downloaded from [singlecell.broadinstitute.org][sc-broadinstitute]:
- [HCA LungMAP COVID-19 Smokers Lung (277k+ cells, genes)][data-smokers]
- [HCA LungMAP COVID-19 Internal Non-Smokers Lung (96k+ cells, 27 genes)][data-nonsmokers]

## How dataset was prepared?

Based on the notebooks and publication, we can deduce that:
1. Multiple dataset were gathered;
2. Divided to **smoking** if `stat` was in [`current`, `smoked`, `active`, `former`, `heavy`, `light`], **nonsmoking** and **NaN**;
3. Dropped **NaN**;

## On what genes the analysis focused?

The surface receptor angiotensin-converting enzyme 2 (ACE2) and the associated proteases, transmembrane protease serine 2 (TMPRSS2) and Cathepsin L (CTSL), were previously identified mediators of SARS-CoV cellular entry. In the [original publication][pub-org] single-cell RNA-seq (scRNA-seq) across diverse tissues to assess the cell-type-specific expression of ACE2, TMPRSS2, and CTSL. Specific subsets of respiratory epithelial cells were identified as putative targets of viral infection, including subsets of epithelial cells in the nasal passages, lung and airways. Additionally, they detected expression in other tissues that may serve as routes of viral transmission, including the gut and corneal epithelia, and in cells potentially associated with COVID-19 clinical pathology including cardiomyocytes, olfactory sustentacular cells, and renal epithelial cells.

# Resources for COVID Smokers data

## Videos

- [StatQuest: A gentle introduction to RNA-seq][yt-statquest-rna-seq]
- [Single Cell Sequencing - Eric Chow (UCSF)][yt-chow-sc-seq]
- [Single cell RNA sequencing Playlist][yt-rna-seq-lst]
- [Stephanie Hicks: Scalable statistical methods and software for single-cell data science][yt-hicks-sc-seq]

## Publications

- [Original publication][pub-org]
- [A practical guide to single cell RNA sequencing][pub-sc-rna]

## Courses

- [Introduction to single cell RNA seq][intro-sc-rna-seq]


<!-- Resources -->

[data-smokers]: https://singlecell.broadinstitute.org/single_cell/study/SCP876/hca-lungmap-covid-19-smokers-lung?scpbr=hca-covid-19-integrated-analysis
[data-nonsmokers]: https://singlecell.broadinstitute.org/single_cell/study/SCP875/hca-lungmap-covid-19-internal-nonsmokers-lung?scpbr=hca-covid-19-integrated-analysis
[intro-sc-rna-seq]: https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html
[sc-broadinstitute]: https://singlecell.broadinstitute.org/
[pub-org]: https://www.biorxiv.org/content/10.1101/2020.04.19.049254v2
[pub-sc-rna]: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4
[yt-statquest-rna-seq]: https://www.youtube.com/watch?v=tlf6wYJrwKY
[yt-chow-sc-seq]: https://www.youtube.com/watch?v=k9VFNLLQP8c
[yt-rna-seq-lst]: https://www.youtube.com/playlist?list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN
[yt-hicks-sc-seq]: https://www.youtube.com/watch?v=Sqr2UFpJKkM