# sc_h2_analysis


## General

For the purpose of this analysis, the instructions described in the [CELLECT-LDSC](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial) tutorial were followed. CELLECT aims at the quantification of the association between GWAS signal (heritability) and cell-type expression specificity of genes. While the cell-type expression can be generated in multiple ways, here we use [CELLEX](https://github.com/perslab/CELLEX) and [MAST](https://github.com/RGLab/MAST). We use established genetic prioritization models kindly implemented by the CELLECT authors to perform genetic prioritization (S-LDSC, Finucane et al., 2015, Nature Genetics; De Leeuw et al., 2015, PLoS). The result of CELLECT is a list of prioritized etiologic cell-types for a given human complex disease or trait.

## GWAS

Next, we use the IBD, CD and UC summary statistics from De Lange et al. which we have provided in this repository, as well as Educational Attainment (Lee et al.) and Height (Yengo et al.) GWAS sumstats (as negative controls). Precuse instructions for how to download the EA other GWAS can be found in the [CELLECT GitHub repo](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial#step-1-download-and-munge-gwas). The script `mtag_munge.py` can be found in this repo `sc_heritability_analysis/CELLECT/ldsc/mtag_munge.py`.

* **Ulcerative colitis**
  * [GCST004133](https://www.ebi.ac.uk/gwas/studies/GCST004133)
  * Sample size: N = 45975
* **Crohn's Disease**
  * [GCST004132](https://www.ebi.ac.uk/gwas/studies/GCST004132)
  * Sample size: N = 40266
* **Inflammatory Bowel Disease**
  * [GCST004131](https://www.ebi.ac.uk/gwas/studies/GCST004131)
  * Sample size: N = 59927
* **Height**
  * [GWAS Download package](https://github.com/mikegloudemans/gwas-download)
* **Educational attainment**
  * [GWAS Download package](https://github.com/mikegloudemans/gwas-download)

## Sources

* [GWAS Summary Statistics](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
* [CELLECT](https://github.com/perslab/CELLECT)
* [Nextflow](https://www.nextflow.io/)
* [CELLEX](https://github.com/perslab/CELLEX)
* [LDSC](https://github.com/bulik/ldsc)
* [MAGMA](https://ctg.cncr.nl/software/magma)

***

# Quick start
If you aren't fussed by the details and just wish to perform the analysis see below:

Install dependencies 

```bash
conda install -c conda-forge git-lfs
# Install repo somewhere with minimum 15GB storage spare
git clone --recurse-submodules https://github.com/andersonlab/sc_heritability_analysis.git 
conda env create -f sc_heritability_env.yml

conda activate sc_heritability
```

Set up environment

```bash
export PROJECT_HOME=~/single-cell/sc_heritability_analysis
export DATASET=freeze03-ti-cd_healthy
export STUDY_DIR=$PROJECT_HOME/studies/gut-freeze003/ti-cd_healthy

mkdir -p $STUDY_DIR/{plots/{pdf,png},results,data}
```


Run CELLEX
```bash
cd $PROJECT_HOME
./src/run_CELLEX.py \
--h5_anndata Discovery.h5ad \
--output_file $STUDY_DIR/results/$DATASET \
--annotation_columns Cell-type \
--verbose True
```


Run CELLECT
```bash
cd $PROJECT_HOME/CELLECT
# Run these one at a time
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile $PROJECT_HOME/cellect_config.yml
snakemake --use-conda -j -s cellect-magma.snakefile --configfile $PROJECT_HOME/cellect_config.yml
snakemake --use-conda -j -s cellect-genes.snakefile --configfile $PROJECT_HOME/cellect_config.yml
```

Your results should be in your ```STUDY_DIR```. Feel free to plot them yourself or see how we chose to plot them below.

***

## References

* Timshel et al., 2020, eLife <https://elifesciences.org/articles/55851>
* Bulik-Sullivan et al., 2015, Nature Genetics <https://www.nature.com/articles/ng.3211>
* Finucane et al., 2015, Nature Genetics <https://www.nature.com/articles/ng.3404>
* De Leeuw et al., 2015, PLOS <https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004219>
* Finucane et al., 2018, Nature Genetics <https://www.nature.com/articles/s41588-018-0081-4>
* De Lange et al., 2017, Nature Genetics <https://www.nature.com/articles/ng.3760>
* Finak et al., 2015, Genome Biology <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5>
* Skene et al., 2018, Nature Genetics <https://www.nature.com/articles/s41588-018-0129-5>
* Corridoni et al., 2020, Nature Medicine <https://www.nature.com/articles/s41591-020-1003-4>
* Martin et al., 2019, Cell <https://www.sciencedirect.com/science/article/pii/S0092867419308967>

***

## Authors

Tobi Alegbe, and Moritz Przybilla