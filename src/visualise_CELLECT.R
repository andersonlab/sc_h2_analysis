suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(patchwork))

############################################################################
##              VARIABLES
############################################################################

if (interactive()){
    # Swap out variable names here if running interactively
    h2.partition.method <- 'LDSC'
    PROJECT_HOME <- '~/single-cell/sc_heritability_analysis'
    STUDY_DIR <- file.path(PROJECT_HOME,'/studies/gut-freeze003/ti-cd_healthy/')
    
    # Data variables
    cellect.results <- paste0(PROJECT_HOME, '/CELLECT/CELLECT-OUT/CELLECT-',h2.partition.method,'/results/prioritization.csv')
    # cellect.results <- paste0(PROJECT_HOME, '/CELLECT/CELLECT-OUT/CELLECT-',h2.partition.method,'/results/heritability.csv')

  
    # dataset.list <- c("freeze003-ti-cd_healthy.Cell-type"="SEG", "freeze03-ti-cd_healthy.mast.Cell-type.zscore"="DEG")
    dataset.list <- c("freeze003-ti-cd_healthy.NMF.k75.gep"="PEG")
    # dataset.list <- c("freeze003-blood-cd_healthy.Azimuth"="SEG")
    
    gwas.list <- c(paste0(c('CD'), '_DeLange2017'),"Height_Yengo2018", "EA3_Lee2018") # "CD_DeLange2017"  "UC_DeLange2017"  "IBD_DeLange2017" "Height_Yengo2018" "EA3_Lee2018"
    
    heatmap.type <- 'GWAS' # GWAS or DATASET or COND
    
    # heatmap.show <- c("freeze003-ti-cd_healthy.Cell-type","freeze003-ti-cd_healthy.Cell-type_Healthy","freeze003-ti-cd_healthy.Cell-type_CD") # "freeze003-ti-cd_healthy.Cell-type" "freeze003-ti-cd_healthy.Cell-type_Healthy" "freeze003-ti-cd_healthy.Cell-type_CD"
    
    
} else {
    # Reads in flags when ran as a script
    option_list = list(make_option(c("-r", "--results"), type="character", default=NULL, help="CELLECT results file", metavar="character"),
                       make_option(c("-d", "--dataset"), type="character", default=NULL, help="Name of datasets provided to CELLECT. N.b. If providing multiple these must have the same annotations", metavar="character"),
                       make_option(c("-g", "--gwas"), type="character", default=NULL, help="Names of GWAS provided to CELLECT", metavar="character"),
                       make_option(c("-h", "--heritability"), type="character", default='LDSC', help="Use LDSC or MAGMA results. Default is LDSC", metavar="character"),
                       make_option(c("-p", "--plot"), type="character", default='FOREST', help="Main plot is FOREST plot of betas and s.e. or SCATTER plot of p-values. Default is FOREST", metavar="character"),
                       make_option(c("-hm", "--heatmap"), type="character", default='GWAS', help="Heatmap will show other GWAS or DATASETs. Default is GWAS", metavar="character"),
                       make_option(c("-hms", "--heatmapshow"), type="character", default=NULL, help="The values to plot in the heatmap, dependent on the --heatmap flag. Default is to plot everything", metavar="character"))
    
    opt_parser = OptionParser(option_list=option_list)
    opt = parse_args(opt_parser)
    
    PROJECT_HOME = Sys.getenv('PROJECT_HOME')
    STUDY_DIR = Sys.getenv('STUDY_DIR')
    
    if (!is.na(opt$results)) {
        cellect.results <- trimws(opt$results) # cellect.results <- paste0(PROJECT_HOME, '/CELLECT/CELLECT-OUT/CELLECT-',h2.partition.method,'/results/prioritization.csv')
    } else {
        stop("CELLECT results must be provided")
    }
    if (!is.na(opt$dataset)) {
        dataset <- trimws(opt$dataset) # dataset <- "freeze003-ti-cd_healthy.Cell-type"
    } else {
        stop("CELLEX dataset name must be provided")
    }
    if (!is.na(opt$gwas)) {
        main.gwas <- strsplit(trimws(opt$gwas), ',') # main.gwas <- "CD_DeLange2017,EA3_Lee2018" becomes c("CD_DeLange2017EA3_Lee2018")
    } else {
        stop("GWAS sum stats name must be provided")
    }
    if (dir.exists(PROJECT_HOME) & dir.exists(STUDY_DIR)) {
    } else {
        stop("PROJECT_HOME and STUDY_DIR directories must exist")
    }
    if (opt$heritability %in% c('LDSC', 'MAGMA')) {
        h2.partition.method <- trimws(opt$heritability) 
    } else {
        stop("Heritability partioning method must be either LDSC or MAGMA")
    }
    if (opt$plot %in% c('FOREST', 'SCATTER')) {
      plot.type <- trimws(opt$plot) 
    } else {
        stop("Main plot type must be either FOREST or SCATTER")
    }
    if (opt$heatmap %in% c('GWAS', 'DATASET')) {
      heatmap.type <- trimws(opt$heatmap.type) 
    } else {
        stop("The heatmap option provided must be either GWAS or DATASET")
    }
    if (!is.na(opt$gwas)) {
      main.gwas <- trimws(opt$gwas) # main.gwas <- "CD_DeLange2017"
    } else {
        stop("GWAS sum stats name must be provided")
    }
}

out.dir <- file.path(STUDY_DIR,'plots/')
relabeller.file <- file.path(STUDY_DIR, 'data/name_mappings.txt')

############################################################################
##              CODE
############################################################################

# Import some common functions for plotting
source(file.path(PROJECT_HOME,'src/sc_heritability_lib.R'))

# Read in dataframe containing name remappings then make named vectors from the
# key-value column pairings
remap.labels <- read_tsv(relabeller.file)
gwas.labeller <- deframe(drop_na(remap.labels[paste0('GWAS', c('', '.labels'))]))
annot.labeller <- c()
for (ds in names(dataset.list)){
  # Make new temporary named vector
  tmp.annot.labeller <- deframe(drop_na(remap.labels[paste0(ds, c('', '.labels'))]))
  
  # Drop labels from the temp named vector that are already in the existing annot.labeller
  tmp.annot.labeller <- tmp.annot.labeller[!(names(tmp.annot.labeller) %in% names(annot.labeller))]
  annot.labeller <- c(annot.labeller, tmp.annot.labeller)
}

# Read data and mutate labels to be more readable
cellect.results.df <- read_csv(cellect.results) %>%
  mutate(gwas.label = case_when(gwas %in% names(gwas.labeller) ~ unname(gwas.labeller[gwas]),
                                TRUE ~ gwas),
         # annot.label = annotation)
         annot.label = case_when(annotation %in% names(annot.labeller) ~ unname(annot.labeller[annotation]),
                                TRUE ~ annotation),
         SEG.DEG.PEG = dataset.list[specificity_id])


  # Filter results for those matching the GWAS and dataset provided
filtered.cellect.results.df <- cellect.results.df %>%
  filter(gwas %in% gwas.list, specificity_id %in% names(dataset.list))

# Bonferonni alpha corrected by number of cell-types
bonf.alpha <- 0.05 / length(unique(filtered.cellect.results.df$annot.label))

# Add annotation groupings (if present)
annotation.groupings <- read.csv('~/single-cell/ti_freeze003-annot_grouping.csv')

filtered.cellect.results.df <- filtered.cellect.results.df %>% 
  left_join(annotation.groupings)


# Make label columns into factor for ordering
filtered.cellect.results.df$gwas.label <- factor(filtered.cellect.results.df$gwas.label,
                                        levels=unique(c(gwas.labeller,
                                                        (filtered.cellect.results.df$gwas))))
filtered.cellect.results.df$annot.label <- factor(filtered.cellect.results.df$annot.label,
                                         levels=rev(unique(c(annot.labeller,
                                                             (filtered.cellect.results.df$annotation)))))
filtered.cellect.results.df$SEG.DEG.PEG <- factor(filtered.cellect.results.df$SEG.DEG.PEG,
                                                  levels=c("SEG", "DEG", "PEG"))

# Make the heatmap
p <-heatmapCELLECTpval(filtered.cellect.results.df, alpha, heatmap.type)
# Custom labels for heatmap
p <-   p + facet_grid(grouping.label ~ SEG.DEG.PEG, scales = 'free', space = 'free') 
p

# Decide plot size based on number of rows in df
if (nrow(filtered.cellect.results.df) > 25) {
  plot.width <- 14
  plot.height <- 11
}
if (nrow(filtered.cellect.results.df) <= 25) {
  plot.width <- 12
  plot.height <- 5
}


# If labels are short, make plot smaller
longest.label <- max(unlist(lapply(filtered.cellect.results.df$annot.label, function(x) nchar(as.character(x)))))

# Simultaneous equation solving for a * longest.label + b = 2.8, 6
plot.width <-  (0.08 * longest.label) + 2.32
# heatmap.or.main <- paste(heatmap.or.main, heatmap.type,sep = '-')

# 
# for (ext in c('pdf', 'png')){
#         plot.path <- paste0(out.dir, ext,'/',
#                             paste0(dataset.list, collapse = '-'), gwas[1],
#                             "heatmap_LDSC",
#                             "_visualiseCELLECT.", ext)
#         print(plot.path)
#         ggsave(plot.path, width = plot.width, height = plot.height)
# }

########################
# Build super plot
########################


# p.seg <- p.final + theme(axis.text.x = element_blank(), legend.position = "none")
# p.deg <- p.final + theme(axis.text.x = element_blank(), legend.position = "none")
# p.peg <- p.final + guides(fill = guide_legend(override.aes = list(fill = "grey24",
#                                                                 alpha = seq(.1, 0.7, length.out = 6))))
# 
# p.all <- p.dendro + p.seg + p.deg + p.peg + plot_layout(widths = c(8,1,1,1))


