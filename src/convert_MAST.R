################################################################################################################################################
##                                             
##  Converts MAST results to CELLECT DEG input
##                                       
##  Date: 22 Sep 2023             
##  
##  Authors: Tobi Alegbe                               
##                                            
##                                              
################################################################################################################################################

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))


############################################################################
##  VARIABLES
############################################################################

option_list = list(make_option(c("-o", "--outpath"), type="character", default=NULL, help="Complete path of where to save output", metavar="character"),
                    make_option(c("-d", "--dge_file"), type="character", default=NULL, help="Complete path of the DGE file provided as input", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
    

out.path <- trimws(opt$outpath)
out.dir <- '/nfs/users/nfs_o/oa3/single-cell/code/sc_heritability_analysis/studies/gut-freeze003/replication_labeled/results/'
dge.file.name <- trimws(opt$dge_file)

############################################################################
##  CODE
############################################################################

# Reading in the dataframe returned by https://github.com/andersonlab/sc_nf_diffexpression
print('Reading in data...')
dge.df <- read_tsv(dge.file.name)
print('Data successfully read.')

# For rows that contain a test stat, min-max normalising the test stat
print('Transforming test statistic...')
transformed.dge.df <- dge.df %>% 
    filter(!is.na(test_statistic)) %>%
    group_by(cell_label) %>%
    mutate(normalised_stat = (abs(test_statistic) - min(abs(test_statistic))) / (max(abs(test_statistic)) - min(abs(test_statistic))),
           normalised_lfc = (abs(log2fc) - min(abs(log2fc))) / (max(abs(log2fc)) - min(abs(log2fc)))) %>%
    ungroup() %>% 
    mutate(cell_label = str_replace_all(cell_label, ' ', '_')) %>%
    mutate(cell_label = str_replace_all(cell_label, '[^0-9a-zA-Z+\\-\\(\\)_]+', ''))
print('Data successfully transformed.')

for (stat in c("normalised_stat", "normalised_lfc")){
    p <- ggplot(transformed.dge.df, aes_string(x=stat)) + 
        geom_histogram() +
        theme_classic() +
        facet_wrap(~cell_label)
    print(p)
}

print('Getting data into correct format for CELLECT...')
cellect.stat.dge.df <- transformed.dge.df %>%
    select(gene, cell_label, normalised_stat) %>%
    spread(key=cell_label, value=normalised_stat,fill = 0)
cellect.lfc.dge.df <- transformed.dge.df %>%
    select(gene, cell_label, normalised_lfc) %>%
    spread(key=cell_label, value=normalised_lfc,fill = 0)
print('Data successfully transposed.')
    

out.stat.name <- paste0(out.name,'-norm_test_stat.deg.csv.gz')
out.lfc.name <- paste0(out.name,'-norm_lfc.deg.csv.gz')
write_csv(cellect.stat.dge.df,file = paste0(out.path,'-norm_test_stat.deg.csv.gz'))
write_csv(cellect.lfc.dge.df,file = paste0(out.path,'-norm_lfc.deg.csv.gz'))
