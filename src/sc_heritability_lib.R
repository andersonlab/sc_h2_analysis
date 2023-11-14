scatterPlotCELLECTpval <- function(tibble, alpha, point.colour){
  # Scatter plot of p values against cell-type labels with bonferonni significant 
  # cell-types coloured and labelled
  
  p <- ggplot() + 
    geom_point(data=tibble, aes(y=annot.label, x=-log10(pvalue)), color="gray") +
    geom_point(data=filter(tibble, pvalue <= alpha),
               aes(y=annot.label, x=-log10(pvalue)), color=point.colour) +
    geom_vline(xintercept=-log10(bonf.alpha), linetype="dashed", color="darkgray") + 
    ggrepel::geom_text_repel(data=filter(tibble, pvalue <= alpha),
                             aes(y=annot.label, x=-log10(pvalue), label=annot.label), hjust = 0, nudge_y = 1.5, show.legend=F) + 
    theme_classic() + 
    theme(axis.text.x=element_text(size=rel(1.5))) + # REF: https://stackoverflow.com/a/47144823/6639640
    theme(legend.position="bottom")
  return(p)
  
}

forestPlotCELLECTbeta <- function(tibble, alpha, point.colour){
  # Forest plot of betas/weights against cell-type labels with bonferonni significant 
  # cell-types coloured and labelled
  p <- ggplot() + 
    geom_pointrange(tibble, mapping=aes(y=annot.label, x=beta, xmin=beta-beta_se, xmax=beta+beta_se), color= 'gray') +
    geom_point(data=dplyr::filter(tibble, pvalue <= alpha),
               mapping=aes(y=annot.label, x=beta), color=point.colour) +
    geom_vline(xintercept=0, linetype="dashed", color="darkgray") +
    ggrepel::geom_text_repel(data=dplyr::filter(tibble, pvalue <= alpha),
                             mapping=aes(y=annot.label, x=beta, label=annot.label), hjust = 0, nudge_y = 1.5, show.legend=F) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(1.5))) + # REF: https://stackoverflow.com/a/47144823/6639640
    theme(legend.position="bottom")
    
  return(p)
}

heatmapCELLECTpval <- function(tibble, alpha, gwas.or.SDP.or.cond){
  # Heatmap of p values with darker colours indicating a lower p value and an 
  # asterisk to show bonferonni significance has been met
  # gwas.or.categ.or.cond refers to where the data being plotted is either
  # gwas - different gwas sum stats
  # SDP - different SEG/DEG/PEG inputs
  # cond - conditioning on SDP categories i.e. cell-types or factors
  
  # Change the plot depending on the gwas.or.SDP.or.cond parameter
  x.axis.vals <- "gwas.label"
  x.axis.title <- "GWAS"
  if (gwas.or.SDP.or.cond == "SDP"){
    x.axis.vals <- "specificity_id"
    x.axis.title <- "Expression dataset"
  }
  if (gwas.or.SDP.or.cond == "COND"){
    x.axis.vals <- "conditional_annotation"
    x.axis.title <- "Conditional annotation"
  }
  
  # Make the plot
  p <- ggplot(tibble , aes_string(x=x.axis.vals, y="annot.label")) + 
    geom_tile(aes(alpha=-log10(pvalue),fill=SEG.DEG.PEG)) + 
    geom_text(data=filter(tibble, pvalue <= bonf.alpha),
              label="*", color='white', hjust=0.5, vjust=0.75, size=8) + # add asterisk if fdr significant
    theme_minimal() + # set this first
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), # remove grid
          legend.position="right",legend.title = element_text(size=7), # position legend
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 90, vjust = 0.8,hjust = 0.95),
          strip.text.y = element_blank(),
          strip.text.x = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 18)) +
    labs(x=x.axis.title, y='',alpha=expression(-log[10](P[CELLECT-LDSC])),
         fill="Single-cell\ninput") +
    scale_fill_manual(values=c("DEG" = "darkblue",
                               "PEG" = "goldenrod",
                               "SEG" = "darkred"),
                      guide='none') +
    scale_alpha_continuous(range = c(0.1, 1),limits=c(0,5))
  
  return(p)
}

findSigAnnots <- function(tibble, alpha, in.datasets){
  # Finds which annotations in a given tibble are significant after bonferonni correction
  # and returns tibble with logical column indicating as suhc
  filtered.tibble <- tibble %>%
    filter(specificity_id %in% in.datasets) %>% 
    group_by(specificity_id, gwas) %>% 
    mutate(bonf.sig = pvalue <= 0.05/n()) %>% 
    ungroup()
  
  return(filtered.tibble)
}
