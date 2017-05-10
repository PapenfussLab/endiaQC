## Using phyloseq data structures to make abundance charts for Phylum, Order or Genus
## Bar charts and box-and-whiskers plots, for publication. Figures 2A, 3A, 4A and S3

## Working at stool-sample level - SampleNumber = 72 samples
## and sometimes at aliquot or biological replicate level = 215 or 216 samples
## Use OTUs with counts >= 20, as proportions of samples
## June 2016

library("phyloseq")
library(ggplot2)
library(reshape2)
library("RColorBrewer")
theme_set(theme_bw(base_size=8))

Ntop <- 10
palette_phyl <- c(brewer.pal(9, "Set1"), "#000000", brewer.pal(9, "Pastel1"))  
## extended for more than 9 phyla. 
## Swap some colours: Red, blue, green, brown, yellow, purple, orange
palette_phyl <- palette_phyl[c(1:3, 7, 6, 4, 5, 8:length(palette_phyl))]
pie(rep(1, length(palette_phyl)), col=palette_phyl)
## Consistent colours needed for 2 data sets, so hard-coding phylOrder. As Baylor seqs gave 12 phyla,
## I am using that order, rather than the local WEHI seqs 10 phyla order
## Top 7 or 9 the WEHI and Baylor sets are equal - top 8 they differ
phylOrder <- c( "Bacteroidetes",  
               "Firmicutes",     
               "Proteobacteria", 
               "Actinobacteria", 
               "Verrucomicrobia",
               "Tenericutes",    
               "Cyanobacteria",  
               "Euryarchaeota",  
               "Lentisphaerae",  
               "Fusobacteria",   
               "TM7",            
               "Deferribacteres")

## For same reason, hardcode orderOrder
## I want a mix of highest average proportion and highest maximum proportion
## The 12 that have a proportion in some sample greater than 2% in Baylor data makes 
## a good list - it includes the WEHI top 10 and Baylor top 9. 
## Union of the 2 top 10s adds Desulfovibrionales, a 5th Proteobacteria

# tax_table(psL4prop)[which(apply(otu_table(psL4prop), 2, max) >=0.02), c(2,4)]
ordOrder <- c( # 1 x bacteriodetes, 2 x firmicutes
  "Bacteroidales", "Clostridiales", "Erysipelotrichales", 
  # 4 x proteobacteria
  "Burkholderiales", "Pasteurellales", "RF32", "Enterobacteriales", 
  # 1 actinobacteria 
  "Bifidobacteriales", 
  # 1 verrucomicrobia 
  "Verrucomicrobiales", 
  # 2 x tenericutes 
  "Anaeroplasmatales", "RF39", 
  #1 cyanobacteria
  "YS2" )
palette_order <- c(palette_phyl[1], #  red for bacter,
                   brewer.pal(3, "Blues")[3:2],   # 2 x blue for firmi
                   brewer.pal(5, "Greens")[5:2],   # 4 shades of green, dark to pale
                   palette_phyl[4:5],  # brown for actino, yellow for verruco
                   brewer.pal(3, "Purples")[3:2], # 2 x purple for tenericutes
                   palette_phyl[7]     # orange for cyanobacteria
)

## List of genus, or higher label where genus undefined,  present in at least 10% in any sample. 
##  Compromise list for showing representative genera in WEHI and Baylor data
genusLabelList <- c( # 4 o__Bacteroidales
  'g__Bacteroides', 'g__Prevotella' , 'f__S24-7', 'o__Bacteroidales'
  # 5 o__Clostridiales 
  , 'g__Lachnospira', 'g__Blautia', 'f__Lachnospiraceae', 'g__Faecalibacterium', 'o__Clostridiales' 
  # 3 p__Proteobacteria 
  , 'g__Sutterella', 'g__Haemophilus', 'Proteobacteria' 
  # 1 p__Actinobacteria
  , 'g__Bifidobacterium'
  # 1 p__Verrucomicrobia
  , 'g__Akkermansia'
  )
palette_genus <- c(brewer.pal(5, "Reds")[c(4,2,3,5)], #  4 x red for bacter,
                   brewer.pal(6, "Blues")[c(4,2,5,3,6)],   # 5 x blue for firmi
                   brewer.pal(4, "Greens")[c(3,2,4)],   # 3 shades of green for proteo
                   palette_phyl[4:5],  # brown for actino, yellow for verruco
                   palette_phyl[10]  # black for 'other'
)


applyNewGenusLabel <- function(tax){   ## Match the lowest possible taxa name to a name in genusLabelList
  tax_char <- as.vector(tax, mode="character")
  checkrank <- 7
  newRank <- tax_char[checkrank]    # 'lastDefRank' in column 7 is starting point 
  while (! newRank %in% genusLabelList & checkrank >2){ # if not in list, return phylum=Rank1
    checkrank <- checkrank - 1
    newRank <- tax_char[checkrank] 
  }
  newRank
}

genusBarplot <- function(phyl_L6props,
                         phylOrd = phylOrder,
                         Nphyla = 6, 
                         palGenus = palette_genus, 
                         genusOrder = genusLabelList
){
  tax_table(phyl_L6props)[,"Rank2"] <- sapply(tax_table(phyl_L6props)[,"Rank2"],
                                          function(s){if (grepl('__', s)) 
                                            unlist(strsplit(s, '__'))[2] else s})
  tax_table(phyl_L6props)[,'Rank2'] <- apply(tax_table(phyl_L6props), 1, applyNewGenusLabel)
  # 'Rank2' column now holds genus-label value instead of phylum
  psL6summary <- tax_glom(phyl_L6props, 'Rank2')
  summaryprop <- otu_table(psL6summary)
  colnames(summaryprop) <- as.data.frame(tax_table(psL6summary))[, "Rank2"]
  summarypropdf <- data.frame(summaryprop[, intersect(genusOrder, colnames(summaryprop))], 
                              # intersect drops items not in list, keeps in list order, and is OK if any missing
                              # if(taxa_are_rows(summaryprop)) this won't work, need to transpose
                              personID=sample_data(psL6summary)$personID,
                              DayXMethod=paste0(sample_data(psL6summary)$Stool,
                                                sample_data(psL6summary)$Method)  )
  summarypropdf$Otherbacteria <- apply(summaryprop[, setdiff(colnames(summaryprop), genusOrder)], 
                                       1-taxa_are_rows(summaryprop), sum)
  longProp <- melt(summarypropdf, value.name='propAbund',
                   id.vars=c("DayXMethod", "personID"),
                   variable.name='genuslabel')
  bars <- ggplot(longProp, aes(x=DayXMethod, y=propAbund, fill=genuslabel))
  bars + geom_bar(stat="identity") + 
    facet_wrap(~personID) + 
    labs(x="Day and method",
         y="Proportional abundance", 
         title="") +
    scale_fill_manual(values = palGenus, 
                      guide=guide_legend(nrow=4, title = "", keyheight = 0.5)) + 
    theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5)
          , strip.background=element_rect(fill=NA) 
          , legend.key = element_rect(size=NA)
          , legend.position = "bottom"
    )

}


orderBoxplot <- function(phyl_L6counts,
                         phylOrd = phylOrder , 
                         palPhyl = palette_phyl) {
  psL4propNo0 <- transform_sample_counts(tax_glom(phyl_L6counts, taxrank="Rank4") , 
                                         function(x) ( (x+1)/sum(x) ))
  tax_table(psL4propNo0)[,"Rank2"] <- sapply(tax_table(psL4propNo0)[,"Rank2"],
                                             function(s){unlist(strsplit(s, '__'))[2]})
  tax_table(psL4propNo0)[,"Rank4"] <- sapply(tax_table(psL4propNo0)[,"Rank4"],
                                             function(s){unlist(strsplit(s, '__'))[2]})
  ## This uses the phyloseq 'plot_bar' function to set the ggplot data and aesthetics. 
  ## Two boxplots are added, 1st with border and point colours from Phylum, including 
  ## the default solid circle outliers. The 2nd has black borders, open circle outliers 
  ## and coloured fill. This has the effect of colouring outliers to match bars.
  ## The original barplot layer is then removed
  pbox <- plot_bar(psL4propNo0, x="Rank4") + scale_y_log10() + 
    geom_boxplot(aes(colour=Rank2)) + 
    geom_boxplot(aes(fill=Rank2), outlier.shape=21) + 
    labs(title="Proportion of bacterial Orders in each sample (log scale)", y="Proportion",
         x=element_blank()) +
    scale_fill_manual(values=palPhyl, name="Phylum") + 
    scale_colour_manual(values=palPhyl, name="Phylum") + 
    theme(axis.text.x = element_text(vjust=0.5))
  pbox$layers <- pbox$layers[-1]
  ## X-axis goes in order of decreasing total sum of proportions
  barOrder <- tax_table(psL4propNo0)[names(sort(taxa_sums(psL4propNo0), TRUE)),"Rank4"]
  pbox$data$Rank4 <- factor(pbox$data$Rank4, levels=barOrder) 
  pbox$data$Rank2 <- factor(pbox$data$Rank2, levels=phylOrd) 
  pbox  
}

avg_log_transform_fn <- function(phylobject, rank = "Rank2", topNum=4, rankOrder=phylOrder){
  ## Agglomerate counts to Rank 2, Phylum, discard least abundant 
  ## and transform to natural log + 1
  rank_glom <- tax_glom(phylobject, taxrank = rank)
  ## Abbreviate phylum names by removing leading "p__"
  tax_table(rank_glom) <- apply(tax_table(rank_glom), 1:2, 
                                        function(s){if (grepl('__', s))
                                          unlist(strsplit(s, '__'))[2] else
                                            s})
  ## Put sample data in a dataframe
  samplesPs <- data.frame(sample_data(rank_glom),
                            stringsAsFactors = FALSE)
  ## Put counts in a dataframe
  if (taxa_are_rows(rank_glom)) {
    phylacounts <- data.frame(t(otu_table(rank_glom)))
  } else {
    phylacounts <- data.frame(otu_table(rank_glom))
  }
  colnames(phylacounts) <- tax_table(rank_glom)[,rank]
  ## Combine counts and sample_data, keep only topNum
  top_rank_df <- data.frame(phylacounts[,rankOrder[1:topNum]],
                            samplesPs
  )
  ## Convert phyla as column names to column of phyla names in long format
  top_rank_long <- melt(top_rank_df, measure.vars =c(1:topNum), variable.name = rank)
  top_rank_long$log_yij <- log1p(top_rank_long$value)  # using ln(value + 1) 
  require(plyr)
  top_rank_summary <- ddply(top_rank_long,
                            .(Rank2, personID, Method), summarise,
                            meanij=mean(value), sdij=sd(value),
                            meanlnij=mean(log_yij), varlnij=var(log_yij),
                            num_reps = length(log_yij))
}

log_sd_xPerson_fPhyla_plot <- function(stats_per_cell_df){
  ggplot(data=stats_per_cell_df,
         aes(x=personID, y=meanlnij, colour=Method, shape=Method)) +
    geom_point(position=position_dodge(width=0.5)
               , size=1
    ) +
    geom_errorbar(aes(ymin=meanlnij-sqrt(varlnij), ymax=meanlnij+sqrt(varlnij))
                  , position=position_dodge(width=0.5)
                  , width=0.5
                  , alpha=0.6
                  , show.legend = FALSE) +
    facet_wrap(~Rank2, scales='free_y') +
    labs(y='log(standardised count)' ) +
    scale_color_brewer(palette = 'Set1') +
    #scale_colour_discrete(name='Method', palette='Set1') +
    theme(legend.key.size = unit(1, "lines")
          , legend.key = element_rect(size=NA)
          , legend.position = "bottom"
          , legend.margin = unit(0, "cm")
          , strip.background=element_rect(fill=NA)
    )
}


############################################################################

### Local-local New trimmed data 29 Aug ###
baseDir <- file.path("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/metagenomics",
                     "ENDIA/ENDIA_QC")
imageDir <- file.path(baseDir, "outputs/Aug16_trimmed_seqs/taxabund_plots") # WEHI trimmed
localRdata_fp <- file.path(baseDir, "analysis_tools/endia_biomR",
                           "local_data/New_WEHI_aligned_trimmed")

## Load the proportional counts for 72 samples: 
## (bar plots sum values by default - using 215 samples gives sum of 2 or 3 values per bar)
load(file.path(localRdata_fp, "Stdcounts_and_props.rdata"), verbose=TRUE) 

## Genus barplot (L6), list above of genus, or higher label where genus undefined,
## present in at least 10% in any sample. Fig 2A
genuslabelbars <- genusBarplot(psL6prop) # other parameters as default
## Save to pdf
plot_name <- file.path(imageDir, "barPropGenus_phylaColours84w.pdf") 
genuslabelbars
ggsave( filename= plot_name,  
        width=84, height=110, units="mm")

### Order (L4) box-and-whiskers, with log-scale ####
## Fig S3A
boxes <- orderBoxplot(phyl_L6counts = l6Std) 
print(boxes)
ggsave( filename= file.path(imageDir, "boxPropOrderColourPhylumLvl2.pdf"), 
        width=297, height=210, units="mm")
##-----------

### Load the proportional counts for 215 samples. ###
## Box-and whiskers above can be re-done with more points if desired. 
## Faceted plots should be done as geom_points if using 72 points 
## as only 3 stools-samples per Method.personID entry
load(file.path(localRdata_fp, "Stdcounts_and_propsLvl2.rdata"), verbose=TRUE) 

### Mean and sd at the Lvl2 (aliquot) level, faceted by Phylum ###
## Top 4 phyla, log-transformed standardised counts for 215 samples
## Fig 3A 
phyla_stats_summary <- avg_log_transform_fn(psOTUstd)
log_sd_plot <- log_sd_xPerson_fPhyla_plot(phyla_stats_summary)
log_sd_plot
ggsave( filename= file.path(imageDir, "ln_std_devColMeth84.85.pdf"),
        width=84, height=85, units="mm")

############################################################################

### Baylor-local ###
baseDir <- file.path("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/metagenomics",
                     "ENDIA/ENDIA_QC/baylordata_and_analysis/")
localRdata_fp <- "baylor_data"
imageDir <- file.path(baseDir, "analysis/local/taxabund_plots/")
## Load the proportional counts for 72 samples
load(file.path(localRdata_fp, "Stdcounts_and_props.rdata"), verbose=TRUE) 

## Genus barplot (L6), list above of genus, or higher label where genus undefined,
## present in at least 10% in any sample. 
## Fig 4B
genuslabelbars <- genusBarplot(psL6prop) # other parameters as default
## Save to pdf
plot_name <- file.path(imageDir, "barPropGenus_phylaColours84w_halfheightkey.pdf") 
genuslabelbars
ggsave( filename= plot_name,  
        width=84, height=110, units="mm")

##-----------
load(file.path(localRdata_fp, "Stdcounts_and_propsLvl2.rdata"), verbose=TRUE) 
### Order (L4) box-and-whiskers, with log-scale. Supplementary Fig. S3B
boxes <- orderBoxplot(phyl_L6counts = l6Std) 
print(boxes)
ggsave( filename= paste0(imageDir, "boxPropOrderColourPhylumLvl2.pdf"), 
        width=297, height=210, units="mm")

### Mean and sd at the Lvl2 (aliquot) level, faceted by Phylum  ###
## Top 4 phyla, log-transformed standardised counts for 215 samples. Suppl Fig. 6C
phyla_stats_summary <- avg_log_transform_fn(l6Std)
log_sd_plot <- log_sd_xPerson_fPhyla_plot(phyla_stats_summary)
log_sd_plot
ggsave( filename= file.path(imageDir, "ln_std_devColMeth84.85.pdf"), 
        width=84, height=85, units="mm")
