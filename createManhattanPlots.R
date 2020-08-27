# This program exists to build a Manhattan plot provided a summary statistics
# input file. Requires that we add positions to the summary stats prior to
# running

library(ggplot2)
library(ggrepel)
library(magrittr)
library(dplyr)

gwasResults <- read.table("SourceData/BMI_2015_summaryStatPositions_noXY.tsv", header = TRUE)

# Remove high p-value results to reduce the time to create these
# plots
gwasResults <- gwasResults %>% 
  filter(-log10(pvalue)>1)

outputFilename <- "2015_BMI_GWAS-Skeleton.png"

# Function to get the cumulative position of the SNPs
don <- gwasResults %>%

  # Compute chromosome size
  group_by(chr) %>% summarise(chr_len=max(pos)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("chr"="chr")) %>%

  # Add a cumulative position of each SNP
  arrange(chr, pos) %>% mutate( poscum=pos+tot)

axisdf <- don %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

png(outputFilename,width=1000,height=500)

ggplot(don, aes(x=poscum, y=-log10(pvalue))) + 
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(c("grey"), 23 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$'chr', breaks= axisdf$center ) +

    ylim(0, 20) + 

    # Custom the theme:
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", color = "black", size = 16, hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.y = element_text(color = "black", size = 14, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) 

dev.off()
