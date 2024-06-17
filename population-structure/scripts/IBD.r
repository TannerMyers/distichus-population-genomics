setwd("/Users/tannermyers/Dropbox/Distichus_Project/ddRADseq_Phylogeography/AnoDist_aln/population-structure/Fst")

library(tidyverse)
library(cowplot)
library(vcfR)
library(ade4)
library(adegenet)
library(geosphere)
library(LEA)
library(MASS)
library(eks)
library(colorspace)
library(multcomp)
library(MuMIn)

## Define plotting function for mantel.randtest output
ggplot_mantel <- function(mant, fill = "gray50") {
        df <- data.frame(x = mant$plot$hist$mids,
                        y = mant$plot$hist$counts)
        ggplot(df, aes(x, y)) + 
        geom_col(orientation = "x", 
                    width = diff(mant$plot$hist$breaks)[1],
                    fill = fill, color = "gray30") +
        labs(x = 'Permutation', y = "Frequency", title = "Mantel Test Significance") +
        scale_x_continuous(limits = mant$plot$xlim) +
        geom_segment(aes(x = mant$obs, xend = mant$obs, y = 0,
                            yend = 0.75 * max(y), colour = "red")) +
        geom_point(aes(x = mant$obs, y = 0.75 * max(y)), size = 5,
                    shape = 18, colour = "red") +
        theme(axis.text.x = element_text(colour = "black")) +
        theme(legend.title = element_blank()) + theme(axis.text.y = element_text(colour = "black", size = 12)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(colour = "black", size = 12)) +
        theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
    }

vcf <- read.vcfR("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/AnoDist_aln/distichus_geno90_filtered_LD_50kb_100_0.1_biallelic_pruned.recode.vcf")  

## Load specimen information into R
columns <- c("specimen", "specimen_pop", "filepath", "subspecies", "population", "complex", "island", "locality", "longitude", "latitude", "note")
RAD_data <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE) 
RAD_data_north <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/north-island-dataset-master.tsv", col_names = FALSE)
RAD_data_south <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/south-island-dataset-master.tsv", col_names = FALSE)
    colnames(RAD_data) <- columns
    colnames(RAD_data_north) <- columns
    colnames(RAD_data_south) <- columns

coords <- distinct(RAD_data[, c("locality", "longitude", "latitude")])

colnames(vcf@gt)[2:192] <- RAD_data$specimen_pop
data <- vcfR2genind(vcf)
pop(data) <- RAD_data$locality

pop_data <- genind2genpop(data)

gen_dist <- dist.genpop(pop_data, method = 1) # 1 for Nei's D, 2 for Edwards' distance, 

geo_dist <- distm(coords[, 2:3], fun = distGeo)
geo_dist <- as.dist(geo_dist)

ibd <- mantel.randtest(m1 = gen_dist, m2 = geo_dist, nrepet = 100000)
p_mantel <- ggplot_mantel(ibd)
col2_1 <- plot_grid(NULL, p_mantel, labels = "", ncol = 1)


# taken directly from adegenet tutorial to create 2-D kernel density estimate
pdf("IBD-distichus-subsp.pdf")
    dens <- MASS::kde2d(geo_dist, gen_dist, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(geo_dist, gen_dist, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(as.numeric(gen_dist) ~ as.numeric(geo_dist)))
    #title("Isolation by distance plot")
dev.off()

dists <- data.frame(cbind(geo_dist, gen_dist))
dens2 <- tidy_kde(dists)
p_dens <- ggplot(dens2, aes(x = geo_dist, y = gen_dist)) +
    geom_raster(aes(fill = estimate), interpolate = TRUE) +
    labs(x = "Geographic Distance", y = "Genetic Distance (Nei's D)") +
    scale_fill_continuous_sequential(palette = "YlGnBu") + geom_point_ks(alpha = 0.1) +
    theme(axis.text.x = element_text(colour = "black")) +
    theme(axis.text.y = element_text(colour = "black", size = 12)) +
    theme(legend.title = element_blank(), legend.position = "none") +
    theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
p_dens + geom_smooth(data = dens2[, 1:2], method = "lm", se = FALSE)

row1 <- plot_grid(p_dens, col2_1, labels = "", rel_widths = c(3,1), rel_heights = c(3,1))


# Now do for north and south paleo island populations
n_vcf <- read.vcfR("../north_geno90_50kb_100_0.1_biallelic.recode.vcf")
    colnames(n_vcf@gt)[2:144] <- RAD_data_north$specimen_pop
s_vcf <- read.vcfR("../south_geno90_50kb_100_0.1_biallelic.recode.vcf")
    colnames(s_vcf@gt)[2:49] <- RAD_data_south$specimen_pop

n_coords <- distinct(RAD_data_north[, c("locality", "longitude", "latitude")])
s_coords <- distinct(RAD_data_south[, c("locality", "longitude", "latitude")])

n_data <- vcfR2genind(n_vcf)
s_data <- vcfR2genind(s_vcf)

pop(n_data) <- RAD_data_north$locality
pop(s_data) <- RAD_data_south$locality

n_pop <- genind2genpop(n_data)
s_pop <- genind2genpop(s_data)

n_gen_dist <- dist.genpop(n_pop, method = 1)
s_gen_dist <- dist.genpop(s_pop, method = 1) # Nei's D

n_geo_dist <- distm(n_coords[, 2:3], fun = distGeo)
n_geo_dist <- as.dist(n_geo_dist)

s_geo_dist <- distm(s_coords[, 2:3], fun = distGeo)
s_geo_dist <- as.dist(s_geo_dist)

# north
n_ibd <- mantel.randtest(m1 = n_gen_dist, m2 = n_geo_dist, nrepet = 100000)
col2_2 <- plot_grid(NULL, ggplot_mantel(n_ibd), labels = "", ncol = 1)

# south
s_ibd <- mantel.randtest(m1 = s_gen_dist, m2 = s_geo_dist, nrepet = 100000)
col2_3 <- plot_grid(NULL, ggplot_mantel(s_ibd), labels = "", ncol = 1)

# north
pdf("IBD-north-distichus-subsp.pdf")
    dens <- kde2d(n_geo_dist, n_gen_dist, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(n_geo_dist, n_gen_dist, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(as.numeric(n_gen_dist) ~ as.numeric(n_geo_dist)))
    #title("Isolation by distance plot")
dev.off()

# Or, in ggplot, do
dists <- data.frame(cbind(n_geo_dist, n_gen_dist))
dens2 <- tidy_kde(dists)
n_dens <- ggplot(dens2, aes(x = n_geo_dist, y = n_gen_dist)) +
    geom_raster(aes(fill = estimate), interpolate = TRUE) +
    labs(x = "Geographic Distance", y = "Genetic Distance (Nei's D)") +
    scale_fill_continuous_sequential(palette = "YlGnBu") + geom_point_ks(alpha = 0.1) +
    theme(axis.text.x = element_text(colour = "black")) +
    theme(axis.text.y = element_text(colour = "black", size = 12)) +
    theme(legend.title = element_blank(), legend.position = "none") +
    theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
#n_dens + geom_smooth(data = dens2[,1:2], method = "lm", se = FALSE)

row2 <- plot_grid(n_dens, col2_2, labels = "", rel_widths = c(3,1), rel_heights = c(3,1))

# south
pdf("IBD-south-distichus-subsp.pdf")
    dens <- kde2d(s_geo_dist, s_gen_dist, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(s_geo_dist, s_gen_dist, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(as.numeric(s_gen_dist) ~ as.numeric(s_geo_dist)))
    #title("Isolation by distance plot")
dev.off()

dists <- data.frame(cbind(s_geo_dist, s_gen_dist))
dens2 <- tidy_kde(dists)
s_dens <- ggplot(dens2, aes(x = s_geo_dist, y = s_gen_dist)) +
    geom_raster(aes(fill = estimate), interpolate = TRUE) +
    labs(x = "Geographic Distance", y = "Genetic Distance (Nei's D)") +
    scale_fill_continuous_sequential(palette = "YlGnBu") + geom_point_ks(alpha = 0.1) +
    theme(axis.text.x = element_text(colour = "black")) +
    theme(axis.text.y = element_text(colour = "black", size = 12)) +
    theme(legend.title = element_blank(), legend.position = "none") +
    theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
#s_dens + geom_smooth(data = dens2[,1:2], method = "lm", se = FALSE)

row3 <- plot_grid(s_dens, col2_3, labels = "", rel_widths = c(3,1), rel_heights = c(3,1))

# Final Plot
plot_grid(row1, NULL, row2, NULL, row3, labels = c("A", "", "B", "", "C"), label_size = 24, ncol = 1, rel_heights = c(1, 0, 0.8, 0, 0.8))
ggsave("~/Dropbox/Apps/Overleaf/Anolis distichus Population Genomics ms/supp-mat/mantel-test-hierarchical.pdf", height = 24, width = 16)