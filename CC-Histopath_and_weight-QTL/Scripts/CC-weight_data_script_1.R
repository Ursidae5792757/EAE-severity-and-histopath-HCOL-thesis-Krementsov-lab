### CC-weight_data_script_1
# This script is designed to take weight data from the CC EAE project and perform QTL mapping on it. This script is heavily based on https://smcclatchy.github.io/mapping/ and Dr. Dimitry Krementsov's QTL mapping pipeline
# Author: Lauren Downs
# Started: 8/31/2025

## Importing the library
library(qtl2)

## Importing CC data
setwd("/Users/Lauren Downs/Documents/UVM Krementsov Lab/QTL Mapping/CC-Histopath_and_weight-QTL") #used to get the following line to act right
cc1 <-read_cross2("Input/cc/cc.json")
names(cc1)
summary(cc1)


## Exploring CC data
#cc1$covar
#summary(cc1$covar)
#cc1$is_female

## Histogram - checking for normal
eae1_pheno <- read.csv("Input/cc/eae1_pheno.csv")
hist(eae1_pheno$Global_CDS, main = "CDS")
hist(eae1_pheno$net_area_under_curve_weight_data, main = "Weigth")

## Making a map with pseudomarkers (as per the tutorial linked above)
cc_map <- insert_pseudomarkers(map=cc1$gmap, step=1)


## Calculating Gene Probabilities
gene_pr <- calc_genoprob(cc1,error_prob=0.002) # the DK way
gene_pr_psdo <- calc_genoprob(cross=cc1, map = cc_map,error_prob=0.002) # the tutorial way

#looking at the probs
names(gene_pr)
names(gene_pr_psdo)
gene_pr$'1'
gene_pr_psdo$'1'


## Calculating Allele Probabilities
all_pr <- genoprob_to_alleleprob(gene_pr)

#looking at the probs
names(all_pr)
all_pr$'1'


## Calculating Kinship (LOCO version only)
kinship_all_pr <- calc_kinship(all_pr, 'loco')
kinship_gene_pr <- calc_kinship(gene_pr,'loco')
kinship_gene_pr_psdo <- calc_kinship(gene_pr_psdo, 'loco')

#looking at kinship
chromosome <- 16
heatmap(kinship_all_pr[[chromosome]]) 
heatmap(kinship_gene_pr[[chromosome]]) 
heatmap(kinship_gene_pr_psdo[[chromosome]]) 


## Sex difference covariance 
x_covar <- get_x_covar(cc1)
head(x_covar) # apparently this doesn't work properly because it returns NULL


## QTL Scanning
all_pr_scan <- scan1(all_pr, pheno = cc1$pheno)
gene_pr_scan <- scan1(gene_pr, pheno = cc1$pheno)
gene_pr_psdo_scan <- scan1(gene_pr_psdo, pheno = cc1$pheno)

#plotting
plot_scan1(all_pr_scan, map= cc_map,lodcolumn = "net_area_under_curve_weight_data",col = 'blue')
plot_scan1(all_pr_scan, map= cc_map,lodcolumn = "Global_CDS",col = 'green', add= TRUE) # using two columns (blue for weight data and green for global CDS)

plot_scan1(gene_pr_scan, map= cc_map)
plot_scan1(gene_pr_psdo_scan, map= cc_map) #these all seem to look the same


## Finding Peaks
scan_output <- all_pr_scan # here you can choose from above methods for getting a scan
pr_of_interest <-all_pr # likewise for choosing the complimentary probability

#doing permutations for comparison against random chance
perm_output <- scan1perm(genoprobs = pr_of_interest, pheno = cc1$pheno, Xcovar = x_covar, n_perm = 1000) # replace 1000 with 10 for expediency

#viewing permutations
column <- 'Global_CDS'
hist(perm_output[,column], breaks = 50, xlab = "LOD", main = sprintf("LOD scores for %s scan with threshold in red", column)) #THIS changes based on what method used
abline(v = summary(perm_output)[,column], col = 'red', lwd = 2)

#actually finding peaks
thr <- summary(perm_output,alpha=c(0.2))
find_peaks(scan1_output = scan_output, map = cc_map, threshold = thr)

#plotting again
ylim <-c(0, maxlod(scan_output)*1.4 ) # borrowed from the internet https://rdrr.io/github/rqtl/qtl2plot/man/plot_scan1.html
plot_scan1(scan_output, map= cc_map,lodcolumn = "Global_CDS",col = 'green', ylim =ylim,main = "LOD scores for scan with threshold in matching colors") # using two columns (blue for weight data and green for global CDS)
plot_scan1(scan_output, map= cc_map,lodcolumn = "net_area_under_curve_weight_data",col = 'blue',add= TRUE)
abline(h = summary(perm_output)[,"Global_CDS"], col = 'darkgreen', lwd = 2)
abline(h = summary(perm_output)[,"net_area_under_curve_weight_data"], col = 'darkblue', lwd = 2)

