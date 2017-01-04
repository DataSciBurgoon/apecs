#apecs_in_vivo.R
#
# Autoencoder Predicting Estrogenic Chemical Substances (APECS) -- in vivo 
#
# Author: Lyle D. Burgoon, Ph.D.
# US Army Research and Development Center
# License: CC0/Public Domain
#
# Description:
# This script predicts if a chemical is estrogenic or not based on the in vivo animal
# estrogenicity database provided in Browne et al 
# (2015, Environmental Science and Technology, 49(14): 8804-8814) using the 
# ToxCast assay data. Data can be used from the same in vitro assays as ToxCast 
# so long as the input data meet the same structure.
#
# Version: 1.0.0  4-JAN-2017

#arg_params:
#1 input_file with toxcast assay data
#2 [optional] memory size (default is 2GB)
#3 [optional] number of threads (default is 2)


#calculates a nonlinear regression on the data
loess_calc <- function(chem_data){
  chem_data.lo <- loess(mean_resp ~ logc, chem_data)
  new_logc <- approx(chem_data$logc, n=50)$y
  new_resp <- predict(chem_data.lo, data.frame(logc=new_logc))
  return(new_resp)
}

#calculates the data for 17-beta estradiol
comparator_estimates <- function(assay, x, comparator = "17beta-Estradiol"){
  chem_data <- filter(x, chnm == as.character(comparator) & assay_name == as.character(assay))
  #print(filter(x, chnm == as.character(comparator), assay_name == as.character(assay)))
  chem_estimate <- loess_calc(chem_data)
  return(chem_estimate)
}

#x: vector containing chemical name and assay name
#comparator: this is the chemical that we're comparing everything to
#data: the toxcast data, per_chem_estrogen_assay_logc
similarity_analysis <- function(x, comparator_estimates, data){
  chem_name <- as.character(x[1])
  assay <- as.character(x[2])
  chem_data <- filter(data, chnm == chem_name & assay_name == assay)
  if(nrow(chem_data) > 3){
    chem_estimate <- loess_calc(chem_data)
    chem_correlation <- cor(chem_estimate, comparator_estimates[, assay])
    return(c(chem_name, assay, chem_correlation))
  }
}

list_of_packages <- c("ggplot2", "h2o")

#Test to see if ggplot2 and h2o are installed.
#If not, then install them.

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='http://archive.linux.duke.edu/cran/')
library(ggplot2)
library(h2o)
#library(gRain)
#library(ROCR)
#library(pracma)
#library(dplyr)
#library(plyr)
#library(tidyr)

arg_params = commandarg_params(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(arg_params)==0) {
  arg_params[1] <- "data.txt"
  arg_params[2] <- "2g"
  arg_params[3] <- 2
} else if(length(arg_params)==1 && class(arg_params[1]) == "character"){
  arg_params[2] <- "2g"
  arg_params[3] <- 2
}

toxcast_full <- read.csv(file=arg_params[1], header=TRUE, sep="\t")

#chem_cas_count <- ddply(toxcast_full, .(chnm, casn), summarise, count=length(casn))

#save(per_chem_assay_logc, file="per_chem_assay_logc.RData")

in_vitro_per_chem_assay_logc <- ddply(toxcast_full, .(chnm, assay_name, logc, gene_name), 
                                      summarise, mean_resp=mean(resp))

#Toxcast assays for the ER positive ground truth
aenm_list <- c("NVS_NR_hER", "OT_ER_ERaERa_0480", "OT_ER_ERaERa_1440", 
               "OT_ER_ERaERb_0480", "OT_ER_ERaERb_1440", "OT_ER_ERbERb_0480", "OT_ER_ERbERb_1440",
               "TOX21_ERa_BLA_Agonist", "ATG_TRANS", "ATG_CIS")

in_vitro_per_chem_estrogen_assay_logc <- filter(in_vitro_per_chem_assay_logc, assay_name %in% aenm_list, 
                                                gene_name == "estrogen receptor 1" | gene_name == "estrogen receptor 2 (ER beta)")

in_vitro_distinct_chem_assay <- distinct(select(in_vitro_per_chem_estrogen_assay_logc, chnm, assay_name))
in_vitro_distinct_assay <- distinct(select(in_vitro_per_chem_estrogen_assay_logc, assay_name))

load("in_vitro_estradiol_curves.RData")

#Fun Fact:
#The first time I ran this apply I got a warning that said "Chernobyl!"
#I looked into why, and finally found the answer:
#..."eigenvalue meltdown occurs and the coded message  Chernobyl  is raised 

#In VITRO
chem_assay_corr_estrogen <- apply(in_vitro_distinct_chem_assay, 1, similarity_analysis, 
                                  comparator_estimates=in_vitro_estradiol_curves, 
                                  data=in_vitro_per_chem_estrogen_assay_logc)

chem_assay_corr_estrogen_df <- do.call(rbind.data.frame, chem_assay_corr_estrogen)
colnames(chem_assay_corr_estrogen_df) <- c("Chemical", "Assay", "Correlation")

chem_assay_corr_estrogen_wide <- spread(chem_assay_corr_estrogen_df, 
                                        Assay, Correlation)

chem_assay_corr_estrogen_wide[is.na(chem_assay_corr_estrogen_wide)] <- 0        #Convert NAs to 0
row_names <- as.character(chem_assay_corr_estrogen_wide[,1])
chem_assay_corr_estrogen_wide <- chem_assay_corr_estrogen_wide[,-1]
in_vitro_chem_assay_corr_estrogen_wide <- as.data.frame(sapply(chem_assay_corr_estrogen_wide, as.numeric))
rownames(in_vitro_chem_assay_corr_estrogen_wide) <- row_names
in_vitro_chem_assay_corr_estrogen_wide[nrow(in_vitro_chem_assay_corr_estrogen_wide) + 1, ] <- rep(1, ncol(in_vitro_chem_assay_corr_estrogen_wide))
rownames(in_vitro_chem_assay_corr_estrogen_wide)[nrow(in_vitro_chem_assay_corr_estrogen_wide)] <- "17-beta estradiol"

#localH2O = h2o.init(max_mem_size = "13g", nthreads=6)
localH2O = h2o.init(max_mem_size = arg_params[2], nthreads=as.numeric(arg_params[3]))
estrogenic.hex<-as.h2o(in_vitro_chem_assay_corr_estrogen_wide, destination_frame="train.hex")
estrogenic_autoencoder <- h2o.loadModel("DeepLearning_model_R_1481832873762_1608")
estrogenic_supervised_features3 <- h2o.deepfeatures(estrogenic_autoencoder, estrogenic.hex, layer=2)


# estrogenic_autoencoder = h2o.deeplearning(
#   x = 1:ncol(in_vitro_chem_assay_corr_estrogen_wide),
#   training_frame = estrogenic.hex,
#   hidden = c(10, 2, 10),
#   activation="Tanh",
#   epochs = 20000,
#   autoencoder = TRUE
# )


#This is the test error

plotdata2 <- as.data.frame(estrogenic_supervised_features3)
plotdata2$chem <- rownames(in_vitro_chem_assay_corr_estrogen_wide[50:65, ])

estradiol_coords <- plotdata2[which(plotdata2$chem=="17-beta estradiol"), c(1:2)]
angle <- seq(-pi, pi, length = 50)
circle_boundary <- data.frame(x = (sin(angle)*1.35)+estradiol_coords$DF.L2.C1, 
                              y = (cos(angle)*1.35)+estradiol_coords$DF.L2.C2)

plotdata2$id <- c(1:nrow(plotdata2))

png("all_assays_autoencoder_plot_colors_in_vitro_with_boundary.png", width=9500, height=9500, res=1200)

ggplot(plotdata2) +
  xlim(c(-2, 2)) +
  geom_point(aes(x=DF.L2.C1, y=DF.L2.C2, label=chem, size=1)) +
  geom_text(aes(x=DF.L2.C1+0.03, y=DF.L2.C2+0.03, label=chem),hjust=0, vjust=0) +
  #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
  #             data = ae_angle_drivers,
  #             arrow=arrow()) +
  geom_path(aes(x, y), data = circle_boundary, inherit.aes = F) +
  guides(size=FALSE)


dev.off()


#Distance in Autoencoder 2-D space
estradiol_row <- which(plotdata2$chem == "17-beta estradiol")
in_vitro_dist <- sqrt((plotdata2[,1] - plotdata2[estradiol_row,1])^2 + (plotdata2[, 2] - plotdata2[estradiol_row,2])^2)
in_vitro_dist_df <- data.frame(plotdata2, dist=in_vitro_dist)

in_vitro_active_dist <- subset(in_vitro_dist_df, dist <= 1.35)
in_vitro_inactive_dist <- subset(in_vitro_dist_df, dist > 1.35)

in_vitro_active_dist$estrogenic <- "estrogenic"
in_vitro_inactive_dist$estrogenic <- "not estrogenic"

in_vitro_activity <- rbind(in_vitro_active_dist, in_vitro_inactive_dist)
write.table(in_vitro_activity[, c(3,5,6)], file="in_vitro_activity.txt", sep="\t", row.names=FALSE)
h2o.shutdown()
