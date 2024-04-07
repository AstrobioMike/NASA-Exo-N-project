library(tidyverse)
library(KEGGREST)
library(stringi)

## helper functions at bottom ##

setwd("~/Documents/NASA/Nitrogen-project/large-exp/metagenomics/combined-outputs")

all_cov <- read.table("All-combined-KO-function-coverages.tsv", sep = "\t", header = TRUE, quote = "")
head(all_cov)
dim(all_cov)

# reading in sample info table
sample_info_tab <- read.table("../sample-and-group-info.tsv", sep = "\t", header = TRUE, comment.char = "")


# making KEGG summary table with all info for the KOs we have
all_our_KOs <- all_cov %>% pull(KO_ID) %>% as.vector()

Master_KEGG_tab <- make_KO_summary_tab(all_our_KOs)
write.table(Master_KEGG_tab, "Master-kegg-annotations.tsv", sep="\t", quote=F, row.names=F)



# choosing based on archaea and bacteria from KEGG's page here: https://www.genome.jp/kegg/annotation/br01610.html
target_SCGs <- c("L3" = "K02906",
                 "L23" = "K02892",
                 "L2" = "K02886",
                 "L22" = "K02890",
                 "L29" = "K02904",
                 "L14" = "K02874",
                 "L24" = "K02895",
                 "L5" = "K02931",
                 "L6" = "K02933",
                 "L18" = "K02881",
                 "L30" = "K02907",
                 "L15" = "K02876",
                 "L13" = "K02871",
                 "L10" = "K02864",
                 "L1" = "K02863",
                 "L11" = "K02867",
                 "S10" = "K02946",
                 "S19" = "K02965",
                 "S3" = "K02982",
                 "S17" = "K02961",
                 "S14" = "K02954",
                 "S8" = "K02994",
                 "S5" = "K02988",
                 "S13" = "K02952",
                 "S11" = "K02948",
                 "S4" = "K02986",
                 "S9" = "K02996",
                 "S7" = "K02992",
                 "S12" = "K02950",
                 "S2" = "K02967",
                 "S15" = "K02956")

length(target_SCGs)

# median of all
apply(all_cov[,c(3:dim(all_cov)[2])], 2, median)

# median of target SCGs
SCGs_cov <- all_cov %>% filter(KO_ID %in% target_SCGs)
apply(SCGs_cov[,-1], 2, median)

apply(SCGs_cov[,-1], 2, summary)


# median of individual target SCGs across samples (though this would be directly affected by library size, maybe the distributions would look similar for good candidate final SCGs to use?)

plot_hist_of_SCG_covs <- function(target_SCG, cov_tab = all_cov) {

    curr_covs <- cov_tab %>% filter(KO_ID == target_SCG) %>% select(-c(KO_ID, KO_function)) %>% as.numeric
    curr_summary <- summary(curr_covs)
    print(curr_summary)

    return(hist(curr_covs, main = target_SCG, xlab = "Coverages"))
}

   # those with the same number of "#"s have similar distributions
plot_hist_of_SCG_covs("K02906") #
plot_hist_of_SCG_covs("K02892") #
plot_hist_of_SCG_covs("K02886") #
plot_hist_of_SCG_covs("K02890") #
plot_hist_of_SCG_covs("K02904") #
plot_hist_of_SCG_covs("K02874") ##
plot_hist_of_SCG_covs("K02895") ##
plot_hist_of_SCG_covs("K02931") ##
plot_hist_of_SCG_covs("K02881") #
plot_hist_of_SCG_covs("K02907") #
plot_hist_of_SCG_covs("K02876") ##
plot_hist_of_SCG_covs("K02871") #
plot_hist_of_SCG_covs("K02864") #
plot_hist_of_SCG_covs("K02863") ##
plot_hist_of_SCG_covs("K02867") #

plot_hist_of_SCG_covs("K02946") #
plot_hist_of_SCG_covs("K02965") #
plot_hist_of_SCG_covs("K02982") ##
plot_hist_of_SCG_covs("K02961") ##
plot_hist_of_SCG_covs("K02954") ###
plot_hist_of_SCG_covs("K02994") #
plot_hist_of_SCG_covs("K02988") #
plot_hist_of_SCG_covs("K02952") #
plot_hist_of_SCG_covs("K02948") ##
plot_hist_of_SCG_covs("K02986") #
plot_hist_of_SCG_covs("K02996") #
plot_hist_of_SCG_covs("K02992") #
plot_hist_of_SCG_covs("K02950") #
plot_hist_of_SCG_covs("K02967") #
plot_hist_of_SCG_covs("K02956") #


SCGs_cov <- all_cov %>% filter(KO_ID %in% target_SCGs) %>% select(-2)

# want to see histograms of all SCGs in a single sample?
hist(SCGs_cov %>% pull(F4N1))
hist(SCGs_cov %>% pull(F4N5))

summary(SCGs_cov %>% pull(F4N1))
summary(SCGs_cov %>% pull(F4N5))



# summaries of all
apply(SCGs_cov[, -1], 2, summary)

SCGs_cov_long <- SCGs_cov %>% pivot_longer(!KO_ID, names_to = "Sample", values_to = "coverage")

ggplot(SCGs_cov_long, aes(coverage)) + geom_histogram() + facet_wrap(~Sample, scales = "free")


## want some summary stats including coefficient of variation of all SCGs in each individual sample
## thinking something like
# stat  F4N1    F4N5
# Min
# 1st Q
# Med
# Mean
# 3rd Q
# Max
# SD
# CoV

## Or maybe transposed
# Sample    Min     1stQ    Med     Mean    3rdQ    Max     SD      CoV
# F4N1
# F4N5

##########
all_samples <- names(SCGs_cov)[-1]

SCG_cov_summary_tab <- data.frame(row.names = c("Min", "1stQ", "Med", "Mean", "3rdQ", "Max", "SD", "CoV"))

for ( sample in all_samples ) {

    # getting current sample SCG coverages
    curr_vec <- SCGs_cov %>% pull(sample)

    # generating summary values and adding to list
    curr_list <- list()
    curr_list[["Min"]] <- min(curr_vec)
    curr_list[["1stQ"]] <- quantile(curr_vec, 0.25)
    curr_list[["Med"]] <- median(curr_vec)
    curr_list[["Mean"]] <- mean(curr_vec)
    curr_list[["3rdQ"]] <- quantile(curr_vec, 0.75)
    curr_list[["Max"]] <- max(curr_vec)
    curr_list[["SD"]] <- sd(curr_vec)
    curr_list[["CoV"]] <- curr_list[["SD"]] / curr_list[["Mean"]]

    # adding to building table
    SCG_cov_summary_tab <- SCG_cov_summary_tab %>% add_column("{sample}" := unlist(curr_list))
}



# transposing and moving rownames to column
SCG_cov_summary_tab <- SCG_cov_summary_tab %>% t() %>% data.frame(check.names = FALSE)
SCG_cov_summary_tab <- SCG_cov_summary_tab %>% rownames_to_column("Sample")

# sorting by CoV
SCG_cov_summary_tab <- SCG_cov_summary_tab %>% arrange(desc(CoV))

#########
SCG_cov_summary_tab


### let's randomly pick the same number of KOs and look at their CoV (then maybe do it 1,000 times? /shrug)

all_KOs_in_data <- all_cov$KO_ID %>% as.vector


random_KOs <- sample(all_KOs_in_data, length(target_SCGs), replace = FALSE)
random_KOs_cov <- all_cov %>% filter(KO_ID %in% random_KOs) %>% select(-2)
random_KOs_cov_summary_tab <- data.frame(row.names = c("Min", "1stQ", "Med", "Mean", "3rdQ", "Max", "SD", "CoV"))

for ( sample in all_samples ) {

    # getting current sample SCG coverages
    curr_vec <- random_KOs_cov %>% pull(sample)

    # generating summary values and adding to list
    curr_list <- list()
    curr_list[["Min"]] <- min(curr_vec)
    curr_list[["1stQ"]] <- quantile(curr_vec, 0.25)
    curr_list[["Med"]] <- median(curr_vec)
    curr_list[["Mean"]] <- mean(curr_vec)
    curr_list[["3rdQ"]] <- quantile(curr_vec, 0.75)
    curr_list[["Max"]] <- max(curr_vec)
    curr_list[["SD"]] <- sd(curr_vec)
    curr_list[["CoV"]] <- curr_list[["SD"]] / curr_list[["Mean"]]

    # adding to building table
    random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% add_column("{sample}" := unlist(curr_list))
}



# transposing and moving rownames to column
random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% t() %>% data.frame(check.names = FALSE)
random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% rownames_to_column("Sample")

# sorting by CoV
# random_KOs_cov_summary_tab %>% arrange(desc(CoV))

# t1 <- summary(random_KOs_cov_summary_tab$CoV)
t2 <- summary(random_KOs_cov_summary_tab$CoV)

#####################
# FINAL NORMALIZATION PROCESS STARTS HERE
#####################

# MG-DCA
#     - 31 ribosomal proteins as SCGs from KEGG: https://www.genome.jp/kegg/annotation/br01610.html
#     - for all samples in a dataset, for each sample, get coverage of those 31
#     - get list of those that are in interquartile range of those coverages in most or all samples
#     - get counts of how many samples each candidate SCG is in the IQR of coverages of the 31 candidates
#     - take the SCGs that are within the IQR in at least half of the total samples
#     - get the median coverage of those for each sample, and divide the values of all by that
#     - maybe multiply to scale it to something more intuitive

head(SCGs_cov)
F4N1_vec <- SCGs_cov$F4N1
first_Q <- quantile(F4N1_vec, 0.25)
third_Q <- quantile(F4N1_vec, 0.75)

SCGs_cov %>% filter(F4N1 >= first_Q, F4N1 <= third_Q) %>% pull(KO_ID) %>% as.vector


## make table that says in how many samples each KO term is in the IQR of coverage of the candidate SCGs, e.g.
# KO_ID     num_IQR
# K02863    2
# K02864    0
# K02867    8


for ( sample in all_samples ) {

    curr_tab <- SCGs_cov %>% select(KO_ID, all_of(sample))
    first_Q <- quantile(curr_tab[,2], 0.25)
    third_Q <- quantile(curr_tab[,2], 0.75)
    curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector
}

curr_tab %>% filter(P72 >= first_Q, P72 <= third_Q)

curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector

median(curr_tab[,2])

curr_sub_tab <- curr_tab %>% filter(KO_ID %in% curr_IQR_KOs)
median(curr_sub_tab[,2])

### doesn't change it for a single sample, but if we do it for all, and then just use those SCGs that are within the IQR for all, it might, let's see

### run all together ###
# making initial vector
IQR_KO_counts <- vector()
for ( KO in target_SCGs ) {
    IQR_KO_counts[KO] <- 0
}
# going through samples, getting KOs that are in IQR of coverage distribution for the 31 total, and adding a count if a given KO is in the IQR
for ( sample in all_samples ) {

    curr_tab <- SCGs_cov %>% select(KO_ID, all_of(sample))
    first_Q <- quantile(curr_tab[,2], 0.25)
    third_Q <- quantile(curr_tab[,2], 0.75)
    curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector

    # incrementing those in this list
    for ( KO in curr_IQR_KOs) {
        IQR_KO_counts[KO] <- IQR_KO_counts[KO] + 1
    }

}
###                   ###

IQR_KO_counts
# these are pretty low...
# K02906 K02892 K02886 K02890 K02904 K02874 K02895 K02931 K02933 K02881 K02907 K02876 K02871 K02864 K02863 K02867 K02946 K02965 K02982 K02961 K02954
#     32     40     36     42     29     42     56     47     65     60      3     51     39     37     32     40     27     46     42     55     18
# K02994 K02988 K02952 K02948 K02986 K02996 K02992 K02950 K02967 K02956
#    47     50     33     44     33     41     45     40     11     32

# say we want to keep what is in the IQR in at least half of the samples
length(all_samples)
IQR_KO_counts[IQR_KO_counts >= length(all_samples) / 2]
# K02890 K02874 K02895 K02931 K02933 K02881 K02876 K02965 K02982 K02961 K02994 K02988 K02948 K02996 K02992
#     42     42     56     47     65     60     51     46     42     55     47     50     44     41     45

#    L22    L14    L24     L5     L6    L18    L15    S19     S3    S17     S8     S5    S11     S9     S7


wanted_KOs <- names(IQR_KO_counts[IQR_KO_counts >= length(all_samples) / 2])

sub_SCGs_cov <- SCGs_cov %>% filter(KO_ID %in% wanted_KOs)

# median of all SCGs
apply(SCGs_cov[,-1], 2, median)
# median after subsetting to more stable SCGs
apply(sub_SCGs_cov[,-1], 2, median)

data.frame(apply(SCGs_cov[,-1], 2, median), apply(sub_SCGs_cov[,-1], 2, median))
   # it increases in most, stays the same in some, doesn't decrease in any i don't think

normalization_factors <- apply(sub_SCGs_cov[,-1], 2, median)

head(all_cov)

all_cov_norm <- all_cov
all_cov_norm[, -c(1,2)] <- all_cov_norm[, -c(1,2)] / normalization_factors

all_cov_norm %>% filter(KO_ID == "K02890")

summary(all_cov_norm[, -c(1,2)])
range(all_cov_norm[, -c(1,2)])

   # might want to scale this by 1000 or something, or maybe now normalize to coverage per million (maybe just for visualizations, and not for any stats)

#####################

## what if we do DESeq2's approach based on just the candidate SCGs? (again, getting help from here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)
## so now that we have our wanted ones based on the above
wanted_KOs

## now we do:
# 1. create pseudo-reference for those (row-wise geometric mean)
# 2. get ratio of each sample to the pseudo-reference
# 3. get normalization factor (median value of all ratios for each sample)
# 4. normalize full coverage table based on these normalization factors

# 1.
working_sub_SCGs_cov <- sub_SCGs_cov[, -1]
library(psych)
geomeans <- apply(working_sub_SCGs_cov, 1, geometric.mean)

# 2.
ratios_tab <- t( t(working_sub_SCGs_cov) / geomeans )

# 3.
MR_normalization_factors <- apply(ratios_tab, 2, median)

# 4.
all_cov_partMR_norm <- all_cov
all_cov_partMR_norm[, -c(1,2)] <- all_cov_partMR_norm[, -c(1,2)] / MR_normalization_factors

all_cov_partMR_norm %>% filter(KO_ID == "K02890")

summary(all_cov_partMR_norm[, -c(1,2)])
range(all_cov_partMR_norm[, -c(1,2)])

library(psych)
vec <- c(1,2,3,4,0,5)
geometric.mean(vec[vec>0])
exp(mean(log(vec[vec>0])))

# i don't know which of these is more sensible...

# original way chooses ribosomal proteins based on them having coverages within IQR of 31 candidates rib prots, in at least half the samples, then gets ratio of original coverage to the median of those coverages for each sample
# MR-hybrid way chooses the rib prots the same way, then gets MR for each sample of just those, then gets ratio of original coverage to that median-ratio for each sample
  # intuitively, I like the first approach better

head(all_cov_partMR_norm)
head(all_cov_norm)

# writing out
write.table(all_cov_norm, "R-normalized-tables/Combined-KO-SCG-norm-coverages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_cov_partMR_norm, "R-normalized-tables/Combined-KO-SCG-MR-norm-coverages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## clustering all
all_cov_norm_dist <- vegdist(t(all_cov_norm[, -c(1,2)]))
all_cov_norm_hclust <- hclust(all_cov_norm_dist, method = "ward.D2")
all_cov_norm_dendro <- as.dendrogram(all_cov_norm_hclust)

all_cov_norm_dendro_w_timepoint_colors <- all_cov_norm_dendro
all_cov_norm_dendro_w_treatment_colors <- all_cov_norm_dendro
all_cov_norm_dendro_w_flume_colors <- all_cov_norm_dendro

# coloring labels based on timepoint
all_cov_norm_dendro_timepoint_col_vec <- vector()
for ( label in labels(all_cov_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(timepoint_color) %>% as.vector
    all_cov_norm_dendro_timepoint_col_vec <- c(all_cov_norm_dendro_timepoint_col_vec, curr_col)

}
labels_colors(all_cov_norm_dendro_w_timepoint_colors) <- all_cov_norm_dendro_timepoint_col_vec


# coloring labels based on treatment
all_cov_norm_dendro_treatment_col_vec <- vector()
for ( label in labels(all_cov_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(treatment_color) %>% as.vector
    all_cov_norm_dendro_treatment_col_vec <- c(all_cov_norm_dendro_treatment_col_vec, curr_col)

}
labels_colors(all_cov_norm_dendro_w_treatment_colors) <- all_cov_norm_dendro_treatment_col_vec


# coloring labels based on flume
all_cov_norm_dendro_flume_col_vec <- vector()
for ( label in labels(all_cov_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(flume_color) %>% as.vector
    all_cov_norm_dendro_flume_col_vec <- c(all_cov_norm_dendro_flume_col_vec, curr_col)

}
labels_colors(all_cov_norm_dendro_w_flume_colors) <- all_cov_norm_dendro_flume_col_vec



plot(all_cov_norm_dendro_w_timepoint_colors, main="Clustering based on normalized KO coverages (timepoints colored)", ylab="Bray-Curtis dissimilarity")
plot(all_cov_norm_dendro_w_treatment_colors, main="Clustering based on normalized KO coverages (treatments colored)", ylab="Bray-Curtis dissimilarity")
plot(all_cov_norm_dendro_w_flume_colors, main="Clustering based on normalized KO coverages (flumes colored)", ylab="Bray-Curtis dissimilarity")

  ## across all KOs there is no consistency based on anything (this norm or full-MR done next)

## checking against straight median-ratio
straight_MR_norm_covs <- read.table("All-combined-KO-function-MR-normalized-coverages.tsv", sep = "\t", header = TRUE, quote = "")

all_cov_straight_MR_norm_dist <- vegdist(t(straight_MR_norm_covs[, -c(1,2)]))
all_cov_straight_MR_norm_hclust <- hclust(all_cov_straight_MR_norm_dist, method = "ward.D2")
all_cov_straight_MR_norm_dendro <- as.dendrogram(all_cov_straight_MR_norm_hclust)

all_cov_straight_MR_norm_dendro_w_timepoint_colors <- all_cov_straight_MR_norm_dendro
all_cov_straight_MR_norm_dendro_w_treatment_colors <- all_cov_straight_MR_norm_dendro
all_cov_straight_MR_norm_dendro_w_flume_colors <- all_cov_straight_MR_norm_dendro

# coloring labels based on timepoint
all_cov_straight_MR_norm_dendro_timepoint_col_vec <- vector()
for ( label in labels(all_cov_straight_MR_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(timepoint_color) %>% as.vector
    all_cov_straight_MR_norm_dendro_timepoint_col_vec <- c(all_cov_straight_MR_norm_dendro_timepoint_col_vec, curr_col)

}
labels_colors(all_cov_straight_MR_norm_dendro_w_timepoint_colors) <- all_cov_straight_MR_norm_dendro_timepoint_col_vec


# coloring labels based on treatment
all_cov_straight_MR_norm_dendro_treatment_col_vec <- vector()
for ( label in labels(all_cov_straight_MR_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(treatment_color) %>% as.vector
    all_cov_straight_MR_norm_dendro_treatment_col_vec <- c(all_cov_straight_MR_norm_dendro_treatment_col_vec, curr_col)

}
labels_colors(all_cov_straight_MR_norm_dendro_w_treatment_colors) <- all_cov_straight_MR_norm_dendro_treatment_col_vec


# coloring labels based on flume
all_cov_straight_MR_norm_dendro_flume_col_vec <- vector()
for ( label in labels(all_cov_straight_MR_norm_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(flume_color) %>% as.vector
    all_cov_straight_MR_norm_dendro_flume_col_vec <- c(all_cov_straight_MR_norm_dendro_flume_col_vec, curr_col)

}
labels_colors(all_cov_straight_MR_norm_dendro_w_flume_colors) <- all_cov_straight_MR_norm_dendro_flume_col_vec



plot(all_cov_straight_MR_norm_dendro_w_timepoint_colors, main="Clustering based on MR-normalized KO coverages (timepoints colored)", ylab="Bray-Curtis dissimilarity")
plot(all_cov_straight_MR_norm_dendro_w_treatment_colors, main="Clustering based on MR-normalized KO coverages (treatments colored)", ylab="Bray-Curtis dissimilarity")
plot(all_cov_straight_MR_norm_dendro_w_flume_colors, main="Clustering based on MR-normalized KO coverages (flumes colored)", ylab="Bray-Curtis dissimilarity")

  ## across all KOs there is no consistency (this way or with my new normalization done above)


#############
# JUST PLOT SOMETHING TO START, T3 N-FIXATION VS T3 CONTROL; T3 N-FIXATION VS T0
#############


### pulling out just N genes

  # searching KEGG pathways for Nitrogen
N_pathway_info <- keggFind(database = "pathway", query = "nitrogen")

  # getting "Nitrogen metabolism" pathway ID
N_pathway_ID <- sub(names(N_pathway_info), pattern = "path:", replacement = "")

  # getting info for all KOs associated with "Nitrogen metabolism" pathway
N_pathway_KOs_info <- keggLink("ko", N_pathway_ID)

  # getting KO IDs alone
N_pathway_KOs <- as.character(sub(N_pathway_KOs_info, pattern="ko:", replacement=""))

length(N_pathway_KOs) # 65

Nitrogen_metabolism_KEGG_tab <- make_KO_summary_tab(N_pathway_KOs)

# subsetting our normalized table down to just N-related genes
N_all_cov_norm <- all_cov_norm %>% filter(KO_ID %in% N_pathway_KOs)

library(vegan)
library(dendsort)

##### breakdowns

T3_N_and_control_hclusts <- make_dendros("T3", c("N", "Control"), N_all_cov_norm)


T3_N_and_control_hclusts[[2]]

all_hclusts <- make_dendros("", "", all_cov_norm)

all_hclusts[[2]]

all_hclusts_w_env <- make_dendros("", "", all_cov_norm, include_environmental = TRUE)

all_hclusts_w_env[[2]]




# looking at just N-ammended T1 vs T3
N_ammended_T1_T3_samples <- sample_info_tab %>% filter(timepoint %in% c("T1", "T3"), treatment == "N") %>% pull(sample) %>% as.vector




# looking at just T3, N-ammended vs control
T3_N_wanted_samples <- sample_info_tab %>% filter(timepoint == "T3", treatment %in% c("N", "Control")) %>% pull(sample) %>% as.vector

N_all_cov_norm_T3_only <- N_all_cov_norm %>% select(c(1,2), all_of(T3_N_wanted_samples))

## clustering
N_norm_T3_only_dist <- vegdist(t(N_all_cov_norm_T3_only[, -c(1,2)]))
N_norm_T3_only_hclust <- hclust(N_norm_T3_only_dist, method = "ward.D2")
N_norm_T3_only_dendro <- as.dendrogram(N_norm_T3_only_hclust)

N_norm_T3_only_dendro_w_treatment_colors <- N_norm_T3_only_dendro
N_norm_T3_only_dendro_w_flume_colors <- N_norm_T3_only_dendro

# coloring labels based on treatment
N_norm_T3_only_dendro_treatment_col_vec <- vector()
for ( label in labels(N_norm_T3_only_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(treatment_color) %>% as.vector
    N_norm_T3_only_dendro_treatment_col_vec <- c(N_norm_T3_only_dendro_treatment_col_vec, curr_col)

}
labels_colors(N_norm_T3_only_dendro_w_treatment_colors) <- N_norm_T3_only_dendro_treatment_col_vec


# coloring labels based on flume
N_norm_T3_only_dendro_flume_col_vec <- vector()
for ( label in labels(N_norm_T3_only_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(flume_color) %>% as.vector
    N_norm_T3_only_dendro_flume_col_vec <- c(N_norm_T3_only_dendro_flume_col_vec, curr_col)

}
labels_colors(N_norm_T3_only_dendro_w_flume_colors) <- N_norm_T3_only_dendro_flume_col_vec


plot(N_norm_T3_only_dendro_w_treatment_colors, main="Clustering based on normalized N-related coverages (T3, treatments colored)", ylab="Bray-Curtis dissimilarity")
plot(N_norm_T3_only_dendro_w_flume_colors, main="Clustering based on normalized N-related coverages (T3, flumes colored)", ylab="Bray-Curtis dissimilarity")





# looking at T2
# subsetting down to just T2, N-ammended vs control
T2_N_wanted_samples <- sample_info_tab %>% filter(timepoint == "T2", treatment %in% c("N", "Control")) %>% pull(sample) %>% as.vector

N_all_cov_norm_T2_only <- N_all_cov_norm %>% select(c(1,2), all_of(T2_N_wanted_samples))

## clustering
N_norm_T2_only_dist <- vegdist(t(N_all_cov_norm_T2_only[, -c(1,2)]))
N_norm_T2_only_hclust <- hclust(N_norm_T2_only_dist, method = "ward.D2")
N_norm_T2_only_dendro <- as.dendrogram(N_norm_T2_only_hclust)

# coloring labels based on treatment
N_norm_T2_only_dendro_col_vec <- vector()
for ( label in labels(N_norm_T2_only_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(treatment_color) %>% as.vector
    N_norm_T2_only_dendro_col_vec <- c(N_norm_T2_only_dendro_col_vec, curr_col)

}

labels_colors(N_norm_T2_only_dendro) <- N_norm_T2_only_dendro_col_vec

plot(N_norm_T2_only_dendro, main="Clustering based on normalized N-related coverages (T2)", ylab="Bray-Curtis dissimilarity")

# looking at T1
# subsetting down to just T1, N-ammended vs control
T1_N_wanted_samples <- sample_info_tab %>% filter(timepoint == "T1", treatment %in% c("N", "Control")) %>% pull(sample) %>% as.vector

N_all_cov_norm_T1_only <- N_all_cov_norm %>% select(c(1,2), all_of(T1_N_wanted_samples))

## clustering
N_norm_T1_only_dist <- vegdist(t(N_all_cov_norm_T1_only[, -c(1,2)]))
N_norm_T1_only_hclust <- hclust(N_norm_T1_only_dist, method = "ward.D2")
N_norm_T1_only_dendro <- as.dendrogram(N_norm_T1_only_hclust)

# coloring labels based on treatment
N_norm_T1_only_dendro_col_vec <- vector()
for ( label in labels(N_norm_T1_only_dendro) ) {

    curr_col <- sample_info_tab %>% filter(sample == label) %>% pull(treatment_color) %>% as.vector
    N_norm_T1_only_dendro_col_vec <- c(N_norm_T1_only_dendro_col_vec, curr_col)

}

labels_colors(N_norm_T1_only_dendro) <- N_norm_T1_only_dendro_col_vec

plot(N_norm_T1_only_dendro, main="Clustering based on normalized N-related coverages (T1)", ylab="Bray-Curtis dissimilarity")


######### STATS (generally following here: https://www.datanovia.com/en/lessons/anova-in-r/#check-assumptions)

library(ggpubr)
library(rstatix)

packageVersion("rstatix")
packageVersion("tidyverse") # 1.3.0

# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
N_treatment_N_KOs_cov_subset_tab <- make_subset_cov_tab("", "N", N_pathway_KOs)
# getting long-form
N_treatment_N_KOs_cov_subset_tab_long <- make_longform_tab(N_treatment_N_KOs_cov_subset_tab)



## checking for normality by KO_ID (across all timepoints)
N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% shapiro_test(cov)
    # some fail (p < 0.05), but the site notes the shapiro test can be unreliable misleading with > 50 samples and says qqplot is preferred
ggqqplot(N_treatment_N_KOs_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # hm, across all KOs it looks fine, points all fall on the line, seems we can treat these as normally distributed

## checking for heteroskedasticity
N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ timepoint) %>% filter(p <= 0.05)
    # just one doesn't pass, K00371
make_single_KO_single_treatment_scatter_plot("K00371", "", "N")

## running regular anova

N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ timepoint) %>% data.frame(check.names = FALSE)
N_treatment_N_KOs_anova_tab <- left_join(N_treatment_N_KOs_anova_tab, N_treatment_N_KOs_cov_subset_tab %>% select(1,2)) %>% relocate(KO_function, .after = KO_ID)

anova_sig_N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_anova_tab %>% filter(p <= 0.1)

anova_sig_N_treatment_N_KOs <- anova_sig_N_treatment_N_KOs_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig
# make_multi_KO_multi_treatment_scatter_plot(anova_sig_N_treatment_N_KOs)


# making subset table and running posthoc tukey
N_treatment_N_KOs_cov_subset_tab_long_sig <- N_treatment_N_KOs_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_N_treatment_N_KOs)

N_treatment_N_KOs_posthoc_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% tukey_hsd(cov ~ timepoint) %>% data.frame(check.names = FALSE)

posthoc_sig_N_teratment_N_KOs_posthoc_tab <- N_treatment_N_KOs_posthoc_tab %>% filter(p.adj.signif <= 0.1)

posthoc_sig_N_treatment_N_KOs <- posthoc_sig_N_teratment_N_KOs_posthoc_tab %>% pull(KO_ID) %>% as.vector()

# plotting all with sig
make_KO_treatment_scatterplot(posthoc_sig_N_treatment_N_KOs, title = "Those with sig. diff. in N-treatment across time")
make_KO_treatment_scatterplot("K00371")

make_KO_treatment_scatterplot

make_single_KO_single_treatment_scatter_plot("K00371", "N")






make_single_KO_single_treatment_scatter_plot("K00371", "", "N")
make_single_KO_single_treatment_scatter_plot("K00371", "", "P")
make_single_KO_single_treatment_scatter_plot("K00371", "", "Control")

make_single_KO_single_treatment_scatter_plot("K02588", "", "N")
make_single_KO_single_treatment_scatter_plot("K02588", "", "P")
make_single_KO_single_treatment_scatter_plot("K02588", "", "Control")

make_single_KO_single_timepoint_scatter_plot("K02588", "T1", "")
make_single_KO_single_timepoint_scatter_plot("K02588", "T2", "")
make_single_KO_single_timepoint_scatter_plot("K02588", "T3", "")



save.image("explore.RData")


######### HELPER FUNCTIONS ##########

# FUNCTION TO SUBSET A COVERAGE TABLE
make_subset_cov_tab <- function(wanted_timepoints, wanted_treatments, wanted_KOs, input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab, include_environmental = FALSE) {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }
    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments)) %>% pull(sample) %>% as.vector

    # adding environmental samples if specified
    if ( include_environmental ) {

        curr_tab <- input_cov_tab %>% select(1, 2, all_of(c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))))

    } else {

        curr_tab <- input_cov_tab %>% select(1, 2, all_of(curr_wanted_samples))

    }

    # filtering down to wanted KOs if any were provided
    if ( length(wanted_KOs) == 1 ) {

        if ( wanted_KOs != "" ) {

            curr_tab <- curr_tab %>% filter(KO_ID %in% wanted_KOs)

        }

    } else {

        curr_tab <- curr_tab %>% filter(KO_ID %in% wanted_KOs)

    }

    return(curr_tab)

}

# FUNCTION TO MAKE LONG-FORM TABLE FOR STATS

make_longform_tab <- function(input_tab, sam_info_tab = sample_info_tab) {

    curr_long_tab <- input_tab %>% select(-2) %>% pivot_longer(!KO_ID, names_to = "sample", values_to = "cov")

    curr_merged_tab <- left_join(curr_long_tab, sam_info_tab)

    return(curr_merged_tab)
}

make_longform_tab(N_treatment_N_KOs_cov_subset_tab)


# FUNCTION TO MAKE SINGLE-KO, SINGLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_single_KO_single_treatment_scatter_plot <- function(wanted_KO, wanted_treatment, wanted_timepoints = "", input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoints, wanted_treatment, wanted_KO, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov)) + xlab("Timepoint") +
        ylab("Norm. Coverage (SCG)") + ggtitle(wanted_KO, subtitle = paste0("Treatment: ", wanted_treatment)) + theme_bw()

    return(curr_plot)
}


# FUNCTION TO MAKE SINGLE-KO, SINGLE-TIMEPOINT SCATTERPLOTS BASED ON TREATMENT
make_single_KO_single_timepoint_scatter_plot <- function(wanted_KO, wanted_timepoint, wanted_treatments = "", input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoint, wanted_treatments, wanted_KO, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = treatment, y = cov)) + xlab("Treatment") +
        ylab("Norm. Coverage (SCG)") + ggtitle(wanted_KO, subtitle = paste0("Timepoint: ", wanted_timepoint)) + theme_bw()

    return(curr_plot)
}


# FUNCTION TO MAKE MULTIPLE-KO, SINGLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_multi_KO_single_treatment_scatter_plot <- function(wanted_KOs, wanted_treatment, wanted_timepoints = "", input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoints, wanted_treatment, wanted_KOs, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov)) + xlab("Timepoint") +
        ylab("Norm. Coverage (SCG)") + ggtitle(paste0("Treatment: ", wanted_treatment)) + facet_wrap(~KO_ID, scales = "free_y") +
        theme_bw()

    return(curr_plot)
}

make_multiple_KO_single_treatment_scatter_plot(c("K00371", "K02568"), "N")

# FUNCTION TO MAKE MULTIPLE-KO, MULTIPLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_KO_treatment_scatterplot <- function(wanted_KOs, wanted_treatments = "", wanted_timepoints = "", title = "", color_by_flumes = TRUE, input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoints, wanted_treatments, wanted_KOs, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    # getting KO function for title if there is only one being plotted
    if ( length(wanted_KOs) == 1 ) {

        curr_KO_function <- curr_sub_tab %>% filter(KO_ID %in% wanted_KOs) %>% pull(KO_function) %>% as.vector()

    }

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_depth_colors_vec <- unique(curr_sub_tab_long$depth_color[order(curr_sub_tab_long$depth)]) %>% as.vector()
    curr_flume_colors_vec <- unique(curr_sub_tab_long$flume_color[order(curr_sub_tab_long$flume)]) %>% as.vector()

    if ( title == "" ) {

        if ( color_by_flumes == FALSE ) {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = depth), size = 3) +
                scale_color_manual(values = curr_depth_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (SCG)") + facet_grid(treatment ~ KO_ID, scales = "free_y") +
                theme_bw() + theme(legend.position = "bottom")

        } else {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = flume), size = 3) +
                scale_color_manual(values = curr_flume_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (SCG)") + facet_grid(treatment ~ KO_ID, scales = "free_y") +
                theme_bw() + theme(legend.position = "bottom")

        }

    } else {

        if ( color_by_flumes == FALSE ) {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = depth), size = 3) +
                scale_color_manual(values = curr_depth_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (SCG)") + facet_grid(treatment ~ KO_ID, scales = "free_y") + ggtitle(title) +
                theme_bw() + theme(legend.position = "bottom")

        } else {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = flume), size = 3) +
                scale_color_manual(values = curr_flume_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (SCG)") + facet_grid(treatment ~ KO_ID, scales = "free_y") + ggtitle(title) +
                theme_bw() + theme(legend.position = "bottom")

        }

    }

    return(curr_plot)
}

make_KO_treatment_scatterplot(c("K00371", "K02568"))
make_KO_treatment_scatterplot("K00371")
make_KO_treatment_scatterplot("K00371", color_by_flumes = FALSE)
make_KO_treatment_scatterplot("K02933")

# FUNCTION TO MAKE DENDROGRAMS
make_dendros <- function(wanted_timepoints, wanted_treatments, input_cov_tab, sam_info_tab = sample_info_tab, include_environmental = FALSE) {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }
    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments)) %>% pull(sample) %>% as.vector

    # adding environmental samples if specified
    if ( include_environmental ) {

        curr_tab <- input_cov_tab %>% select(all_of(c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))))

    } else {

        curr_tab <- input_cov_tab %>% select(all_of(curr_wanted_samples))

    }
    # making dendrogram
    curr_dist <- vegdist(t(curr_tab))
    curr_hclust <- hclust(curr_dist, method = "ward.D2")
    curr_dendro <- as.dendrogram(curr_hclust)

    # making colors
    curr_dendro_w_timepoint_cols <- curr_dendro
    curr_dendro_w_treatment_cols <- curr_dendro
    curr_dendro_w_flume_cols <- curr_dendro

    curr_timepoint_cols_vec <- vector()
    curr_treatment_cols_vec <- vector()
    curr_flume_cols_vec <- vector()

    for ( curr_label in labels(curr_dendro) ) {

        curr_timepoint_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(timepoint_color) %>% as.vector
        curr_timepoint_cols_vec <- c(curr_timepoint_cols_vec, curr_timepoint_col)

        curr_treatment_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(treatment_color) %>% as.vector
        curr_treatment_cols_vec <- c(curr_treatment_cols_vec, curr_treatment_col)

        curr_flume_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(flume_color) %>% as.vector
        curr_flume_cols_vec <- c(curr_flume_cols_vec, curr_flume_col)

    }

    labels_colors(curr_dendro_w_timepoint_cols) <- curr_timepoint_cols_vec
    labels_colors(curr_dendro_w_treatment_cols) <- curr_treatment_cols_vec
    labels_colors(curr_dendro_w_flume_cols) <- curr_flume_cols_vec

    # making base title
    curr_base_title <- paste0("Timepoint(s): ", paste0(wanted_timepoints, collapse = ","), "  |  Treatment(s): ", paste0(wanted_treatments, collapse = ", "), "\n")

    plot(curr_dendro_w_timepoint_cols, main = paste(curr_base_title, "(colored by timepoint)", sep = " "), ylab="Bray-Curtis dissimilarity")
    curr_timepoint_plot <- recordPlot()
    plot(curr_dendro_w_treatment_cols, main = paste(curr_base_title, "(colored by treatment)", sep = " "), ylab="Bray-Curtis dissimilarity")
    curr_treatment_plot <- recordPlot()
    plot(curr_dendro_w_flume_cols, main = paste(curr_base_title, "(colored by flume)", sep = " "), ylab="Bray-Curtis dissimilarity")
    curr_flume_plot <- recordPlot()

    return(list(timepoint_colored_plot = curr_timepoint_plot, treatment_colored_plot = curr_treatment_plot, flume_colored_plot = curr_flume_plot))

}

# FUNCTION TO MAKE KO SUMMARY TABLE
make_KO_summary_tab <- function(target_KOs) {

    # removing "Not annotated" if provided, will be added as needed below
    target_KOs <- target_KOs[target_KOs != "Not annotated"]

    # the kegg API limits individual requests to 10, so breaking the total we want into blocks of 10
    list_of_pathway_KO_blocks <- split(target_KOs, ceiling(seq_along(target_KOs)/10))
    num_blocks <- length(list_of_pathway_KO_blocks)

    # initializing some vectors we're going to populate
    entry_vec <- vector()
    name_vec <- vector()
    def_vec <- vector()
    module_ids_vec <- vector()
    module_vec <- vector()
    path_ids_vec <- vector()
    path_vec <- vector()

    # iterating through each block of at most 10 KO terms
    for ( block in seq(1, num_blocks) ) {

        # getting and storing the current block of KO terms' information
        current <- keggGet(list_of_pathway_KO_blocks[[block]])

        # here iterating through that stored block of information, one term at a time to get the info we want for each
        for ( num in seq(1, length(current)) ) {

            # KO ID
            current_entry <- current[[num]]$ENTRY
            if ( length(current_entry) > 0 ) {
                entry_vec <- c(entry_vec, current_entry)
            } else {
                entry_vec <- c(entry_vec, NA)
            }

            # KO name(s)
            current_name <- current[[num]]$NAME
            if ( length(current_name) > 0 ) {
                name_vec <- c(name_vec, current_name)
            } else {
                name_vec <- c(name_vec, NA)
            }

            # KO definition
            current_def <- current[[num]]$DEFINITION
            if ( length(current_def) > 0 ) {
                # getting rid of EC number if present
                current_def <- sub(current_def, pattern=" \\[EC:.*$", replacement="")
                def_vec <- c(def_vec, current_def)
            } else {
                def_vec <- c(def_vec, NA)
            }

            # module ID(s)
            current_module_id <- stri_join(names(current[[num]]$MODULE), collapse="; ")
            if ( length(current_module_id) > 0 ) {
                module_ids_vec <- c(module_ids_vec, current_module_id)
            } else {
                module_ids_vec <- c(module_ids_vec, NA)
            }

            # module definition(s)
            current_module <- stri_join(current[[num]]$MODULE, collapse="; ")
            if ( length(current_module) > 0 ) {
                module_vec <- c(module_vec, current_module)
            } else {
                module_vec <- c(module_vec, NA)
            }

            # pathway ID(s)
            current_path_id <- stri_join(names(current[[num]]$PATHWAY), collapse="; ")
            if ( length(current_path_id) > 0 ) {
                path_ids_vec <- c(path_ids_vec, current_path_id)
            } else {
                path_ids_vec <- c(path_ids_vec, NA)
            }

            # pathway definition(s)
            current_path <- stri_join(current[[num]]$PATHWAY, collapse="; ")
            if ( length(current_path) > 0 ) {
                path_vec <- c(path_vec, current_path)
            } else {
                path_vec <- c(path_vec, NA)
            }

        }
    }

    # now combining into table
    out_tab <- data.frame("KO_ID"=entry_vec, "KO_name"=name_vec, "KO_def"=def_vec, "module_IDs"=module_ids_vec, "module_defs"=module_vec, "pathway_IDs"=path_ids_vec, "pathway_defs"=path_vec, stringsAsFactors=F)

    # some KO terms may have been removed from KEGG since annotation was performed, and therefore wouldn't be found
    # putting a check and reporting if that's the case for any, adding to table to keep track of them, but making all values set to 'Not found at KEGG'
    got_KOs <- out_tab %>% pull(KO_ID)

    missed_KOs <- setdiff(target_KOs, got_KOs)

    if ( length(missed_KOs) != 0 ) {
        num_missed <- length(missed_KOs)

        cat(" ", num_missed, "KOs were not found at KEGG, they were likely removed after the utilized annotation db was created:\n\n")

        for ( term in missed_KOs ) {
            cat("\t", term, "\n")
            out_tab <- rbind(out_tab, c(term, rep("Not found at KEGG", 6)))
        }

        cat("\n")
        cat("  These have been added to the output table with 'Not found at KEGG' in all other fields besides the 'KO_ID' one.\n\n")

    }

    # sorting out table
    out_tab <- out_tab %>% arrange(KO_ID)

    return(out_tab)

}












