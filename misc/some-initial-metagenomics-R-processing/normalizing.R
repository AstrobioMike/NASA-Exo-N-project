### script generating ribosomal-protein normalized coverages
  ## initial exploring at bottom, final at top here
library(tidyverse)

all_cov <- read.table("../combined-outputs/All-combined-KO-function-coverages.tsv", sep = "\t", header = TRUE, quote = "")
head(all_cov)
dim(all_cov)

# reading in sample info table
sample_info_tab <- read.table("../sample-and-group-info.tsv", sep = "\t", header = TRUE, comment.char = "")


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

SCGs_cov <- all_cov %>% filter(KO_ID %in% target_SCGs) %>% select(-2)

all_samples <- names(SCGs_cov)[-1]

#####################
# FINAL NORMALIZATION PROCESS STARTS HERE
#####################

# my MG-DCA approach
#     - 31 ribosomal proteins as SCGs from KEGG: https://www.genome.jp/kegg/annotation/br01610.html
#     - for all samples in a dataset, for each sample, get coverage of those 31
#     - get list of those that are in interquartile range of those coverages in most or all samples
#     - get counts of how many samples for which each candidate SCG is in the IQR of coverages of the 31 candidates
#     - take the SCGs that are within the IQR in at least half of the total samples
#     - get the median coverage of those for each sample, and divide the values of all in that sample by that sample's median coverage of the target rib prots
#     - maybe multiply to scale it to something more intuitive


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
head(all_cov_norm)

write.table(all_cov_norm, "Combined-KO-rib-prot-norm-coverages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)







##############################################################
####################### all initial code #####################
##############################################################

# all_cov <- read.table("../combined-outputs/All-combined-KO-function-coverages.tsv", sep = "\t", header = TRUE, quote = "")
# head(all_cov)
# dim(all_cov)
#
# # reading in sample info table
# sample_info_tab <- read.table("../sample-and-group-info.tsv", sep = "\t", header = TRUE, comment.char = "")
#
#
# # choosing based on archaea and bacteria from KEGG's page here: https://www.genome.jp/kegg/annotation/br01610.html
# target_SCGs <- c("L3" = "K02906",
#                  "L23" = "K02892",
#                  "L2" = "K02886",
#                  "L22" = "K02890",
#                  "L29" = "K02904",
#                  "L14" = "K02874",
#                  "L24" = "K02895",
#                  "L5" = "K02931",
#                  "L6" = "K02933",
#                  "L18" = "K02881",
#                  "L30" = "K02907",
#                  "L15" = "K02876",
#                  "L13" = "K02871",
#                  "L10" = "K02864",
#                  "L1" = "K02863",
#                  "L11" = "K02867",
#                  "S10" = "K02946",
#                  "S19" = "K02965",
#                  "S3" = "K02982",
#                  "S17" = "K02961",
#                  "S14" = "K02954",
#                  "S8" = "K02994",
#                  "S5" = "K02988",
#                  "S13" = "K02952",
#                  "S11" = "K02948",
#                  "S4" = "K02986",
#                  "S9" = "K02996",
#                  "S7" = "K02992",
#                  "S12" = "K02950",
#                  "S2" = "K02967",
#                  "S15" = "K02956")
#
# length(target_SCGs)
#
#
#
#
#
#
# # median of all
# apply(all_cov[,c(3:dim(all_cov)[2])], 2, median)
#
# # median of target SCGs
# SCGs_cov <- all_cov %>% filter(KO_ID %in% target_SCGs)
# apply(SCGs_cov[,-1], 2, median)
#
# apply(SCGs_cov[,-1], 2, summary)
#
#
# # median of individual target SCGs across samples (though this would be directly affected by library size, maybe the distributions would look similar for good candidate final SCGs to use?)
#
# plot_hist_of_SCG_covs <- function(target_SCG, cov_tab = all_cov) {
#
#     curr_covs <- cov_tab %>% filter(KO_ID == target_SCG) %>% select(-c(KO_ID, KO_function)) %>% as.numeric
#     curr_summary <- summary(curr_covs)
#     print(curr_summary)
#
#     return(hist(curr_covs, main = target_SCG, xlab = "Coverages"))
# }
#
#    # those with the same number of "#"s have similar distributions
# plot_hist_of_SCG_covs("K02906") #
# plot_hist_of_SCG_covs("K02892") #
# plot_hist_of_SCG_covs("K02886") #
# plot_hist_of_SCG_covs("K02890") #
# plot_hist_of_SCG_covs("K02904") #
# plot_hist_of_SCG_covs("K02874") ##
# plot_hist_of_SCG_covs("K02895") ##
# plot_hist_of_SCG_covs("K02931") ##
# plot_hist_of_SCG_covs("K02881") #
# plot_hist_of_SCG_covs("K02907") #
# plot_hist_of_SCG_covs("K02876") ##
# plot_hist_of_SCG_covs("K02871") #
# plot_hist_of_SCG_covs("K02864") #
# plot_hist_of_SCG_covs("K02863") ##
# plot_hist_of_SCG_covs("K02867") #
#
# plot_hist_of_SCG_covs("K02946") #
# plot_hist_of_SCG_covs("K02965") #
# plot_hist_of_SCG_covs("K02982") ##
# plot_hist_of_SCG_covs("K02961") ##
# plot_hist_of_SCG_covs("K02954") ###
# plot_hist_of_SCG_covs("K02994") #
# plot_hist_of_SCG_covs("K02988") #
# plot_hist_of_SCG_covs("K02952") #
# plot_hist_of_SCG_covs("K02948") ##
# plot_hist_of_SCG_covs("K02986") #
# plot_hist_of_SCG_covs("K02996") #
# plot_hist_of_SCG_covs("K02992") #
# plot_hist_of_SCG_covs("K02950") #
# plot_hist_of_SCG_covs("K02967") #
# plot_hist_of_SCG_covs("K02956") #
#
#
# SCGs_cov <- all_cov %>% filter(KO_ID %in% target_SCGs) %>% select(-2)
#
# # want to see histograms of all SCGs in a single sample?
# hist(SCGs_cov %>% pull(F4N1))
# hist(SCGs_cov %>% pull(F4N5))
#
# summary(SCGs_cov %>% pull(F4N1))
# summary(SCGs_cov %>% pull(F4N5))
#
#
#
# # summaries of all
# apply(SCGs_cov[, -1], 2, summary)
#
# SCGs_cov_long <- SCGs_cov %>% pivot_longer(!KO_ID, names_to = "Sample", values_to = "coverage")
#
# ggplot(SCGs_cov_long, aes(coverage)) + geom_histogram() + facet_wrap(~Sample, scales = "free")
#
#
# ## want some summary stats including coefficient of variation of all SCGs in each individual sample
# ## thinking something like
# # stat  F4N1    F4N5
# # Min
# # 1st Q
# # Med
# # Mean
# # 3rd Q
# # Max
# # SD
# # CoV
#
# ## Or maybe transposed
# # Sample    Min     1stQ    Med     Mean    3rdQ    Max     SD      CoV
# # F4N1
# # F4N5
#
# ##########
# all_samples <- names(SCGs_cov)[-1]
#
# SCG_cov_summary_tab <- data.frame(row.names = c("Min", "1stQ", "Med", "Mean", "3rdQ", "Max", "SD", "CoV"))
#
# for ( sample in all_samples ) {
#
#     # getting current sample SCG coverages
#     curr_vec <- SCGs_cov %>% pull(sample)
#
#     # generating summary values and adding to list
#     curr_list <- list()
#     curr_list[["Min"]] <- min(curr_vec)
#     curr_list[["1stQ"]] <- quantile(curr_vec, 0.25)
#     curr_list[["Med"]] <- median(curr_vec)
#     curr_list[["Mean"]] <- mean(curr_vec)
#     curr_list[["3rdQ"]] <- quantile(curr_vec, 0.75)
#     curr_list[["Max"]] <- max(curr_vec)
#     curr_list[["SD"]] <- sd(curr_vec)
#     curr_list[["CoV"]] <- curr_list[["SD"]] / curr_list[["Mean"]]
#
#     # adding to building table
#     SCG_cov_summary_tab <- SCG_cov_summary_tab %>% add_column("{sample}" := unlist(curr_list))
# }
#
#
#
# # transposing and moving rownames to column
# SCG_cov_summary_tab <- SCG_cov_summary_tab %>% t() %>% data.frame(check.names = FALSE)
# SCG_cov_summary_tab <- SCG_cov_summary_tab %>% rownames_to_column("Sample")
#
# # sorting by CoV
# SCG_cov_summary_tab <- SCG_cov_summary_tab %>% arrange(desc(CoV))
#
# #########
# SCG_cov_summary_tab
#
#
# ### let's randomly pick the same number of KOs and look at their CoV (then maybe do it 1,000 times? /shrug)
#
# all_KOs_in_data <- all_cov$KO_ID %>% as.vector
#
#
# random_KOs <- sample(all_KOs_in_data, length(target_SCGs), replace = FALSE)
# random_KOs_cov <- all_cov %>% filter(KO_ID %in% random_KOs) %>% select(-2)
# random_KOs_cov_summary_tab <- data.frame(row.names = c("Min", "1stQ", "Med", "Mean", "3rdQ", "Max", "SD", "CoV"))
#
# for ( sample in all_samples ) {
#
#     # getting current sample SCG coverages
#     curr_vec <- random_KOs_cov %>% pull(sample)
#
#     # generating summary values and adding to list
#     curr_list <- list()
#     curr_list[["Min"]] <- min(curr_vec)
#     curr_list[["1stQ"]] <- quantile(curr_vec, 0.25)
#     curr_list[["Med"]] <- median(curr_vec)
#     curr_list[["Mean"]] <- mean(curr_vec)
#     curr_list[["3rdQ"]] <- quantile(curr_vec, 0.75)
#     curr_list[["Max"]] <- max(curr_vec)
#     curr_list[["SD"]] <- sd(curr_vec)
#     curr_list[["CoV"]] <- curr_list[["SD"]] / curr_list[["Mean"]]
#
#     # adding to building table
#     random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% add_column("{sample}" := unlist(curr_list))
# }
#
#
#
# # transposing and moving rownames to column
# random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% t() %>% data.frame(check.names = FALSE)
# random_KOs_cov_summary_tab <- random_KOs_cov_summary_tab %>% rownames_to_column("Sample")
#
# # sorting by CoV
# # random_KOs_cov_summary_tab %>% arrange(desc(CoV))
#
# # t1 <- summary(random_KOs_cov_summary_tab$CoV)
# t2 <- summary(random_KOs_cov_summary_tab$CoV)
#
# #####################
# # FINAL NORMALIZATION PROCESS STARTS HERE
# #####################
#
# # MG-DCA
# #     - 31 ribosomal proteins as SCGs from KEGG: https://www.genome.jp/kegg/annotation/br01610.html
# #     - for all samples in a dataset, for each sample, get coverage of those 31
# #     - get list of those that are in interquartile range of those coverages in most or all samples
# #     - get counts of how many samples each candidate SCG is in the IQR of coverages of the 31 candidates
# #     - take the SCGs that are within the IQR in at least half of the total samples
# #     - get the median coverage of those for each sample, and divide the values of all by that
# #     - maybe multiply to scale it to something more intuitive
#
# head(SCGs_cov)
# F4N1_vec <- SCGs_cov$F4N1
# first_Q <- quantile(F4N1_vec, 0.25)
# third_Q <- quantile(F4N1_vec, 0.75)
#
# SCGs_cov %>% filter(F4N1 >= first_Q, F4N1 <= third_Q) %>% pull(KO_ID) %>% as.vector
#
#
# ## make table that says in how many samples each KO term is in the IQR of coverage of the candidate SCGs, e.g.
# # KO_ID     num_IQR
# # K02863    2
# # K02864    0
# # K02867    8
#
#
# for ( sample in all_samples ) {
#
#     curr_tab <- SCGs_cov %>% select(KO_ID, all_of(sample))
#     first_Q <- quantile(curr_tab[,2], 0.25)
#     third_Q <- quantile(curr_tab[,2], 0.75)
#     curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector
# }
#
# curr_tab %>% filter(P72 >= first_Q, P72 <= third_Q)
#
# curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector
#
# median(curr_tab[,2])
#
# curr_sub_tab <- curr_tab %>% filter(KO_ID %in% curr_IQR_KOs)
# median(curr_sub_tab[,2])
#
# ### doesn't change it for a single sample, but if we do it for all, and then just use those SCGs that are within the IQR for all, it might, let's see
#
# ### run all together ###
# # making initial vector
# IQR_KO_counts <- vector()
# for ( KO in target_SCGs ) {
#     IQR_KO_counts[KO] <- 0
# }
# # going through samples, getting KOs that are in IQR of coverage distribution for the 31 total, and adding a count if a given KO is in the IQR
# for ( sample in all_samples ) {
#
#     curr_tab <- SCGs_cov %>% select(KO_ID, all_of(sample))
#     first_Q <- quantile(curr_tab[,2], 0.25)
#     third_Q <- quantile(curr_tab[,2], 0.75)
#     curr_IQR_KOs <- curr_tab[curr_tab[sample] >= first_Q & curr_tab[sample] <= third_Q, 1] %>% as.vector
#
#     # incrementing those in this list
#     for ( KO in curr_IQR_KOs) {
#         IQR_KO_counts[KO] <- IQR_KO_counts[KO] + 1
#     }
#
# }
# ###                   ###
#
# IQR_KO_counts
# # these are pretty low...
# # K02906 K02892 K02886 K02890 K02904 K02874 K02895 K02931 K02933 K02881 K02907 K02876 K02871 K02864 K02863 K02867 K02946 K02965 K02982 K02961 K02954
# #     32     40     36     42     29     42     56     47     65     60      3     51     39     37     32     40     27     46     42     55     18
# # K02994 K02988 K02952 K02948 K02986 K02996 K02992 K02950 K02967 K02956
# #    47     50     33     44     33     41     45     40     11     32
#
# # say we want to keep what is in the IQR in at least half of the samples
# length(all_samples)
# IQR_KO_counts[IQR_KO_counts >= length(all_samples) / 2]
# # K02890 K02874 K02895 K02931 K02933 K02881 K02876 K02965 K02982 K02961 K02994 K02988 K02948 K02996 K02992
# #     42     42     56     47     65     60     51     46     42     55     47     50     44     41     45
#
# #    L22    L14    L24     L5     L6    L18    L15    S19     S3    S17     S8     S5    S11     S9     S7
#
#
# wanted_KOs <- names(IQR_KO_counts[IQR_KO_counts >= length(all_samples) / 2])
#
# sub_SCGs_cov <- SCGs_cov %>% filter(KO_ID %in% wanted_KOs)
#
# # median of all SCGs
# apply(SCGs_cov[,-1], 2, median)
# # median after subsetting to more stable SCGs
# apply(sub_SCGs_cov[,-1], 2, median)
#
# data.frame(apply(SCGs_cov[,-1], 2, median), apply(sub_SCGs_cov[,-1], 2, median))
#    # it increases in most, stays the same in some, doesn't decrease in any i don't think
#
# normalization_factors <- apply(sub_SCGs_cov[,-1], 2, median)
#
# head(all_cov)
#
# all_cov_norm <- all_cov
# all_cov_norm[, -c(1,2)] <- all_cov_norm[, -c(1,2)] / normalization_factors
#
# all_cov_norm %>% filter(KO_ID == "K02890")
#
# summary(all_cov_norm[, -c(1,2)])
# range(all_cov_norm[, -c(1,2)])
#
#    # might want to scale this by 1000 or something, or maybe now normalize to coverage per million (maybe just for visualizations, and not for any stats)
#
# #####################
#
# ## what if we do DESeq2's approach based on just the candidate SCGs? (again, getting help from here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)
# ## so now that we have our wanted ones based on the above
# wanted_KOs
#
# ## now we do:
# # 1. create pseudo-reference for those (row-wise geometric mean)
# # 2. get ratio of each sample to the pseudo-reference
# # 3. get normalization factor (median value of all ratios for each sample)
# # 4. normalize full coverage table based on these normalization factors
#
# # 1.
# working_sub_SCGs_cov <- sub_SCGs_cov[, -1]
# library(psych)
# geomeans <- apply(working_sub_SCGs_cov, 1, geometric.mean)
#
# # 2.
# ratios_tab <- t( t(working_sub_SCGs_cov) / geomeans )
#
# # 3.
# MR_normalization_factors <- apply(ratios_tab, 2, median)
#
# # 4.
# all_cov_partMR_norm <- all_cov
# all_cov_partMR_norm[, -c(1,2)] <- all_cov_partMR_norm[, -c(1,2)] / MR_normalization_factors
#
# all_cov_partMR_norm %>% filter(KO_ID == "K02890")
#
# summary(all_cov_partMR_norm[, -c(1,2)])
# range(all_cov_partMR_norm[, -c(1,2)])
#
# library(psych)
# vec <- c(1,2,3,4,0,5)
# geometric.mean(vec[vec>0])
# exp(mean(log(vec[vec>0])))
#
# # i don't know which of these is more sensible...
#
# # original way chooses ribosomal proteins based on them having coverages within IQR of 31 candidates rib prots, in at least half the samples, then gets ratio of original coverage to the median of those coverages for each sample
# # MR-hybrid way chooses the rib prots the same way, then gets MR for each sample of just those, then gets ratio of original coverage to that median-ratio for each sample
#   # intuitively, I like the first approach better
#
# head(all_cov_partMR_norm)
# head(all_cov_norm)
#
# # writing out
# write.table(all_cov_norm, "R-normalized-tables/Combined-KO-SCG-norm-coverages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(all_cov_partMR_norm, "R-normalized-tables/Combined-KO-SCG-MR-norm-coverages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
