library(tidyverse)
library(rstatix)
library(ggpubr)


# loading my functions
source("my-functions.R")

# reading in rib. prot. normalized table
all_cov_norm <- read.table("Combined-KO-rib-prot-norm-coverages.tsv", sep = "\t", header = TRUE, quote = "")
head(all_cov_norm)

# reading in sample info table
sample_info_tab <- read.table("../sample-and-group-info.tsv", sep = "\t", header = TRUE, comment.char = "", quote = "")
head(sample_info_tab)

# reading in master kegg tab (just in case we want it for anything)
Master_KEGG_tab <- read.table("kegg-tables/Master-kegg-annotations.tsv", sep = "\t", header = TRUE, quote = "", comment.char = "")
dim(Master_KEGG_tab)

all_our_KOs <- Master_KEGG_tab %>% pull(KO_ID) %>% as.vector()

# reading in Nitrogen-focused tab
Nitrogen_metabolism_KEGG_tab <- read.table("kegg-tables/Kegg-N-annotations.tsv", sep = "\t", header = TRUE, quote = "", comment.char = "")
dim(Nitrogen_metabolism_KEGG_tab)

# getting N-related KOs
N_pathway_KOs <- Nitrogen_metabolism_KEGG_tab %>% pull(KO_ID) %>% as.vector()
length(N_pathway_KOs)

######### STATS (generally following here: https://www.datanovia.com/en/lessons/anova-in-r/#check-assumptions)

# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
N_treatment_N_KOs_cov_subset_tab <- make_subset_cov_tab("", "N", N_pathway_KOs)

# getting long-form
N_treatment_N_KOs_cov_subset_tab_long <- make_longform_tab(N_treatment_N_KOs_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)
N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% shapiro_test(cov)
    # some fail (p < 0.05), but the site notes the shapiro test can be unreliable misleading with > 50 samples and says qqplot is preferred
ggqqplot(N_treatment_N_KOs_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # hm, actually across all KOs it looks fine, points virtually all fall on the line, seems we can treat these as normally distributed

## checking for heteroskedasticity
N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ timepoint) %>% filter(p <= 0.05)
    # just one doesn't pass, K00371

# looking at it a few ways
make_KO_treatment_scatterplot("K00371", wanted_treatments = "N")

make_KO_treatment_scatterplot("K00371")

make_single_KO_single_treatment_scatter_plot("K00371", "N")
make_single_KO_single_treatment_scatter_plot("K00371", "Control")

make_single_KO_single_timepoint_scatter_plot("K00371", "T1")
make_single_KO_single_timepoint_scatter_plot("K00371", "T3")



## running regular anova
N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ timepoint) %>% data.frame(check.names = FALSE)
# adding KO annotations
N_treatment_N_KOs_anova_tab <- left_join(N_treatment_N_KOs_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_N_treatment_N_KOs_anova_tab)

anova_sig_N_treatment_N_KOs <- anova_sig_N_treatment_N_KOs_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_N_treatment_N_KOs, color_by_flumes = FALSE)

# making subset table and running posthoc tukey
N_treatment_N_KOs_cov_subset_tab_long_sig <- N_treatment_N_KOs_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_N_treatment_N_KOs)
N_treatment_N_KOs_posthoc_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% tukey_hsd(cov ~ timepoint) %>% data.frame(check.names = FALSE)

# adding KO annotations
N_treatment_N_KOs_posthoc_tab <- left_join(N_treatment_N_KOs_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_N_treatment_N_KOs_posthoc_tab <- N_treatment_N_KOs_posthoc_tab %>% filter(p.adj <= 0.05)

# adding log2fold changes
posthoc_sig_N_treatment_N_KOs_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_N_treatment_N_KOs_posthoc_tab, "N", N_treatment_N_KOs_cov_subset_tab_long)

# getting KO IDs
posthoc_sig_N_treatment_N_KOs <- posthoc_sig_N_treatment_N_KOs_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_N_treatment_N_KOs)

posthoc_sig_N_treatment_N_KOs_posthoc_tab %>% select(1,2,3,5,6,11,13)
#     KO_ID               KO_name                                                    KO_def group1 group2   p.adj log2fc(group2/group1)
# 1  K00371      narH, narY, nxrB  nitrate reductase / nitrite oxidoreductase, beta subunit     T1     T2 0.03940                 -1.63
# 2  K00371      narH, narY, nxrB  nitrate reductase / nitrite oxidoreductase, beta subunit     T1     T3 0.00595                 -3.27
# 3  K00459             ncd2, npd                                   nitronate monooxygenase     T2     T3 0.02430                  1.37
# 4  K02305                  norC                          nitric oxide reductase subunit C     T1     T3 0.03100                  2.08
# 5  K02568                  napB nitrate reductase (cytochrome), electron transfer subunit     T1     T3 0.00152                 -3.45
# 6  K02575 NRT, narK, nrtP, nasA  MFS transporter, NNP family, nitrate/nitrite transporter     T1     T3 0.04420                  1.94
# 7  K02586                  nifD           nitrogenase molybdenum-iron protein alpha chain     T1     T3 0.04710                  1.41
# 8  K02588                  nifH                             nitrogenase iron protein NifH     T2     T3 0.04530                  0.99
# 9  K02591                  nifK            nitrogenase molybdenum-iron protein beta chain     T1     T3 0.04260                  1.40
# 10 K20933                K20933                                hydrazine synthase subunit     T2     T3 0.03090                 -0.88

# looking at them

make_KO_treatment_scatterplot("K00371", title = get_KO_info("K00371"), point_size = 2)

# writing pdf
pdf("figs/K00371-norm-cov-plot.pdf", width = 3.5, height = 6)
make_KO_treatment_scatterplot("K00371", title = get_KO_info("K00371"), point_size = 2)
dev.off()

make_KO_treatment_scatterplot("K02305", title = get_KO_info("K02305"))
make_KO_treatment_scatterplot("K02588", title = get_KO_info("K02588"))

make_KO_treatment_scatterplot(c("K02588", "K02586", "K02591"), point_size = 2, scale_setting = "fixed")
# writing pdf
pdf("figs/K02586-K02588-K02591-norm-cov-plot.pdf", width = 3.5, height = 6)
make_KO_treatment_scatterplot(c("K02588", "K02586", "K02591"), point_size = 2, scale_setting = "fixed")
dev.off()


#######################################################################################################
########## doing the same for the Control to find which above are NOT also sig in one of them #########
#######################################################################################################
# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
Control_treatment_N_KOs_cov_subset_tab <- make_subset_cov_tab("", "Control", N_pathway_KOs)

# getting long-form
Control_treatment_N_KOs_cov_subset_tab_long <- make_longform_tab(Control_treatment_N_KOs_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)
ggqqplot(Control_treatment_N_KOs_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # most look fine here except for maybe 3 (K01915, K02586, K20932)

## checking for heteroskedasticity
Control_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ timepoint) %>% filter(p <= 0.05)
    # just 3 don't pass, K00370, K01674, K03385

## running regular anova
Control_treatment_N_KOs_anova_tab <- Control_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ timepoint) %>% data.frame(check.names = FALSE)
# adding KO annotations
Control_treatment_N_KOs_anova_tab <- left_join(Control_treatment_N_KOs_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_Control_treatment_N_KOs_anova_tab <- Control_treatment_N_KOs_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_Control_treatment_N_KOs_anova_tab) # 5 of them

anova_sig_Control_treatment_N_KOs <- anova_sig_Control_treatment_N_KOs_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_Control_treatment_N_KOs)

# making subset table and running posthoc tukey
Control_treatment_N_KOs_cov_subset_tab_long_sig <- Control_treatment_N_KOs_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_Control_treatment_N_KOs)
Control_treatment_N_KOs_posthoc_tab <- Control_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% tukey_hsd(cov ~ timepoint) %>% data.frame(check.names = FALSE)

# adding KO annotations
Control_treatment_N_KOs_posthoc_tab <- left_join(Control_treatment_N_KOs_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_Control_treatment_N_KOs_posthoc_tab <- Control_treatment_N_KOs_posthoc_tab %>% filter(p.adj.signif <= 0.1)

# adding log2fold changes
posthoc_sig_Control_treatment_N_KOs_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_Control_treatment_N_KOs_posthoc_tab, "Control", Control_treatment_N_KOs_cov_subset_tab_long)

# getting KO IDs
posthoc_sig_Control_treatment_N_KOs <- posthoc_sig_Control_treatment_N_KOs_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_Control_treatment_N_KOs)

posthoc_sig_Control_treatment_N_KOs_posthoc_tab %>% select(1,2,3,5,6,11,13)
#    KO_ID          KO_name                                                    KO_def group1 group2   p.adj log2fc(group2/group1)
# 1 K00370 narG, narZ, nxrA nitrate reductase / nitrite oxidoreductase, alpha subunit     T1     T2 0.04370                 -2.68
# 2 K00370 narG, narZ, nxrA nitrate reductase / nitrite oxidoreductase, alpha subunit     T1     T3 0.02550                 -3.78
# 3 K15864             nirS  nitrite reductase (NO-forming) / hydroxylamine reductase     T1     T2 0.00079                 -2.08
# 4 K15864             nirS  nitrite reductase (NO-forming) / hydroxylamine reductase     T1     T3 0.00040                 -2.42
# 5 K20935              hdh                                   hydrazine dehydrogenase     T1     T3 0.01800                 -1.32

make_KO_treatment_scatterplot("K00370", title = get_KO_info("K00370"), point_size = 2)
make_KO_treatment_scatterplot("K15864", title = get_KO_info("K15864"), point_size = 2)
make_KO_treatment_scatterplot("K20935", title = get_KO_info("K20935"), point_size = 2)

###################################################################################################################
########## doing the same (N KOs) for the P treatment to find which above are NOT also sig in one of them #########
###################################################################################################################

# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
P_treatment_N_KOs_cov_subset_tab <- make_subset_cov_tab("", "P", N_pathway_KOs)

# getting long-form
P_treatment_N_KOs_cov_subset_tab_long <- make_longform_tab(P_treatment_N_KOs_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)
ggqqplot(P_treatment_N_KOs_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # most look fine here except for maybe K20932

## checking for heteroskedasticity
P_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ timepoint) %>% filter(p <= 0.05)
    # just 3 don't pass, K00371, K01501, K03385

## running regular anova
P_treatment_N_KOs_anova_tab <- P_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ timepoint) %>% data.frame(check.names = FALSE)
# adding KO annotations
P_treatment_N_KOs_anova_tab <- left_join(P_treatment_N_KOs_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_P_treatment_N_KOs_anova_tab <- P_treatment_N_KOs_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_P_treatment_N_KOs_anova_tab) # 9 of them

anova_sig_P_treatment_N_KOs <- anova_sig_P_treatment_N_KOs_anova_tab %>% pull(KO_ID) %>% as.vector() %>% unique()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_P_treatment_N_KOs)

# making subset table and running posthoc tukey
P_treatment_N_KOs_cov_subset_tab_long_sig <- P_treatment_N_KOs_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_P_treatment_N_KOs)
P_treatment_N_KOs_posthoc_tab <- P_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% tukey_hsd(cov ~ timepoint) %>% data.frame(check.names = FALSE)

# adding KO annotations
P_treatment_N_KOs_posthoc_tab <- left_join(P_treatment_N_KOs_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_P_treatment_N_KOs_posthoc_tab <- P_treatment_N_KOs_posthoc_tab %>% filter(p.adj.signif <= 0.1)

# adding log2fold changes
posthoc_sig_P_treatment_N_KOs_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_P_treatment_N_KOs_posthoc_tab, "P", P_treatment_N_KOs_cov_subset_tab_long)

# getting KO IDs
posthoc_sig_P_treatment_N_KOs <- posthoc_sig_P_treatment_N_KOs_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_P_treatment_N_KOs)

posthoc_sig_P_treatment_N_KOs_posthoc_tab %>% select(1,2,3,5,6,11,13)
#    KO_ID    KO_name               KO_def group1 group2   p.adj log2fc(group2/group1)
# 1 K01725       cynS        cyanate lyase     T1     T2 0.00169                 -2.06
# 2 K01915 glnA, GLUL glutamine synthetase     T1     T2 0.02060                 -0.79

## plotting those two
make_KO_treatment_scatterplot("K01725", title = get_KO_info("K01725"), point_size = 2)
make_KO_treatment_scatterplot("K01915", title = get_KO_info("K01915"), point_size = 2)



##### running all on just Nitrogen treatments


# making subset of all KOs, and Nitrogen treatments, to compare across timepoints
N_treatment_all_KOs_cov_subset_tab <- make_subset_cov_tab("", "N", all_our_KOs)

# getting long-form
N_treatment_all_KOs_cov_subset_tab_long <- make_longform_tab(N_treatment_all_KOs_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)

# ggqqplot(all_samples_N_KOs_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # hm, actually across all KOs it looks fine, except maybe K20932 points virtually all fall on the line, seems we can treat these as normally distributed

## checking for heteroskedasticity
N_treatment_all_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ timepoint) %>% filter(p <= 0.05) %>% length()
    # several don't pass

# looking at it a few ways
make_KO_treatment_scatterplot("K00371", wanted_treatments = "N")

make_KO_treatment_scatterplot("K00371")

make_single_KO_single_treatment_scatter_plot("K00371", "N")
make_single_KO_single_treatment_scatter_plot("K00371", "Control")

make_single_KO_single_timepoint_scatter_plot("K00371", "T1")
make_single_KO_single_timepoint_scatter_plot("K00371", "T3")



## running regular anova
N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ timepoint) %>% data.frame(check.names = FALSE)
# adding KO annotations
N_treatment_N_KOs_anova_tab <- left_join(N_treatment_N_KOs_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_N_treatment_N_KOs_anova_tab <- N_treatment_N_KOs_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_N_treatment_N_KOs_anova_tab)

anova_sig_N_treatment_N_KOs <- anova_sig_N_treatment_N_KOs_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_N_treatment_N_KOs, color_by_flumes = FALSE)

# making subset table and running posthoc tukey
N_treatment_N_KOs_cov_subset_tab_long_sig <- N_treatment_N_KOs_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_N_treatment_N_KOs)
N_treatment_N_KOs_posthoc_tab <- N_treatment_N_KOs_cov_subset_tab_long %>% group_by(KO_ID) %>% tukey_hsd(cov ~ timepoint) %>% data.frame(check.names = FALSE)

# adding KO annotations
N_treatment_N_KOs_posthoc_tab <- left_join(N_treatment_N_KOs_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_N_treatment_N_KOs_posthoc_tab <- N_treatment_N_KOs_posthoc_tab %>% filter(p.adj.signif <= 0.05)

# adding log2fold changes
posthoc_sig_N_treatment_N_KOs_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_N_treatment_N_KOs_posthoc_tab, "N", N_treatment_N_KOs_cov_subset_tab_long)

# getting KO IDs
posthoc_sig_N_treatment_N_KOs <- posthoc_sig_N_treatment_N_KOs_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_N_treatment_N_KOs)

posthoc_sig_N_treatment_N_KOs_posthoc_tab %>% select(1,2,3,5,6,11,13)
#     KO_ID               KO_name                                                    KO_def group1 group2   p.adj log2fc(group2/group1)



###### checking N-related functions across all timepoints and treatments













# plotting all with sig
make_KO_treatment_scatterplot(posthoc_sig_N_treatment_N_KOs, title = "Those with sig. diff. in N-treatment across time")
make_KO_treatment_scatterplot("K00371")


K00371_N_plot <- make_single_KO_single_treatment_scatter_plot("K00371", "N", color_by_flumes = FALSE)
K00371_P_plot <- make_single_KO_single_treatment_scatter_plot("K00371", "P", color_by_flumes = FALSE)
K00371_Control_plot <- make_single_KO_single_treatment_scatter_plot("K00371", "Control", color_by_flumes = FALSE)


ggarrange(K00371_N_plot, K00371_P_plot, K00371_Control_plot, ncol = 1)


make_single_KO_single_treatment_scatter_plot("K00371", "N")
make_single_KO_single_treatment_scatter_plot("K00371", "P")
make_single_KO_single_treatment_scatter_plot("K00371", "Control")

make_KO_treatment_scatterplot("K02588")
make_single_KO_single_treatment_scatter_plot("K02588", "N")
make_single_KO_single_treatment_scatter_plot("K02588", "P")
make_single_KO_single_treatment_scatter_plot("K02588", "Control")

make_single_KO_single_timepoint_scatter_plot("K02588", "T1", color_by_flumes = FALSE)
make_single_KO_single_timepoint_scatter_plot("K02588", "T2", color_by_flumes = FALSE)
make_single_KO_single_timepoint_scatter_plot("K02588", "T3", color_by_flumes = FALSE)

make_KO_treatment_scatterplot("K02586")
make_KO_treatment_scatterplot("K02591")


# K00371 is sig diff between T1 and T3 in nitrogen treatments













##############
# N vs Control at T2
# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
N_and_Control_N_KOs_T2_cov_subset_tab <- make_subset_cov_tab("T2", c("N", "Control"), N_pathway_KOs)

# getting long-form
N_and_Control_N_KOs_T2_cov_subset_tab_long <- make_longform_tab(N_and_Control_N_KOs_T2_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)
N_and_Control_N_KOs_T2_cov_subset_tab_long %>% group_by(KO_ID) %>% shapiro_test(cov)
    # some fail (p < 0.05), but the site notes the shapiro test can be unreliable misleading with > 50 samples and says qqplot is preferred
ggqqplot(N_and_Control_N_KOs_T2_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # hm, actually across all KOs it looks fine, except maybe K20932

## checking for heteroskedasticity
N_and_Control_N_KOs_T2_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ treatment) %>% filter(p <= 0.05)
    # all pass


## running regular anova
N_and_Control_N_KOs_T2_anova_tab <- N_and_Control_N_KOs_T2_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ treatment) %>% data.frame(check.names = FALSE)
# adding KO annotations
N_and_Control_N_KOs_T2_anova_tab <- left_join(N_and_Control_N_KOs_T2_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_N_and_Control_N_KOs_T2_anova_tab <- N_and_Control_N_KOs_T2_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_N_and_Control_N_KOs_T2_anova_tab)

anova_sig_N_and_Control_N_KOs_T2_KOs <- anova_sig_N_and_Control_N_KOs_T2_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_N_and_Control_N_KOs_T2_KOs)

# making subset table and running posthoc tukey
N_and_Control_N_KOs_T2_cov_subset_tab_long_sig <- N_and_Control_N_KOs_T2_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_N_and_Control_N_KOs_T2_KOs)
N_and_Control_N_KOs_T2_posthoc_tab <- N_and_Control_N_KOs_T2_cov_subset_tab_long_sig %>% group_by(KO_ID) %>% tukey_hsd(cov ~ treatment) %>% data.frame(check.names = FALSE)

# adding KO annotations
N_and_Control_N_KOs_T2_posthoc_tab <- left_join(N_and_Control_N_KOs_T2_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_N_and_Control_N_KOs_T2_posthoc_tab <- N_and_Control_N_KOs_T2_posthoc_tab %>% filter(p.adj <= 0.05)

# adding log2fold changes
posthoc_sig_N_and_Control_N_KOs_T2_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_N_and_Control_N_KOs_T2_posthoc_tab, "N", N_and_Control_N_KOs_T2_cov_subset_tab_long_sig)

# getting KO IDs
posthoc_sig_N_and_Control_N_KOs_T2_KOs <- posthoc_sig_N_and_Control_N_KOs_T2_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_N_and_Control_N_KOs_T2_KOs)

posthoc_sig_N_and_Control_N_KOs_T2_posthoc_tab %>% select(1,2,3,5,6,11,13)

# looking at them

make_KO_treatment_scatterplot("K00368", title = get_KO_info("K00368"), point_size = 2)
make_KO_treatment_scatterplot("K15864", title = get_KO_info("K15864"), point_size = 2)

#######

# N vs Control at T3
# making subset of just N-related funtions, and Nitrogen treatments, to compare across timepoints
N_and_Control_N_KOs_T3_cov_subset_tab <- make_subset_cov_tab("T3", c("N", "Control"), N_pathway_KOs)

# getting long-form
N_and_Control_N_KOs_T3_cov_subset_tab_long <- make_longform_tab(N_and_Control_N_KOs_T3_cov_subset_tab)

## checking for normality by KO_ID (across all timepoints)
ggqqplot(N_and_Control_N_KOs_T3_cov_subset_tab_long, "cov", facet.by = "KO_ID")
    # hm, actually across all KOs it looks fine, except maybe K20932
get_KO_info("K20932")

## checking for heteroskedasticity
N_and_Control_N_KOs_T3_cov_subset_tab_long %>% group_by(KO_ID) %>% levene_test(cov ~ treatment) %>% filter(p <= 0.05)
    # all pass


## running regular anova
N_and_Control_N_KOs_T3_anova_tab <- N_and_Control_N_KOs_T3_cov_subset_tab_long %>% group_by(KO_ID) %>% anova_test(cov ~ treatment) %>% data.frame(check.names = FALSE)
# adding KO annotations
N_and_Control_N_KOs_T3_anova_tab <- left_join(N_and_Control_N_KOs_T3_anova_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering to those below 0.1
anova_sig_N_and_Control_N_KOs_T3_anova_tab <- N_and_Control_N_KOs_T3_anova_tab %>% filter(p <= 0.1)
dim(anova_sig_N_and_Control_N_KOs_T3_anova_tab)

anova_sig_N_and_Control_N_KOs_T3_KOs <- anova_sig_N_and_Control_N_KOs_T3_anova_tab %>% pull(KO_ID) %>% as.vector()

# plotting all sig # useless with so many as the y-axes can't be free fo rall plots with facet_grid
make_KO_treatment_scatterplot(anova_sig_N_and_Control_N_KOs_T3_KOs)

# making subset table and running posthoc tukey
N_and_Control_N_KOs_T3_cov_subset_tab_long_sig <- N_and_Control_N_KOs_T3_cov_subset_tab_long %>% filter(KO_ID %in% anova_sig_N_and_Control_N_KOs_T3_KOs)
N_and_Control_N_KOs_T3_posthoc_tab <- N_and_Control_N_KOs_T3_cov_subset_tab_long_sig %>% group_by(KO_ID) %>% tukey_hsd(cov ~ treatment) %>% data.frame(check.names = FALSE)

# adding KO annotations
N_and_Control_N_KOs_T3_posthoc_tab <- left_join(N_and_Control_N_KOs_T3_posthoc_tab, Master_KEGG_tab %>% select(1,2,3)) %>% relocate(KO_name, KO_def, .after = KO_ID)

# filtering down to those with adj. p <= 0.1
posthoc_sig_N_and_Control_N_KOs_T3_posthoc_tab <- N_and_Control_N_KOs_T3_posthoc_tab %>% filter(p.adj <= 0.1)

# adding log2fold changes
posthoc_sig_N_and_Control_N_KOs_T3_posthoc_tab <- add_log2fc_to_posthoc_tab(posthoc_sig_N_and_Control_N_KOs_T3_posthoc_tab, "N", N_and_Control_N_KOs_T3_cov_subset_tab_long_sig)

# getting KO IDs
posthoc_sig_N_and_Control_N_KOs_T3_KOs <- posthoc_sig_N_and_Control_N_KOs_T3_posthoc_tab %>% pull(KO_ID) %>% as.vector() %>% unique()
length(posthoc_sig_N_and_Control_N_KOs_T3_KOs)

posthoc_sig_N_and_Control_N_KOs_T3_posthoc_tab %>% select(1,2,3,5,6,11,13)

# looking at them

make_KO_treatment_scatterplot("K00368", title = get_KO_info("K00368"), point_size = 2)
make_KO_treatment_scatterplot("K15864", title = get_KO_info("K15864"), point_size = 2)


