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

# FUNCTION TO MAKE LONG-FORM TABLE
make_longform_tab <- function(input_tab, sam_info_tab = sample_info_tab) {

    curr_long_tab <- input_tab %>% select(-2) %>% pivot_longer(!KO_ID, names_to = "sample", values_to = "cov")

    curr_merged_tab <- left_join(curr_long_tab, sam_info_tab)

    return(curr_merged_tab)
}

# make_longform_tab(N_treatment_N_KOs_cov_subset_tab)

# FUNCTION TO CALCULATE LOG2FOLD CHANGE BETWEEN 2 TIMEPOINTS OF A SPECIFIC TREATMENT
calc_log2fc_between_timepoints_for_a_treatment <- function(wanted_KO, wanted_treatment, long_form_coverage_tab, contrast_timepoint_numerator = "T3", contrast_timepoint_denominator = "T1") {

    # getting numerator target mean coverage
    curr_numerator_mean_cov <- long_form_coverage_tab %>% filter(KO_ID == wanted_KO, treatment == wanted_treatment, timepoint == contrast_timepoint_numerator) %>% pull(cov) %>% mean()
    curr_denominator_mean_cov <- long_form_coverage_tab %>% filter(KO_ID == wanted_KO, treatment == wanted_treatment, timepoint == contrast_timepoint_denominator) %>% pull(cov) %>% mean()

    curr_foldchange <- curr_numerator_mean_cov / curr_denominator_mean_cov
    curr_log2fc <- log2(curr_foldchange) %>% round(2)

    return(curr_log2fc)

}

# FUNCTION TO ADD LOG2FOLD CHANGES TO POSTHOC TABLE
add_log2fc_to_posthoc_tab <- function(posthoc_tab, treatment, long_form_coverage_tab) {

    # initalizing new vector
    curr_log2fc_vec <- c()

    # iterating through table and adding log2fold changes to it
    for ( row in 1:dim(posthoc_tab)[1] ) {

        curr_row <- posthoc_tab[row, ]
        curr_KO <- curr_row %>% pull(KO_ID) %>% as.vector()
        curr_group2 <- curr_row %>% pull(group2) %>% as.vector()
        curr_group1 <- curr_row %>% pull(group1) %>% as.vector()

        curr_log2fc <- calc_log2fc_between_timepoints_for_a_treatment(curr_KO, treatment, long_form_coverage_tab, curr_group2, curr_group1)
        curr_log2fc_vec <- c(curr_log2fc_vec, curr_log2fc)

    }

    # adding new column of log2fold changes
    posthoc_tab$"log2fc(group2/group1)" <- curr_log2fc_vec

    return(posthoc_tab)

}

# FUNCTION TO MAKE SINGLE-KO, SINGLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_single_KO_single_treatment_scatter_plot <- function(wanted_KO, wanted_treatment, wanted_timepoints = "", color_by_flumes = FALSE, input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab, kegg_info_tab = Master_KEGG_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoints, wanted_treatment, wanted_KO, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_KO_name <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_name) %>% as.vector()
    curr_KO_def <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_def) %>% as.vector()

    curr_title <- paste0(wanted_KO, " (Treatment: ", wanted_treatment, ")")
    curr_subtitle <- paste0(curr_KO_name, " | ", curr_KO_def)

    curr_depth_colors_vec <- unique(curr_sub_tab_long$depth_color[order(curr_sub_tab_long$depth)]) %>% as.vector()
    curr_flume_colors_vec <- unique(curr_sub_tab_long$flume_color[order(curr_sub_tab_long$flume)]) %>% as.vector()

    if ( color_by_flumes == FALSE ) {

        curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = depth), size = 3) +
            scale_color_manual(values = curr_depth_colors_vec) + xlab("Timepoint") +
             ylab("Norm. Coverage (Rib. Prots.)") + ggtitle(curr_title, subtitle = curr_subtitle) + theme_bw() + theme(legend.position = "bottom")

    } else {

        curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = flume), size = 3) +
            scale_color_manual(values = curr_flume_colors_vec) + xlab("Timepoint") +
             ylab("Norm. Coverage (Rib. Prots.)") + ggtitle(curr_title, subtitle = curr_subtitle) + theme_bw() + theme(legend.position = "bottom")

    }

    return(curr_plot)

}


# FUNCTION TO MAKE SINGLE-KO, SINGLE-TIMEPOINT SCATTERPLOTS BASED ON TREATMENT
make_single_KO_single_timepoint_scatter_plot <- function(wanted_KO, wanted_timepoint, wanted_treatments = "", color_by_flumes = FALSE, input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab, kegg_info_tab = Master_KEGG_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoint, wanted_treatments, wanted_KO, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_KO_name <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_name) %>% as.vector()
    curr_KO_def <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_def) %>% as.vector()

    curr_title <- paste0(wanted_KO, " (Timepoint: ", wanted_timepoint, ")")
    curr_subtitle <- paste0(curr_KO_name, " | ", curr_KO_def)

    curr_depth_colors_vec <- unique(curr_sub_tab_long$depth_color[order(curr_sub_tab_long$depth)]) %>% as.vector()
    curr_flume_colors_vec <- unique(curr_sub_tab_long$flume_color[order(curr_sub_tab_long$flume)]) %>% as.vector()

    if ( color_by_flumes == FALSE ) {

        curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = treatment, y = cov, shape = depth, color = depth), size = 3) +
            scale_color_manual(values = curr_depth_colors_vec) + xlab("Treatment") +
             ylab("Norm. Coverage (Rib. Prots.)") + ggtitle(curr_title, subtitle = curr_subtitle) + theme_bw() + theme(legend.position = "bottom")

    } else {

        curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = treatment, y = cov, shape = depth, color = flume), size = 3) +
            scale_color_manual(values = curr_flume_colors_vec) + xlab("Treatment") +
             ylab("Norm. Coverage (Rib. Prots.)") + ggtitle(curr_title, subtitle = curr_subtitle) + theme_bw() + theme(legend.position = "bottom")

    }

    return(curr_plot)

}

# FUNCTION FOR GETTING KO NAME AND DEF
get_KO_info <- function(wanted_KO, kegg_info_tab = Master_KEGG_tab) {

    curr_KO_name <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_name) %>% as.vector()
    curr_KO_def <- kegg_info_tab %>% filter(KO_ID == wanted_KO) %>% pull(KO_def) %>% as.vector()

    curr_title <- paste0(wanted_KO, " (", curr_KO_name, " | ", curr_KO_def, ")")

    return(curr_title)

}

# FUNCTION TO MAKE MULTIPLE-KO, SINGLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_multi_KO_single_treatment_scatter_plot <- function(wanted_KOs, wanted_treatment, wanted_timepoints = "", input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab) {

    curr_sub_tab <- make_subset_cov_tab(wanted_timepoints, wanted_treatment, wanted_KOs, input_cov_tab = input_cov_tab, sam_info_tab = sam_info_tab)

    curr_sub_tab_long <- make_longform_tab(curr_sub_tab)

    curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov)) + xlab("Timepoint") +
        ylab("Norm. Coverage (Rib. Prots.)") + ggtitle(paste0("Treatment: ", wanted_treatment)) + facet_wrap(~KO_ID, scales = "free_y") +
        theme_bw()

    return(curr_plot)
}

# make_multiple_KO_single_treatment_scatter_plot(c("K00371", "K02568"), "N")

# FUNCTION TO MAKE MULTIPLE-KO, MULTIPLE-TREATMENT SCATTERPLOTS BASED ON TIMEPOINT
make_KO_treatment_scatterplot <- function(wanted_KOs, wanted_treatments = "", wanted_timepoints = "", title = "", color_by_flumes = FALSE, input_cov_tab = all_cov_norm, sam_info_tab = sample_info_tab, point_size = 3, scale_setting = "free_y") {

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

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = depth), size = point_size) +
                scale_color_manual(values = curr_depth_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (Rib. Prots.)") + facet_grid(treatment ~ KO_ID, scales = scale_setting) +
                theme_bw() + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

        } else {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = flume), size = point_size) +
                scale_color_manual(values = curr_flume_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (Rib. Prots.)") + facet_grid(treatment ~ KO_ID, scales = scale_setting) +
                theme_bw() + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

        }

    } else {

        if ( color_by_flumes == FALSE ) {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = depth), size = point_size) +
                scale_color_manual(values = curr_depth_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (Rib. Prots.)") + facet_grid(treatment ~ KO_ID, scales = scale_setting) + ggtitle(title) +
                theme_bw() + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

        } else {

            curr_plot <- ggplot(curr_sub_tab_long) + geom_point(aes(x = timepoint, y = cov, shape = depth, color = flume), size = point_size) +
                scale_color_manual(values = curr_flume_colors_vec) + xlab("Timepoint") +
                ylab("Norm. Coverage (Rib. Prots.)") + facet_grid(treatment ~ KO_ID, scales = scale_setting) + ggtitle(title) +
                theme_bw() + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

        }

    }

    return(curr_plot)
}

# make_KO_treatment_scatterplot(c("K00371", "K02568"))
# make_KO_treatment_scatterplot("K00371")
# make_KO_treatment_scatterplot("K00371", color_by_flumes = FALSE)
# make_KO_treatment_scatterplot("K02933")

# FUNCTION TO MAKE DENDROGRAMS
make_dendros <- function(input_cov_tab, wanted_timepoints = "", wanted_treatments = "", wanted_depths = "", wanted_KOs = "", sam_info_tab = sample_info_tab, include_environmental = FALSE, include_unannotated = FALSE, additional_subtitle = "") {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }

    if ( length(wanted_depths) == 1 && wanted_depths == "" ) {
        wanted_depths <- c("0-1mm", "1-2mm", "2-3mm", "3-4mm")
    }

    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments), depth %in% all_of(wanted_depths)) %>% pull(sample) %>% as.vector

    # dropping unannotated row unless specified to keep
    if ( include_unannotated == FALSE ) {

        curr_tab <- input_cov_tab %>% filter(KO_ID != "Not annotated")

    } else {

        curr_tab <- input_cov_tab
    }

    # if wanted KOs were provided, subsetting here
    if ( length(wanted_KOs == 1) && ! wanted_KOs == "" ) {

        # subsetting table
        curr_tab <- curr_tab %>% filter(KO_ID %in% wanted_KOs)

    }

    # adding environmental samples if specified
    if ( include_environmental ) {

        curr_tab <- curr_tab %>% select(all_of(c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))))
        curr_wanted_samples <- c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))

    } else {

        curr_tab <- curr_tab %>% select(all_of(curr_wanted_samples))

    }

    # making dendrogram
    curr_dist <- vegdist(t(curr_tab))
    curr_hclust <- hclust(curr_dist, method = "ward.D2")
    curr_dendro <- dendsort(as.dendrogram(curr_hclust))

    # making colors and depth shapes vectors
    curr_dendro_w_timepoint_cols <- curr_dendro
    curr_dendro_w_treatment_cols <- curr_dendro
    curr_dendro_w_depth_cols <- curr_dendro
    curr_dendro_w_flume_cols <- curr_dendro

    curr_timepoint_cols_vec <- vector()
    curr_treatment_cols_vec <- vector()
    curr_depth_cols_vec <- vector()
    curr_flume_cols_vec <- vector()

    curr_depth_shapes_vec <- vector()

    for ( curr_label in labels(curr_dendro) ) {

        curr_timepoint_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(timepoint_color) %>% as.vector
        curr_timepoint_cols_vec <- c(curr_timepoint_cols_vec, curr_timepoint_col)

        curr_treatment_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(treatment_color) %>% as.vector
        curr_treatment_cols_vec <- c(curr_treatment_cols_vec, curr_treatment_col)

        curr_depth_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(depth_color) %>% as.vector
        curr_depth_cols_vec <- c(curr_depth_cols_vec, curr_depth_col)

        curr_flume_col <- sam_info_tab %>% filter(sample == curr_label) %>% pull(flume_color) %>% as.vector
        curr_flume_cols_vec <- c(curr_flume_cols_vec, curr_flume_col)

        curr_depth_shape <- sam_info_tab %>% filter(sample == curr_label) %>% pull(depth_pch) %>% as.vector
        curr_depth_shapes_vec <- c(curr_depth_shapes_vec, curr_depth_shape)

    }

    labels_colors(curr_dendro_w_timepoint_cols) <- curr_timepoint_cols_vec
    labels_colors(curr_dendro_w_treatment_cols) <- curr_treatment_cols_vec
    labels_colors(curr_dendro_w_depth_cols) <- curr_depth_cols_vec
    labels_colors(curr_dendro_w_flume_cols) <- curr_flume_cols_vec

    # making base title
    if ( length(wanted_depths) == 4 ) {

        curr_base_title <- paste0("Timepoint(s): ", paste0(wanted_timepoints, collapse = ","), "  |  Treatment(s): ", paste0(wanted_treatments, collapse = ", "), "\n")

    } else {
        # adding which depths if not all
        curr_base_title <- paste0("Timepoint(s): ", paste0(wanted_timepoints, collapse = ","), "  |  Treatment(s): ", paste0(wanted_treatments, collapse = ", "), "  |  Depth(s): ", paste0(wanted_depths, collapse = ", "), "\n")
    }

    # adding additional subtitle if one was specified
    if ( additional_subtitle != "" ) {
        curr_base_title <- paste0(curr_base_title, "(", additional_subtitle, ")")
    }

    # bottom title for depth-shape relationships
    curr_bottom_text <- "Shapes indicate depths: Black Square = 0-1mm; Blue Circle = 1-2mm; Red Triangle = 2-3mm; Brown Diamond = 3-4mm\n"

    # making individual plot subtitles
    curr_timepoint_subtitle <- paste0("Labels colored by timepoint: Black = T1; Blue = T2; Red = T3\n", curr_bottom_text)
    curr_treatment_subtitle <- paste("Labels colored by treatment: Black = Control; Red = P; Blue = N\n", curr_bottom_text)
    curr_depth_subtitle <- paste("Labels colored by depth: Black = 0-1mm; Blue = 1-2mm; Red = 2-3mm; Brown = 3-4mm\n", curr_bottom_text)
    curr_flume_subtitle <- paste("Labels colored by flume: Black = F1; Blue = F2; Red = F3; Brown = F4; Purple = F5; Teal = F6\n", curr_bottom_text)

    # making plots
    curr_dendro_w_timepoint_cols %>% set("leaves_pch", curr_depth_shapes_vec) %>% set("leaves_cex", 1.5) %>% set("leaves_col", curr_depth_cols_vec) %>%
        plot(main = curr_base_title, sub = curr_timepoint_subtitle, ylab="Bray-Curtis dissimilarity")
    curr_timepoint_plot <- recordPlot()

    curr_dendro_w_treatment_cols %>% set("leaves_pch", curr_depth_shapes_vec) %>% set("leaves_cex", 1.5) %>% set("leaves_col", curr_depth_cols_vec) %>%
        plot(main = curr_base_title, sub = curr_treatment_subtitle, ylab="Bray-Curtis dissimilarity")
    curr_treatment_plot <- recordPlot()

    curr_dendro_w_depth_cols %>% set("leaves_pch", curr_depth_shapes_vec) %>% set("leaves_cex", 1.5) %>% set("leaves_col", curr_depth_cols_vec) %>%
        plot(main = curr_base_title, sub = curr_depth_subtitle, ylab="Bray-Curtis dissimilarity")
    curr_depth_plot <- recordPlot()

    curr_dendro_w_flume_cols %>% set("leaves_pch", curr_depth_shapes_vec) %>% set("leaves_cex", 1.5) %>% set("leaves_col", curr_depth_cols_vec) %>%
        plot(main = curr_base_title, sub = curr_flume_subtitle, ylab="Bray-Curtis dissimilarity")
    curr_flume_plot <- recordPlot()

    return(list(timepoint_colored_plot = curr_timepoint_plot, treatment_colored_plot = curr_treatment_plot, depth_colored_plot = curr_depth_plot, flume_colored_plot = curr_flume_plot))

}

# FUNCTION TO MAKE ORDINATIONS
make_ordinations <- function(input_cov_tab, wanted_timepoints = "", wanted_treatments = "", wanted_depths = "", wanted_KOs = "", sam_info_tab = sample_info_tab, include_environmental = FALSE, include_unannotated = FALSE, additional_subtitle = "", shape_size = 2.5) {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }

    if ( length(wanted_depths) == 1 && wanted_depths == "" ) {
        wanted_depths <- c("0-1mm", "1-2mm", "2-3mm", "3-4mm")
    }

    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments), depth %in% all_of(wanted_depths)) %>% pull(sample) %>% as.vector

    # dropping unannotated row unless specified to keep
    if ( include_unannotated == FALSE ) {

        curr_tab <- input_cov_tab %>% filter(KO_ID != "Not annotated")

    } else {

        curr_tab <- input_cov_tab
    }

    # if wanted KOs were provided, subsetting here
    if ( length(wanted_KOs == 1) && ! wanted_KOs == "" ) {

        # subsetting table
        curr_tab <- curr_tab %>% filter(KO_ID %in% wanted_KOs)

    }

    # adding environmental samples if specified
    if ( include_environmental ) {

        curr_tab <- curr_tab %>% select(all_of(c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))))
        curr_wanted_samples <- c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))

    } else {

        curr_tab <- curr_tab %>% select(all_of(curr_wanted_samples))

    }

    # subsetting sample info tab to what's being focused on here
    sam_info_tab <- sam_info_tab %>% filter(sample %in% curr_wanted_samples)
    # changing factors to characters so they don't mess with the plotting below
    sam_info_tab <- sam_info_tab %>% mutate(across(where(is.factor), as.character))

    # getting dist
    curr_dist <- vegdist(t(curr_tab))

    # running pcoa
    curr_pcoa <- pcoa(curr_dist)

    curr_tab <- data.frame("Axis.1" = curr_pcoa$vectors[,1], "Axis.2" = curr_pcoa$vectors[,2], sam_info_tab)

    # getting values of first and second axis eigenvalues
    curr_eig_1 <- curr_pcoa$values$Eigenvalues[1]
    curr_eig_2 <- curr_pcoa$values$Eigenvalues[2]

    # this allows us to scale the axes according to their magnitude of separating apart the samples
    curr_coord_adj <- sqrt(curr_eig_2 / curr_eig_2)

    # getting values of axis contribution for labels on plot
    curr_eigen_sum <- sum(curr_pcoa$values$Eigenvalues)
    curr_eig_1_contr <- round(curr_eig_1 / curr_eigen_sum * 100, 1)
    curr_eig_2_contr <- round(curr_eig_2 / curr_eigen_sum * 100, 1)

    # making axis labels
    curr_x_axis_label <- paste0("Axis.1 (", curr_eig_1_contr, "%)")
    curr_y_axis_label <- paste0("Axis.2 (", curr_eig_2_contr, "%)")


    # making title
    if ( length(wanted_depths) == 4 ) {

        curr_title <- paste0("Timepoint(s): ", paste0(wanted_timepoints, collapse = ","), "  |  Treatment(s): ", paste0(wanted_treatments, collapse = ", "))

    } else {
        # adding which depths if not all
        curr_title <- paste0("Timepoint(s): ", paste0(wanted_timepoints, collapse = ","), "  |  Treatment(s): ", paste0(wanted_treatments, collapse = ", "), "  |  Depth(s): ", paste0(wanted_depths, collapse = ", "))

    }

    curr_plot_with_timepoint_cols <- ggplot(curr_tab, aes(x = Axis.1, y = Axis.2, color = timepoint, shape = depth)) +
        geom_point(size = shape_size) + scale_color_manual(values = unique(curr_tab$timepoint_color)) +
        scale_shape_manual(values = curr_tab$depth_pch) + xlab(curr_x_axis_label) +
        ylab(curr_y_axis_label) + coord_fixed(curr_coord_adj) + theme_bw() +
        ggtitle(curr_title, subtitle = additional_subtitle) + theme(legend.position = "bottom")


    curr_plot_with_treatment_cols <- ggplot(curr_tab, aes(x = Axis.1, y = Axis.2, color = treatment, shape = depth)) +
        geom_point(size = shape_size) + scale_color_manual(values = rev(unique(curr_tab$treatment_color))) +
        scale_shape_manual(values = curr_tab$depth_pch) + xlab(curr_x_axis_label) +
        ylab(curr_y_axis_label) + coord_fixed(curr_coord_adj) + theme_bw() +
        ggtitle(curr_title, subtitle = additional_subtitle) + theme(legend.position = "bottom")

    curr_plot_with_depth_cols <- ggplot(curr_tab, aes(x = Axis.1, y = Axis.2, color = depth, shape = depth)) +
        geom_point(size = shape_size) + scale_color_manual(values = curr_tab$depth_color) +
        scale_shape_manual(values = curr_tab$depth_pch) + xlab(curr_x_axis_label) +
        ylab(curr_y_axis_label) + coord_fixed(curr_coord_adj) + theme_bw() +
        ggtitle(curr_title, subtitle = additional_subtitle) + theme(legend.position = "bottom")


    curr_plot_with_flume_cols <- ggplot(curr_tab, aes(x = Axis.1, y = Axis.2, color = flume, shape = depth)) +
        geom_point(size = shape_size) + scale_color_manual(values = unique(curr_tab$flume_color)) +
        scale_shape_manual(values = curr_tab$depth_pch) + xlab(curr_x_axis_label) +
        ylab(curr_y_axis_label) + coord_fixed(curr_coord_adj) + theme_bw() +
        ggtitle(curr_title, subtitle = additional_subtitle) + theme(legend.position = "bottom")

    return(list("timepoint_colored_plot" = curr_plot_with_timepoint_cols, "treatment_colored_plot" = curr_plot_with_treatment_cols, "depth_colored_plot" = curr_plot_with_depth_cols, "flume_colored_plot" = curr_plot_with_flume_cols))

}

# FUNCTION TO GET DISTANCE MATRIX
make_dist_mat <- function(input_cov_tab, wanted_timepoints = "", wanted_treatments = "", wanted_KOs = "", sam_info_tab = sample_info_tab, include_environmental = FALSE, include_unannotated = FALSE) {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }
    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments)) %>% pull(sample) %>% as.vector

    # dropping unannotated row unless specified to keep
    if ( include_unannotated == FALSE ) {

        curr_tab <- input_cov_tab %>% filter(KO_ID != "Not annotated")

    } else {

        curr_tab <- input_cov_tab
    }

    # if wanted KOs were provided, subsetting here
    if ( length(wanted_KOs == 1) && ! wanted_KOs == "" ) {

        # subsetting table
        curr_tab <- curr_tab %>% filter(KO_ID %in% wanted_KOs)

    }

    # adding environmental samples if specified
    if ( include_environmental ) {

        curr_tab <- curr_tab %>% select(all_of(c(curr_wanted_samples, c("F4N1", "F4N5", "F5N4", "F6", "F8", "FA1", "FA2", "FA3", "FPS"))))

    } else {

        curr_tab <- curr_tab %>% select(all_of(curr_wanted_samples))

    }

    # making dendrogram
    curr_dist <- vegdist(t(curr_tab))

    return(curr_dist)

}

# FUNCTION TO MAKE SUBSET SAMPLE INFO TAB
make_subset_sample_tab <- function(wanted_timepoints = "", wanted_treatments = "", sam_info_tab = sample_info_tab, include_environmental = FALSE) {

    # getting wanted sample IDs
    # if entered empty, using all available (doesn't capture Environmental ones)
    if ( length(wanted_timepoints) == 1 && wanted_timepoints == "" ) {
        wanted_timepoints <- c("T1", "T2", "T3")
    }

    if ( length(wanted_treatments) == 1 && wanted_treatments == "" ) {
        wanted_treatments <- c("P", "N", "Control")
    }
    curr_wanted_samples <- sam_info_tab %>% filter(timepoint %in% all_of(wanted_timepoints), treatment %in% all_of(wanted_treatments)) %>% pull(sample) %>% as.vector

    # subsetting sample info tab to what's being focused on here
    sam_info_tab <- sam_info_tab %>% filter(sample %in% curr_wanted_samples)

    # changing factors to characters as right now i think characters are less likely to be problematic in whatever i'm doing with these subset tables
    sam_info_tab <- sam_info_tab %>% mutate(across(where(is.factor), as.character))

    return(sam_info_tab)

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