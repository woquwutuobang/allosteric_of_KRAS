library(krasddpcams)

# Function to merge predicted and observed fitness data for 3 experimental blocks
K13_19_get_ob_pre_fitness_binding_correlation_3blocks <- function(
    prediction = "path/to/predicted_phenotypes_all.txt",
    block1_dimsum_df = "path/to/block1_fitness.RData",
    block2_dimsum_df = "path/to/block2_fitness.RData", 
    block3_dimsum_df = "path/to/block3_fitness.RData",
    assay_sele = "RAF1", 
    wt_aa_input = wt_aa) {
  
  # Read prediction data
  pre <- fread(prediction)
  pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
  
  # Load all 3 block datasets
  load(block1_dimsum_df)
  block1 <- as.data.table(all_variants)
  load(block2_dimsum_df)
  block2 <- as.data.table(all_variants)
  load(block3_dimsum_df)
  block3 <- as.data.table(all_variants)
  
  # Merge all block data
  data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3,
                           idcol = "block", fill = TRUE)
  
  # Calculate weighted fitness intermediates
  data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness/(data_before_nor$sigma)^2
  data_before_nor$one_over_fitness_sigmasquared <- 1/(data_before_nor$sigma)^2
  
  # Calculate stop codon fitness for each block (weighted average)
  calculate_stop_fitness <- function(block_name) {
    dead_fitness <- data_before_nor[STOP == TRUE & block == block_name, ]
    stop_fitness <- sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(stop_fitness)
  }
  
  stop1_fitness <- calculate_stop_fitness("block1")
  stop2_fitness <- calculate_stop_fitness("block2")
  stop3_fitness <- calculate_stop_fitness("block3")
  
  # Calculate wild-type fitness for each block (weighted average)
  calculate_wt_fitness <- function(block_name) {
    wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
    wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(wt_fitness_block$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(wt_fitness)
  }
  
  wt1_fitness <- calculate_wt_fitness("block1")
  wt2_fitness <- calculate_wt_fitness("block2")
  wt3_fitness <- calculate_wt_fitness("block3")
  
  # Create scaling data for normalization
  scaling_data_fitness <- data.frame(
    block1 = c(stop1_fitness, wt1_fitness),
    block2 = c(stop2_fitness, wt2_fitness),
    block3 = c(stop3_fitness, wt3_fitness)
  )
  
  # Calculate scaling parameters relative to block1
  calculate_scaling_params <- function(target_block) {
    lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
    return(list(
      slope = lm_model$coefficients[[2]],
      intercept = lm_model$coefficients[[1]]
    ))
  }
  
  params_block2 <- calculate_scaling_params("block2")
  params_block3 <- calculate_scaling_params("block3")
  
  d2 <- params_block2$slope
  e2 <- params_block2$intercept
  d3 <- params_block3$slope
  e3 <- params_block3$intercept
  
  # Process prediction data
  pre_nor <- pre_pos
  
  # Extract predicted fitness
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # Assay to phenotype mapping
  assay_sele_df <- data.table(assay = c("K13", "K19", "K27", "K55", "PI3", "RAF1", "RAL", "SOS"), 
                              phenotype_base = c(6, 11, 16, 21, 26, 31, 36, 41))
  
  # Validate assay selection
  if (!assay_sele %in% assay_sele_df$assay) {
    stop("Assay '", assay_sele, "' not found. Available assays: ", 
         paste(assay_sele_df$assay, collapse = ", "))
  }
  
  # Get base phenotype number
  base_pheno <- assay_sele_df[assay == assay_sele, phenotype_base]
  
  # Apply normalization transformations for 3 blocks
  # Block 1
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_mean_fitness, mean)]
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness_sigma, std)]
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness, fitness)]
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness_sigma, sigma)]
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness, predicted_fitness)]
  
  # Block 2
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness_sigma, std * d2)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness, fitness * d2 + e2)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness_sigma, sigma * d2)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
  
  # Block 3
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness_sigma, std * d3)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness, fitness * d3 + e3)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness_sigma, sigma * d3)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
  
  return(pre_nor)
}

######################################

####################################
# 绘图函数保持不变
krasddpcams__plot2d_ddGb_ob_pre_fitness_perblock <- function(pre_nor = pre_nor, phenotypen = phenotypen, rotate_x_axis = TRUE) 
{
  pre_nor
  lm_mochi <- lm(pre_nor_fitness ~ ob_nor_fitness, pre_nor[phenotype == phenotypen, ])
  
  
  p <- ggplot2::ggplot() + 
    ggplot2::stat_binhex(data = pre_nor[phenotype == phenotypen, ], 
                         ggplot2::aes(x = ob_nor_fitness, y = pre_nor_fitness), 
                         bins = 50, size = 0, color = "black") + 
    ggplot2::scale_fill_gradient(low = "white", 
                                 high = "black", 
                                 trans = "log10", 
                                 guide = ggplot2::guide_colorbar(barwidth = 0.5, 
                                                                 barheight = 1.5)) + 
    ggplot2::geom_hline(yintercept = 0) + 
    ggplot2::geom_vline(xintercept = 0) + 
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    ggplot2::annotate("text", 
                      x = -0.8, y = 0.3, 
                      label = paste0("R² = ", round(summary(lm_mochi)$r.squared, 2)), 
                      size = 8 * 0.35) + 
    ggplot2::theme_classic() + 
    ggplot2::xlab("Observed fitness") + 
    ggplot2::ylab("Predicted fitness") + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 8), 
      axis.text = ggplot2::element_text(size = 8), 
      legend.text = ggplot2::element_text(size = 8), 
      legend.key.size = ggplot2::unit(1, "cm"), 
      plot.title = ggplot2::element_text(size = 8)
    ) + 
    ggplot2::coord_fixed()
  
  # Determines whether to rotate the X axis according to the parameters
  if (rotate_x_axis) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  
  return(p)
}

# Wild-type KRAS sequence (保持不变)
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 使用3块版本的数据处理函数
# 注意：这里需要根据您的实际数据调整路径和block数量

# 示例：RAF1 数据处理（调整为3个block）
# RAF1_pre_ob_fitness <- K13_19_get_ob_pre_fitness_binding_correlation_3blocks(
#   prediction = "path/to/predicted_phenotypes_all.txt",
#   block1_dimsum_df = "path/to/RAF_block1.RData",
#   block2_dimsum_df = "path/to/RAF_block2.RData", 
#   block3_dimsum_df = "path/to/RAF_block3.RData",
#   assay_sele = "RAF1",
#   wt_aa_input = wt_aa)

# 示例：绘图调用（调整为3个图）
# krasddpcams__plot2d_ddGb_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 31, rotate_x_axis = TRUE)
# krasddpcams__plot2d_ddGb_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 32, rotate_x_axis = TRUE)
# krasddpcams__plot2d_ddGb_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 33, rotate_x_axis = TRUE)