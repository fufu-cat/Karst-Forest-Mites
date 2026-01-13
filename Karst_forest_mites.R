# Title: (Research article)
# Description: The R code was used to generate all statistical results and figures in the main text.
# Date: November 2025
# R4.3.2




library(vegan) 
library(performance) 
library(agricolae)
library(rstatix)
library(tidyverse)
library(reshape2) 
library(Rmisc) 
library(ggplot2)
library(RColorBrewer)
library(multifunc)
library(ggpmisc)
library(sjPlot)
library(lme4)
library(ggeffects)
library(cowplot)
library(pheatmap)
library(adespatial)
library(NST)
library(ggtern)
library(spaa)
library(Hmisc)   
library(igraph)
library(MetaNet)
library(MuMIn)
library(piecewiseSEM)
library(car)



#### 1 Multitrophic diversity within mite communities - alpha diversity (Figure 1) ####
alpha_diversity1 <- function(data, id_cols = 1:2) {
  require(vegan)
  result <- data.frame(
    data[, id_cols, drop = FALSE],
    Richness = apply(data[, -id_cols], 1, function(x) sum(x > 0)),
    Shannon = diversity(data[, -id_cols], index = "shannon")
  )
  return(result)
}

m_total <- read.csv("mite.csv")
alpha_total <- alpha_diversity1(m_total)
m_fungi <- read.csv("m_fungi.csv")
alpha_fungi <- alpha_diversity1(m_fungi)
m_predator <- read.csv("m_predator.csv")
alpha_predator <- alpha_diversity1(m_predator)
m_decomp <- read.csv("m_decomp.csv")
alpha_decomp <- alpha_diversity1(m_decomp)

alpha_diversity <- data.frame(alpha_total[, 1:2], Richness_total = alpha_total$Richness, Richness_fungi = alpha_fungi$Richness, Richness_predator = alpha_predator$Richness, Richness_decomp = alpha_decomp$Richness, 
                              Shannon_total = alpha_total$Shannon, Shannon_fungi = alpha_fungi$Shannon, Shannon_predator = alpha_predator$Shannon, Shannon_decomp = alpha_decomp$Shanno)

alpha_diversity2 <- alpha_diversity %>% melt(id.vars = 1:2, measure.vars = 3:ncol(alpha_diversity), variable.name = "diversity", value.name = "value")
alpha_diversity3 <- summarySE(alpha_diversity2, measurevar = "value", groupvars = c("diversity", "Forest")) %>% mutate(Forest = factor(Forest, levels = c("RBS", "DBF", "EBF", "CBF")))

theme_bw_facet <- function() {
  theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 1, color = "black", fill = NA),
      strip.background = element_rect(linewidth = 1, color = "#999999", fill = "#999999"),
      strip.text = element_text(size = 9, color = "white", face = "bold"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.ticks = element_line(size = 0.6, color = "black"),
      axis.title = element_blank(),
      legend.position = "none"
    )
}

display.brewer.all()
brewer.pal(n = 9, name = 'Set1')
color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

f_alpha_diversity <- ggplot(data = alpha_diversity3, aes(x = Forest, y = value, color = Forest)) +
  facet_wrap(~ diversity, scales = "free", ncol = 4) +
  geom_point(stat = "identity", position = position_dodge(1), shape = 16, size = 4) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge(1), linewidth = 0.6, width = 0, colour = "black") +
  scale_color_manual(values = color) +
  theme_bw_facet(); f_alpha_diversity

# ggsave(filename = "alpha_diversity.pdf", plot = f_alpha_diversity, device = "pdf", width = 16, height = 8.8, units = "cm", dpi = 300)
## (Figure 1A-H) end



analyze_diversity <- function(data, response_var, group_var = "Forest", trans_method = c("none", "log", "sqrt")) {
  require(performance)
  require(agricolae)
  require(rstatix)
  trans_method <- match.arg(trans_method)
  if (trans_method == "log") {
    data[[response_var]] <- log(data[[response_var]] + 1)
  } else if (trans_method == "sqrt") {
    data[[response_var]] <- sqrt(data[[response_var]])
  }
  formula <- as.formula(paste(response_var, "~", group_var))
  fit <- lm(formula, data = data)
  norm_test <- check_normality(fit)
  homog_test <- check_homogeneity(fit)
  norm_p <- as.numeric(norm_test)
  homog_p <- as.numeric(homog_test)
  assumptions_met <- norm_p > 0.05 && homog_p > 0.05
  results <- list(
    assumptions = list(
      normality = norm_test,
      homogeneity = homog_test,
      assumptions_met = assumptions_met
    )
  )
  if (assumptions_met) {
    results$ols <- summary(fit)
    results$posthoc <- LSD.test(fit, trt = group_var, console = FALSE)
  } else {
    results$kruskal <- kruskal.test(formula, data = data)
    results$posthoc <- pairwise_wilcox_test(
      data = data,
      formula = formula,
      p.adjust.method = "bonferroni"
    )
  }
  return(results)
}

analyze_diversity(alpha_diversity, "Richness_total", group_var = "Forest")
analyze_diversity(alpha_diversity, "Richness_fungi", group_var = "Forest")
analyze_diversity(alpha_diversity, "Richness_predator", group_var = "Forest")
analyze_diversity(alpha_diversity, "Richness_decomp", group_var = "Forest")
analyze_diversity(alpha_diversity, "Shannon_total", group_var = "Forest")
analyze_diversity(alpha_diversity, "Shannon_fungi", group_var = "Forest")
analyze_diversity(alpha_diversity, "Shannon_predator", group_var = "Forest")
analyze_diversity(alpha_diversity, "Shannon_decomp", group_var = "Forest")

## (linear regression or non-parametric tests) end



normalize <- function(data, vars) {
  require(tidyverse)
  normalize <- function(v) {
    (v - min(v)) / (max(v) - min(v))
  }
  ret <- dplyr::mutate(data, dplyr::across(vars, normalize, .names = "{.col}.std")) %>%
    dplyr::select(paste0(vars, ".std"))
  ret$meanFunction <- rowSums(ret) / ncol(ret)
  return(ret)
}

allVars_fungi1 <- qw(Richness_fungi, Shannon_fungi) 
alpha_diversity_fungi <- cbind(alpha_diversity, normalize(alpha_diversity, allVars_fungi1))
allVars_predator1 <- qw(Richness_predator, Shannon_predator) 
alpha_diversity_predator <- cbind(alpha_diversity, normalize(alpha_diversity, allVars_predator1))
allVars_decomp1 <- qw(Richness_decomp, Shannon_decomp) 
alpha_diversity_decomp <- cbind(alpha_diversity, normalize(alpha_diversity, allVars_decomp1))

alpha_multifun <- data.frame(Forest = alpha_diversity$Forest,
                             alpha_diversity_fungi = alpha_diversity_fungi$meanFunction, 
                             alpha_diversity_predator = alpha_diversity_predator$meanFunction, 
                             alpha_diversity_decomp = alpha_diversity_decomp$meanFunction)

fit_alpha <- lmer(alpha_diversity_predator ~ alpha_diversity_fungi + (1|Forest), data = alpha_multifun)
summary(fit_alpha)
check_normality(fit_alpha)
car::Anova(fit_alpha, test.statistic = "F")
r2(fit_alpha)
tab_model(fit_alpha)

summary(alpha_multifun)
rdm3 <- ggpredict(fit_alpha, type = "fixed", back.transform = FALSE, terms = "alpha_diversity_fungi[0.0399,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]")

theme_bw_no_facet <- function() {
  theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 1, color = "black", fill = NA),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.ticks = element_line(size = 0.6, color = "black"),
      axis.title = element_blank(),
      legend.position = "none"
    )
}

alpha_multifun$Forest <- factor(alpha_multifun$Forest, levels = c("RBS", "DBF", "EBF", "CBF"))
alpha_pf <- ggplot(alpha_multifun, aes(x = alpha_diversity_fungi, y = alpha_diversity_predator)) +
  geom_point(aes(color = Forest), size = 2.5, shape = 16) +
  scale_color_manual(values = color) +
  geom_line(data = rdm3, aes(x = x, y = predicted), color = "black", linewidth = 0.8) +
  geom_ribbon(data = rdm3, aes(x = x, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.2, fill = "grey50") +
  geom_segment(aes(x = alpha_diversity_fungi, y = alpha_diversity_predator, xend = alpha_diversity_fungi, yend = predict(fit_alpha, re.form = NA)), color = "#999999", linewidth = 0.5, alpha = 0.6) +
  theme_bw_no_facet(); alpha_pf

# ggsave(filename = "alpha_pf.pdf", plot = alpha_pf, device = "pdf", width = 4.5, height = 4.5, units = "cm", dpi = 300)
## (Figure 1I) end



fit_alpha1 <- lmer(alpha_diversity_predator ~ alpha_diversity_decomp + (1|Forest), data = alpha_multifun)
summary(fit_alpha1)
check_normality(fit_alpha1)
car::Anova(fit_alpha1, test.statistic = "F")
r2(fit_alpha1)
tab_model(fit_alpha1)

summary(alpha_multifun)
rdm4 <- ggpredict(fit_alpha1, type = "fixed", back.transform = FALSE, terms = "alpha_diversity_decomp[0.0000,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]")

alpha_pf1 <- ggplot(alpha_multifun, aes(x = alpha_diversity_decomp, y = alpha_diversity_predator)) +
  geom_point(aes(color = Forest), size = 2.5, shape = 16) +
  scale_color_manual(values = color) +
  geom_line(data = rdm4, aes(x = x, y = predicted), color = "black", linewidth = 0.8) +
  geom_ribbon(data = rdm4, aes(x = x, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.2, fill = "grey50") +
  geom_segment(aes(x = alpha_diversity_decomp, y = alpha_diversity_predator, xend = alpha_diversity_decomp, yend = predict(fit_alpha1, re.form = NA)), color = "#999999", linewidth = 0.5, alpha = 0.6) +
  theme_bw_no_facet(); alpha_pf1

# ggsave(filename = "alpha_pf1.pdf", plot = alpha_pf1, device = "pdf", width = 4.5, height = 4.5, units = "cm", dpi = 300)
## (Figure 1j) end



correlation <- function(community, environment) {
  require(psych)
  require(pheatmap)
  cor_result <- corr.test(community, environment, method = "pearson", adjust = "none")
  r_mat <- cor_result$r
  p_mat <- cor_result$p
  sig_symbols <- c("***" = 0.001, "**" = 0.01, "*" = 0.05)
  sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  for (i in order(sig_symbols, decreasing = FALSE)) {
    sig_mat[p_mat < sig_symbols[i] & sig_mat == ""] <- names(sig_symbols)[i]
  }
  r_mat <- scales::rescale(r_mat, to = c(-1, 1))
  heatmap <- pheatmap(
    r_mat,
    scale = "none",
    color = colorRampPalette(c("#548235", "white", "#CE4368"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = F, cluster_cols = T,
    clustering_method = "average",
    cutree_rows = 4, cutree_cols = 3,
    border_color = "grey100",
    cellwidth = 13, cellheight = 13,
    fontsize_row = 9, fontsize_col = 9,
    display_numbers = sig_mat,
    fontsize_number = 9, number_color = "black",
    legend_breaks = c(-1, -0.5, 0, 0.5, 1), legend_labels = c(-1, -0.5, 0, 0.5, 1),
    treeheight_row = 25, treeheight_col = 15
  )
  return(list(
    correlation_matrix = r_mat,
    p_value_matrix = p_mat,
    heatmap = heatmap
  ))
}

alpha_corr <- correlation(alpha_multifun[,-1], alpha_diversity[,-c(1:2, 3, 7)])

# ggsave(filename = "alpha_corr.pdf", plot = alpha_corr$heatmap, device = "pdf", width = 10, height = 11, units = "cm", dpi = 300)
## (Figure 1k) end






#### 2 Multitrophic diversity within mite communities - beta diversity (Figure 2) ####
mite_nmds <- function(data, Forest = "Forest") {
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(purrr)
  nmds_result <- metaMDS(data[, -(1:2)], distance = "bray", k = 2)
  nmds_points <- data.frame(NMDS1 = nmds_result$points[,1], NMDS2 = nmds_result$points[,2], data[, 1:2])
  centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Forest, data = nmds_points, FUN = mean)
  set.seed(123)
  anosim_overall <- anosim(data[, -(1:2)], grouping = data$Forest, permutations = 999, distance = "bray")
  pairwise_results <- combn(unique(data$Forest), 2, simplify = FALSE) %>% 
    map_dfr(function(pair) {
      subset_data <- filter(data, Forest %in% pair)
      anosim_res <- anosim(subset_data[, 3:ncol(subset_data)], grouping = subset_data$Forest, permutations = 999, distance = "bray")
      tibble(Group1 = pair[1], Group2 = pair[2], R = anosim_res$statistic, p_value = anosim_res$signif)
    }) %>%
    arrange(desc(R)) %>%
    mutate(Significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ "ns"))
  color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  names(color) <- c("RBS", "DBF", "EBF", "CBF")
  Forest = factor(Forest, levels = c("RBS", "DBF", "EBF", "CBF"))
  nmds_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Forest)) +
    geom_point(shape = 20) +
    scale_color_manual(values = color) +
    geom_point(data = centroids, aes(x = NMDS1, y = NMDS2, color = Forest), shape = 3, size = 2, stroke = 1) +
    geom_text(data = centroids, aes(x = NMDS1, y = NMDS2, label = Forest), vjust = -1, size = 3, color = "black") +
    stat_ellipse(geom = "polygon", level = 0.5, linetype = 1, size = 0.5, aes(color = Forest), fill = NA, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5, color = "grey", alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5, color = "grey", alpha = 0.3) +
    theme_bw() +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank(), panel.border = element_rect(linewidth = 1, color = "black", fill = NA), axis.text.x = element_text(size = 9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.ticks = element_line(size = 0.6, color = "black"), axis.title = element_blank()) +
    guides(colour = "none") +
    xlab("") + ylab("")
  return(list(nmds_stress = nmds_result$stress, nmds_data = nmds_points, centroids = centroids, anosim_summary = anosim_overall, pairwise_comparisons = pairwise_results, nmds_plot = nmds_plot))
}

nmds_m_total <- mite_nmds(m_total)
nmds_m_fungi <- mite_nmds(m_fungi)
nmds_m_predator <- mite_nmds(m_predator)
nmds_m_decomp <- mite_nmds(m_decomp)

nmds_p <- plot_grid(nmds_m_total$nmds_plot, nmds_m_fungi$nmds_plot, nmds_m_predator$nmds_plot, nmds_m_decomp$nmds_plot, ncol = 4, nrow = 1, align = "hv"); nmds_p

# ggsave(filename = "beta_diversity.pdf", plot = nmds_p, device = "pdf", width = 21.5, height = 4.9, units = "cm", dpi = 300)
## (Figure 2A-D) end



dist_bray_total <- beta.div.comp(m_total[, -c(1,2)], coef = "J", quant = F)
dist_bray_fungi <- beta.div.comp(m_fungi[, -c(1,2)], coef = "J", quant = F)
dist_bray_predator <- beta.div.comp(m_predator[, -c(1,2)], coef = "J", quant = F)
dist_bray_decomp <- beta.div.comp(m_decomp[, -c(1,2)], coef = "J", quant = F)

beta_part <- dist.3col(dist_bray_total$D)
colnames(beta_part)[3] <- "dist_total"
beta_part$dist_fungi <- dist.3col(dist_bray_fungi$D)[,3]
beta_part$dist_predator <- dist.3col(dist_bray_predator$D)[,3]
beta_part$dist_decomp <- dist.3col(dist_bray_decomp$D)[,3]

beta_part$repl_total <- dist.3col(dist_bray_total$repl)[,3]
beta_part$repl_fungi <- dist.3col(dist_bray_fungi$repl)[,3]
beta_part$repl_predator <- dist.3col(dist_bray_predator$repl)[,3]
beta_part$repl_decomp <- dist.3col(dist_bray_decomp$repl)[,3]
beta_part$riff_total <- dist.3col(dist_bray_total$rich)[,3]

beta_part$riff_fungi <- dist.3col(dist_bray_fungi$rich)[,3]
beta_part$riff_predator <- dist.3col(dist_bray_predator$rich)[,3]
beta_part$riff_decomp <- dist.3col(dist_bray_decomp$rich)[,3]

beta_part_mean <- as.data.frame(t(colMeans(beta_part[,-(1:2)], na.rm = TRUE)))

theme_bvbw_c <- function(base_size = 10, base_family = "sans") {
  col.T = "#A65628"
  col.L = rgb2hex(0, 0, 0)
  col.R = "#999999"
  theme_ggtern(base_size, base_family) %+replace%
    theme(
      tern.plot.background = element_rect(fill = "white", linewidth = NA, color = NA),
      tern.panel.background = element_rect(fill = "white"),
      tern.panel.grid.major = element_line(linewidth = 0.2, linetype = 2),
      tern.panel.grid.major.T = element_line(color = col.T),
      tern.panel.grid.major.L = element_line(color = col.L),
      tern.panel.grid.major.R = element_line(color = col.R),
      tern.panel.grid.minor = element_blank(),
      tern.axis.arrow.show = TRUE,
      tern.axis.arrow.T = element_line(color = col.T),
      tern.axis.arrow.L = element_line(color = col.L),
      tern.axis.arrow.R = element_line(color = col.R),
      tern.axis.title.T = element_blank(),
      tern.axis.title.L = element_blank(),
      tern.axis.title.R = element_blank(),
      tern.axis.arrow.text.T = element_text(color = col.T),
      tern.axis.arrow.text.L = element_text(color = col.L),
      tern.axis.arrow.text.R = element_text(color = col.R),
      tern.axis.line.T = element_line(color = col.T, linewidth = 0.45),
      tern.axis.line.L = element_line(color = col.L, linewidth = 0.45),
      tern.axis.line.R = element_line(color = col.R, linewidth = 0.45),
      tern.axis.ticks.major.T = element_line(color = col.T, linewidth = 1),
      tern.axis.ticks.major.L = element_line(color = col.L, linewidth = 1),
      tern.axis.ticks.major.R = element_line(color = col.R, linewidth = 1),
      tern.axis.ticks.minor.T = element_blank(),
      tern.axis.ticks.minor.L = element_blank(),
      tern.axis.ticks.minor.R = element_blank(),
      legend.position = "none"
    )
}

beta_p_total <- ggtern(beta_part, aes(1 - dist_total, repl_total, riff_total)) +
  geom_point(alpha = 0.6, size = 1.5, color = "grey70", shape = 16) +
  stat_density_tern(geom = 'polygon', bdl.val = NA, aes(fill = after_stat(level), alpha = after_stat(level))) +
  scale_fill_gradient(low = "#f3f421", high = "#007410") +
  geom_point(data = beta_part_mean, color = "#E41A1C", size = 3, shape = 20) +
  labs(x = "similarity(1-βJaccard)", y = "Turnover(βRepl)", z = "Nestedness(βRichDif)", title = "") +
  theme_bvbw_c(); beta_p_total

beta_p_fungi <- ggtern(beta_part, aes(1 - dist_fungi, repl_fungi, riff_fungi)) +
  geom_point(alpha = 0.6, size = 1.5, color = "#984EA3", shape = 16) +
  stat_density_tern(geom = 'polygon', bdl.val = NA, aes(fill = after_stat(level), alpha = after_stat(level))) +
  scale_fill_gradient(low = "#f3f421", high = "#007410") +
  geom_point(data = beta_part_mean, color = "#E41A1C", size = 3, shape = 20) +
  labs(x = "similarity(1-βJaccard)", y = "Turnover(βRepl)", z = "Nestedness(βRichDif)", title = "") +
  theme_bvbw_c(); beta_p_fungi

beta_p_predator <- ggtern(beta_part, aes(1 - dist_predator, repl_predator, riff_predator)) +
  geom_point(alpha = 0.6, size = 1.5, color = "#FF7F00", shape = 16) +
  stat_density_tern(geom = 'polygon', bdl.val = NA, aes(fill = after_stat(level), alpha = after_stat(level))) +
  scale_fill_gradient(low = "#f3f421", high = "#007410") +
  geom_point(data = beta_part_mean, color = "#E41A1C", size = 3, shape = 20) +
  labs(x = "similarity(1-βJaccard)", y = "Turnover(βRepl)", z = "Nestedness(βRichDif)", title = "") +
  theme_bvbw_c(); beta_p_predator

beta_p_decomp <- ggtern(beta_part, aes(1 - dist_decomp, repl_decomp, riff_decomp)) +
  geom_point(alpha = 0.6, size = 1.5, color = "#377EB8", shape = 16) +
  stat_density_tern(geom = 'polygon', bdl.val = NA, aes(fill = after_stat(level), alpha = after_stat(level))) +
  scale_fill_gradient(low = "#f3f421", high = "#007410") +
  geom_point(data = beta_part_mean, color = "#E41A1C", size = 3, shape = 20) +
  labs(x = "similarity(1-βJaccard)", y = "Turnover(βRepl)", z = "Nestedness(βRichDif)", title = "") +
  theme_bvbw_c(); beta_p_decomp

# ggsave(filename = "beta_p_decomp.pdf", plot = beta_p_decomp, device = "pdf", width = 6, height = 6, units = "cm", dpi = 300)
## (Figure 2E-F) end






#### 3 Assembly processes of mite communities across trophic groups (Figure 3) ####
tnst_total <- tNST(comm = m_total[,-c(1:2)], group = m_total[,c(1:2)], dist.method = 'jaccard', abundance.weighted = F, null.model = 'PF', rand = 1000, nworker = 8)
tnst_fungi <- tNST(comm = m_fungi[,-c(1:2)], group = m_fungi[,c(1:2)], dist.method = 'jaccard', abundance.weighted = F, null.model = 'PF', rand = 1000, nworker = 8)
tnst_predator <- tNST(comm = m_predator[,-c(1:2)], group = m_predator[,c(1:2)], dist.method = 'jaccard', abundance.weighted = F, null.model = 'PF', rand = 1000, nworker = 8)
tnst_decomp <- tNST(comm = m_decomp[,-c(1:2)], group = m_decomp[,c(1:2)], dist.method = 'jaccard', abundance.weighted = F, null.model = 'PF', rand = 1000, nworker = 8)

tnst1 <- data.frame(Forest = tnst_total$index.pair.grp[,1],
                    tnst_total$index.pair.grp[,2:3],
                    tnst_total = tnst_total$index.pair.grp$MST.ij.jaccard,
                    tnst_fungi = tnst_fungi$index.pair.grp$MST.ij.jaccard,
                    tnst_predator = tnst_predator$index.pair.grp$MST.ij.jaccard,
                    tnst_decomp = tnst_decomp$index.pair.grp$MST.ij.jaccard)

analyze_diversity(tnst1, "tnst_total", group_var = "Forest", trans_method = "none")
analyze_diversity(tnst1, "tnst_fungi", group_var = "Forest", trans_method = "none")
analyze_diversity(tnst1, "tnst_predator", group_var = "Forest", trans_method = "none")
analyze_diversity(tnst1, "tnst_decomp", group_var = "Forest", trans_method = "none")

tnst11 <- tnst1 %>%
  tidyr::pivot_longer(cols = 4:ncol(tnst1), names_to = "group", values_to = "value") %>%
  mutate(Forest = factor(Forest, levels = c("RBS", "DBF", "EBF", "CBF")))

assembly <- ggplot(data = tnst11, aes(x = Forest, y = value, fill = Forest)) +
  facet_wrap(~ factor(group, levels = c("tnst_total", "tnst_fungi", "tnst_predator", "tnst_decomp")), scales = "free", ncol = 4) +
  geom_boxplot(width = 0.6, linewidth = 0.6, alpha = 0.6, fatten = 1, outlier.shape = NA) +
  geom_jitter(aes(colour = Forest), position = position_jitterdodge(jitter.width = 1), shape = 16, alpha = 1, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 25, size = 1, color = "black", fill = "black") +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  xlab("") + ylab("Normalized stochasticity ratio") +
  theme_bw_facet(); assembly

# ggsave(filename = "community assembly.pdf", plot = assembly, device = "pdf", width = 19, height = 5, units = "cm", dpi = 300)
## (Figure 3A-D) end



community_mean_niche <- function(data, method = "levins") {
  require(spaa)
  B <- niche.width(data, method = method)
  B_com <- numeric(nrow(data))
  for(i in seq_len(nrow(data))) {
    present_species <- data[i, ] > 0
    species_names <- names(data[i, present_species])
    if(length(species_names) > 0) {
      B_com[i] <- mean(as.numeric(B[, species_names]), na.rm = TRUE)
    } else {
      B_com[i] <- NA  
    }
  }
  return(list(
    species_niche = data.frame(species_id = colnames(B), 
                               species_niche = as.numeric(B[1, ])),
    community_mean_niche = data.frame(sample_id = rownames(data),
                                      mean_niche_width = B_com),
    overall_mean_niche = mean(as.numeric(B), na.rm = TRUE)
  ))
}

niche_total <- community_mean_niche(m_total[,-c(1:2)])
niche_fungi <- community_mean_niche(m_fungi[,-c(1:2)])
niche_predator <- community_mean_niche(m_predator[,-c(1:2)])
niche_decomp <- community_mean_niche(m_decomp[,-c(1:2)])

niche1 <- data.frame(Forest = m_total$Forest,
                     niche_total = niche_total$community_mean_niche$mean_niche_width,
                     niche_fungi = niche_fungi$community_mean_niche$mean_niche_width,
                     niche_predator = niche_predator$community_mean_niche$mean_niche_width,
                     niche_decomp = niche_decomp$community_mean_niche$mean_niche_width)

analyze_diversity(niche1, "niche_fungi", group_var = "Forest", trans_method = "none")
analyze_diversity(niche1, "niche_predator", group_var = "Forest", trans_method = "none")
analyze_diversity(niche1, "niche_decomp", group_var = "Forest", trans_method = "none")

niche11 <- niche1 %>%
  tidyr::pivot_longer(cols = 2:ncol(niche1), names_to = "group", values_to = "value") %>%
  mutate(Forest = factor(Forest, levels = c("RBS", "DBF", "EBF", "CBF")))

niche_breadth <- ggplot(data = niche11, aes(x = Forest, y = value, fill = Forest)) +
  facet_wrap(~ factor(group, levels = c("niche_total", "niche_fungi", "niche_predator", "niche_decomp")), scales = "free", ncol = 4) +
  geom_boxplot(width = 0.6, linewidth = 0.6, alpha = 0.6, fatten = 1, outlier.shape = NA) +
  geom_jitter(aes(colour = Forest), position = position_jitterdodge(jitter.width = 0.6), shape = 16, alpha = 1, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 25, size = 1, color = "black", fill = "black") +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  theme_bw_facet() +
  xlab("Forest") + ylab("Habitat niche breadth (Bcom)"); niche_breadth

# ggsave(filename = "niche breadth.pdf", plot = niche_breadth, device = "pdf", width = 19, height = 5, units = "cm", dpi = 300)
## (Figure 3E-H) end






#### 4 Assembly processes of mite communities across trophic groups (Figure 4) ####
## Although the ecological network analysis was performed using the following R code, the co-occurrence network graph was visualized in the Gephi platform (https://gephi.org).
net <- function(data, p_threshold = 0.05, r_threshold = 0.5) {
  data <- as.matrix(data)
  require(Hmisc)   
  require(igraph)
  cor_result <- rcorr(data, type = "spearman")
  p_mite <- cor_result$P
  r_mite <- cor_result$r
  r_mite[p_mite > p_threshold | abs(r_mite) < r_threshold] = 0
  diag(r_mite) <- 0 
  g <- graph_from_adjacency_matrix(
    r_mite,
    weighted = TRUE,
    mode = "undirected",
  )
  g <- simplify(g)
  return(g)
}

# m_total_n[, c("X.21..Camisia", "X.27..Archegozetes", "X.73..Phauloppia")] # Phytophagous and lichen feeders were excluded from this study as they comprised only a few species.
g_m_total <- net(m_total[,-c(1:2, 67, 73, 119)])
g_m_fungi <- net(m_fungi[,-c(1:2)])
g_m_predator <- net(m_predator[,-c(1:2)])
g_m_decomp <- net(m_decomp[,-c(1:2)])

process_net <- function(g) {
  require(igraph)
  require(dplyr)
  g_abs <- g
  if (!is.null(E(g_abs)$weight)) {
    E(g_abs)$weight_abs <- abs(E(g_abs)$weight)  
    E(g_abs)$cor_sign <- sign(E(g_abs)$weight)
  }
  edge <- igraph::as_data_frame(g, what = "edges") %>%
    mutate(weight_abs = abs(weight),
           cor_sign = sign(weight)
    )
  node <- data.frame(id = names(V(g)),
                     degree = degree(g),
                     weight_degree = strength(g),
                     betweenness_centrality = betweenness(g_abs, weights = E(g_abs)$weight_abs),
                     closeness_centrality = closeness(g_abs, weights = E(g_abs)$weight_abs),
                     eigenvector_centrality = eigen_centrality(g_abs, weights = E(g_abs)$weight_abs)$vector
  )
  basic <- data.frame(
    nodes = length(V(g)),
    edges = length(E(g)),
    positive = sum(E(g)$weight>0),
    negative = sum(E(g)$weight<0),
    average_degree = mean(degree(g)),
    network_density = edge_density(g),
    network_diameter = diameter(g_abs, weights = E(g_abs)$weight_abs),
    average_path_length = mean_distance(g_abs, weights = E(g_abs)$weight_abs),
    clustering_coefficient = transitivity(g),
    betweenness_centralization = centr_betw(g)$centralization
  )
  set.seed(123)
  comm <- cluster_fast_greedy(g_abs, weights = E(g_abs)$weight_abs)
  community <- list(
    communities = comm,
    modularity = modularity(comm),
    mod_size = sizes(comm)
  )
  node$membership <-  membership(comm)
  plot(comm, g_abs,
       layout = layout_with_fr(g_abs, weights = E(g_abs)$weight_abs),
       vertex.size = degree(g_abs),
       edge.width = 0.5,
       edge.color = adjustcolor("gray60", alpha.f = 0.5),
       mark.groups = communities(comm)
  )
  return(list(
    edges = edge,
    nodes = node,
    basic = basic,
    community = community,
    original_graph = g,
    abs_weight_graph = g_abs
  ))
}

process_net(g_m_total)
process_net(g_m_fungi)
process_net(g_m_predator)
process_net(g_m_decomp)

m_total_edges <- (process_net(g_m_total)$edges) 
m_total_nodes <- (process_net(g_m_total)$nodes)
m_total_topological <- (process_net(g_m_total)$basic)
# write.csv(m_total_edges, "m_total_edges.csv") # Save the edge list to a CSV file for network visualization in Gephi.
# write.csv(m_total_nodes, "m_total_nodes.csv") # Save the node list to a CSV file for network visualization in Gephi.
# write.graph(process_net(g_m_total)$abs_weight_graph, "g_m_total.graphml", "graphml") # Save the network in GraphML format for immediate visualization in Gephi.

m_fungi_edges <- (process_net(g_m_fungi)$edges)
m_fungi_nodes <- (process_net(g_m_fungi)$nodes)
m_fungi_topological <- (process_net(g_m_fungi)$basic)
# write.csv(m_fungi_edges, "m_fungi_edges.csv")
# write.csv(m_fungi_nodes, "m_fungi_nodes.csv")
# write.graph(process_net(g_m_predator)$abs_weight_graph, "g_m_predator.graphml", "graphml")

m_predator_edges <- (process_net(g_m_predator)$edges)
m_predator_nodes <- (process_net(g_m_predator)$nodes)
# write.csv(m_predator_edges, "m_predator_edges.csv")
# write.csv(m_predator_nodes, "m_predator_nodes.csv")
# write.graph(process_net(g_m_fungi)$abs_weight_graph, "g_m_fungi.graphml", "graphml")

m_decomp_edges <- (process_net(g_m_decomp)$edges)
m_decomp_nodes <- (process_net(g_m_decomp)$nodes)
# write.csv(m_decomp_edges, "m_decomp_edges.csv")
# write.csv(m_decomp_nodes, "m_decomp_nodes.csv")
# write.graph(process_net(g_m_decomp)$abs_weight_graph, "g_m_decomp.graphml", "graphml")

## (Figure 4A) end



m_RBS <- m_total[1:6,-c(1:2, 67, 73, 119)]
m_RBS <- m_RBS[, colSums(m_RBS) > 0]
g_m_RBS <- net(m_RBS)
process_net(g_m_RBS)
# write.graph(process_net(g_m_RBS)$abs_weight_graph, "g_m_RBS.graphml", "graphml")
m_RBS_edges <- (process_net(g_m_RBS)$edges)
m_RBS_nodes <- (process_net(g_m_RBS)$nodes)
# write.csv(m_RBS_edges, "m_RBS_edges.csv")
# write.csv(m_RBS_nodes, "m_RBS_nodes.csv")

m_DBF <- m_total[7:12,-c(1:2, 67, 73, 119)]
m_DBF <- m_DBF [, colSums(m_DBF ) > 0]
names(m_DBF)
g_m_DBF <- net(m_DBF)
process_net(g_m_DBF)
# write.graph(process_net(g_m_DBF)$abs_weight_graph, "g_m_DBF.graphml", "graphml")
m_DBF_edges <- (process_net(g_m_DBF)$edges)
m_DBF_nodes <- (process_net(g_m_DBF)$nodes)
# write.csv(m_DBF_edges, "m_DBF_edges.csv")
# write.csv(m_DBF_nodes, "m_DBF_nodes.csv")

m_EBF <- m_total[13:18,-c(1:2, 67, 73, 119)]
m_EBF <- m_EBF [, colSums(m_EBF ) > 0]
names(m_EBF)
g_m_EBF <- net(m_EBF)
process_net(g_m_EBF)
# write.graph(process_net(g_m_EBF)$abs_weight_graph, "g_m_EBF.graphml", "graphml")
m_EBF_edges <- (process_net(g_m_EBF)$edges)
m_EBF_nodes <- (process_net(g_m_EBF)$nodes)
# write.csv(m_EBF_edges, "m_EBF_edges.csv")
# write.csv(m_EBF_nodes, "m_EBF_nodes.csv")

m_CBF <- m_total[19:24,-c(1:2, 67, 73, 119)]
m_CBF <- m_CBF [, colSums(m_CBF ) > 0]
names(m_CBF)
g_m_CBF <- net(m_CBF)
process_net(g_m_CBF)
# write.graph(process_net(g_m_CBF)$abs_weight_graph, "g_m_CBF.graphml", "graphml")
m_CBF_edges <- (process_net(g_m_CBF)$edges)
m_CBF_nodes <- (process_net(g_m_CBF)$nodes)
# write.csv(m_CBF_edges, "m_CBF_edges.csv")
# write.csv(m_CBF_nodes, "m_CBF_nodes.csv")

## (Figure 4D-G) end


sample_g_m_fungi <- extract_sample_net(g_m_fungi, t(m_fungi[,-c(1:2)]))
sample_g_m_predator <- extract_sample_net(g_m_predator, t(m_predator[,-c(1:2)]))
sample_g_m_decomp <- extract_sample_net(g_m_decomp, t(m_decomp[,-c(1:2)]))
## Topological parameters of the subnetworks were analyzed to construct the network complexity of the multitrophic mite groups.

allVars_mite <- qw(Node_number, Edge_number, Average_path_length, Average_degree, Average_weighted_degree, Diameter, Natural_connectivity)
sample_g_m_fungi1 <- data.frame(Forest = alpha_diversity$Forest, Node_number = sample_g_m_fungi$Node_number, Edge_number = sample_g_m_fungi$Edge_number, Average_path_length = -1 * sample_g_m_fungi$Average_path_length, Average_degree = sample_g_m_fungi$Average_degree, Average_weighted_degree = sample_g_m_fungi$Average_weighted_degree, Diameter = -1 * sample_g_m_fungi$Diameter, Natural_connectivity = sample_g_m_fungi$Natural_connectivity)
network_complexity_fungi <- cbind(sample_g_m_fungi1, normalize(sample_g_m_fungi1, allVars_mite))

sample_g_m_predator1 <- data.frame(Forest = alpha_diversity$Forest, Node_number = sample_g_m_predator$Node_number, Edge_number = sample_g_m_predator$Edge_number, Average_path_length = -1 * sample_g_m_predator$Average_path_length, Average_degree = sample_g_m_predator$Average_degree, Average_weighted_degree = sample_g_m_predator$Average_weighted_degree, Diameter = -1 * sample_g_m_predator$Diameter, Natural_connectivity = sample_g_m_predator$Natural_connectivity)
network_complexity_predator <- cbind(sample_g_m_predator1, normalize(sample_g_m_predator1, allVars_mite))

sample_g_m_decomp1 <- data.frame(Forest = alpha_diversity$Forest, Node_number = sample_g_m_decomp$Node_number, Edge_number = sample_g_m_decomp$Edge_number, Average_path_length = -1 * sample_g_m_decomp$Average_path_length, Average_degree = sample_g_m_decomp$Average_degree, Average_weighted_degree = sample_g_m_decomp$Average_weighted_degree, Diameter = -1 * sample_g_m_decomp$Diameter, Natural_connectivity = sample_g_m_decomp$Natural_connectivity)
network_complexity_decomp <- cbind(sample_g_m_decomp1, normalize(sample_g_m_decomp1, allVars_mite))

complexity_multifun <- data.frame(Forest = alpha_diversity$Forest, network_complexity_fungi = network_complexity_fungi$meanFunction, network_complexity_predator = network_complexity_predator$meanFunction, network_complexity_decomp = network_complexity_decomp$meanFunction)

fit_complexity <- lmer(network_complexity_predator ~ network_complexity_fungi + (1|Forest), data = complexity_multifun)
summary(fit_complexity)
check_normality(fit_complexity)
car::Anova(fit_complexity, test.statistic = "F")
r2(fit_complexity)
tab_model(fit_complexity)

summary(complexity_multifun)
rdm1 <- ggpredict(fit_complexity, type = "fixed", back.transform = FALSE, terms = "network_complexity_fungi[0.1950,0.2,0.3,0.4,0.5,0.6,0.7,0.7661]")
complexity_pf <- ggplot(complexity_multifun, aes(x = network_complexity_fungi, y = network_complexity_predator)) +
  geom_point(aes(color = Forest), size = 2.5, shape = 16) +
  scale_color_manual(values = color) +
  geom_line(data = rdm1, aes(x = x, y = predicted), color = "black", linewidth = 0.8) +
  geom_ribbon(data = rdm1, aes(x = x, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.2, fill = "grey50") +
  geom_segment(aes(x = network_complexity_fungi, y = network_complexity_predator, xend = network_complexity_fungi, yend = predict(fit_complexity, re.form = NA)), color = "#999999", linewidth = 0.5, alpha = 0.6) +
  theme_bw_no_facet(); complexity_pf

# ggsave(filename = "complexity_pf.pdf", plot = complexity_pf, device = "pdf", width = 4.5, height = 4.5, units = "cm", dpi = 300)
## (Figure 4B) end



fit_complexity2 <- lmer(network_complexity_predator ~ network_complexity_decomp + (1|Forest), data = complexity_multifun)
summary(fit_complexity2)
check_normality(fit_complexity2)
car::Anova(fit_complexity2, test.statistic = "F")
r2(fit_complexity2)
tab_model(fit_complexity2)

summary(complexity_multifun)
rdm2 <- ggpredict(fit_complexity2, type = "fixed", back.transform = FALSE, terms = "network_complexity_decomp[0.3535,0.4,0.5,0.6,0.7,0.7490]")

complexity_pf1 <- ggplot(complexity_multifun, aes(x = network_complexity_decomp, y = network_complexity_predator)) +
  geom_point(aes(color = Forest), size = 2.5, shape = 16) +
  scale_color_manual(values = color) +
  geom_line(data = rdm2, aes(x = x, y = predicted), color = "black", linewidth = 0.8) +
  geom_ribbon(data = rdm2, aes(x = x, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.2, fill = "grey50") +
  geom_segment(aes(x = network_complexity_decomp, y = network_complexity_predator, xend = network_complexity_decomp, yend = predict(fit_complexity2, re.form = NA)), color = "#999999", linewidth = 0.5, alpha = 0.6) +
  theme_bw_no_facet(); complexity_pf1

# ggsave(filename = "complexity_pf1.pdf", plot = complexity_pf1, device = "pdf", width = 4.5, height = 4.5, units = "cm", dpi = 300) 
## (Figure 4C) end






#### 5 Relative importance of different soil variables (Figure 5) ####
beta_multifun <- data.frame(Forest = alpha_diversity$Forest, beta_diversity_fungi = nmds_m_fungi$nmds_data[,1], beta_diversity_predator = nmds_m_predator$nmds_data[,1], beta_diversity_decomp = nmds_m_decomp$nmds_data[,1])

envir <- read.csv("environment.csv")
diversity_multifun1 <- data.frame(envir[,1:2], scale(envir[,6:17]), alpha_multifun[,-1], beta_multifun[,-1], complexity_multifun[,-1])

fit11 <- lm(alpha_diversity_fungi ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1) # Full model
check_normality(fit11)
summary(fit11)

alias(fit11)
check_collinearity(lm(alpha_diversity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)) # we iteratively removed the variable with the highest VIF from the set (BD, C.P, SOM, N.P)

fit12 <- lm(alpha_diversity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1) # Reduced model (VIf<10)
check_normality(fit12)
summary(fit12)

options(na.action = 'na.fail')
fit13 <- dredge(fit12) # Model averaging (ΔAIC < 2)
subset(fit13, delta < 2) # 1

fit14 <- lm(alpha_diversity_fungi ~ SWC + TN + TK + C.N, data = diversity_multifun1) # Final model
check_normality(fit14)
summary(fit14)
model_performance(fit11); model_performance(fit12); model_performance(fit14)
r2(fit11); r2(fit12); r2(fit14)
tab_model(fit14)

fit_data11 <- summary(fit14)$coefficients[, 1] %>% as.data.frame() %>% 
  mutate(lower = confint(fit14, full = TRUE)[, 1], upper = confint(fit14, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "TK", "TN", "C.N")),
         type = case_when(variable == "SWC" ~ "Soil structure",
                          variable %in% c("TN", "TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))

f_alpha_diversity_fungi <- ggplot(fit_data11) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_alpha_diversity_fungi

## (Figure 5A) end


fit31 <- lm(alpha_diversity_predator ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit31)
summary(fit31)
vif(fit31)
alias(fit31)
vif(lm(alpha_diversity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit32 <- lm(alpha_diversity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit32)
summary(fit32)
fit33 <- dredge(fit32)
subset(fit33, delta < 2) # 7
best33 <- model.avg(fit33, subset = delta < 2) # full average
summary(best33)
sw(best33)
fit34 <- lm(alpha_diversity_predator ~ SWC + TK + Por + pH + C.N, data = diversity_multifun1)
check_normality(fit34)
summary(fit34)
AICc(fit31, fit32, fit34)
r2(fit31); r2(fit32); r2(fit34)
fit_data21 <- summary(best33)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best33, full = TRUE)[, 1], upper = confint(best33, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "Por", "TK", "pH", "C.N")),
         type = case_when(variable %in% c("SWC", "Por") ~ "Soil structure",
                          variable %in% c("pH", "TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_alpha_diversity_predator <- ggplot(fit_data21) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_alpha_diversity_predator

## (Figure 5B) end



fit41 <- lm(alpha_diversity_decomp ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit41)
summary(fit41)
alias(fit41)
vif(lm(alpha_diversity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit42 <- lm(alpha_diversity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit42)
summary(fit42)
fit43 <- dredge(fit42)
subset(fit43, delta < 2) # 3
best43 <- model.avg(fit43, subset = delta < 2)
summary(best43)
sw(best43)
fit44 <- lm(alpha_diversity_decomp ~ TK + pH + C.N, data = diversity_multifun1)
check_normality(fit44)
AICc(fit41, fit32, fit44)
summary(fit44)
r2(fit41); r2(fit42); r2(fit44)
fit_data31 <- summary(best43)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best43, full = TRUE)[, 1], upper = confint(best43, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("pH", "TK", "C.N")),
         type = case_when(variable %in% c("pH", "TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_alpha_diversity_decomp <- ggplot(fit_data31) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_alpha_diversity_decomp

## (Figure 5C) end



fit51 <- lm(network_complexity_fungi ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit51)
summary(fit51)
alias(fit51)
vif(lm(network_complexity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit52 <- lm(network_complexity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit52)
summary(fit52)
fit53 <- dredge(fit52)
subset(fit53, delta < 2)
best53 <- model.avg(fit53, subset = delta < 2)
summary(best53)
sw(best53)
fit54 <- lm(network_complexity_fungi ~ SWC + NWC + TN + TK + TP + C.N, data = diversity_multifun1)
check_normality(fit54)
AICc(fit51, fit52, fit54)
summary(fit54)
r2(fit51); r2(fit52); r2(fit54)
fit_data41 <- summary(best53)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best53, full = TRUE)[, 1], upper = confint(best53, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = 1) %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "NWC", "TN", "TK", "TP", "C.N")),
         type = case_when(variable %in% c("SWC", "NWC") ~ "Soil structure",
                          variable %in% c("TN", "TK", "TP") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_network_complexity_fungi <- ggplot(fit_data41) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_network_complexity_fungi

## (Figure 5G) end


fit61 <- lm(network_complexity_predator ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit61)
summary(fit61)
alias(fit61)
vif(lm(network_complexity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit62 <- lm(network_complexity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit62)
summary(fit62)
fit63 <- dredge(fit62)
subset(fit63, delta < 2)
best63 <- model.avg(fit63, subset = delta < 2)
summary(best63)
sw(best63)
fit64 <- lm(network_complexity_predator ~ SWC + TK + NWC + C.N + pH, data = diversity_multifun1)
check_normality(fit64)
AICc(fit61, fit62, fit64)
summary(fit64)
r2(fit61); r2(fit62); r2(fit64)
fit_data51 <- summary(best63)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best63, full = TRUE)[, 1], upper = confint(best63, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "NWC", "TK", "pH", "C.N")),
         type = case_when(variable %in% c("SWC", "NWC") ~ "Soil structure",
                          variable %in% c("TK", "pH") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_network_complexity_predator <- ggplot(fit_data51) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_network_complexity_predator

## (Figure 5H) end


fit71 <- lm(network_complexity_decomp ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit71)
summary(fit71)
alias(fit71)
vif(lm(network_complexity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit72 <- lm(network_complexity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit72)
summary(fit72)
fit73 <- dredge(fit72)
subset(fit73, delta < 2)
best73 <- model.avg(fit73, subset = delta < 2)
summary(best73)
sw(best73)
fit74 <- lm(network_complexity_decomp ~ C.N + TK, data = diversity_multifun1)
check_normality(fit74)
AICc(fit71, fit72, fit74)
summary(fit74)
r2(fit71); r2(fit72); r2(fit74)
fit_data61 <- summary(best73)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best73, full = TRUE)[, 1], upper = confint(best73, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("TK", "C.N")),
         type = case_when(variable %in% c("TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_network_complexity_decomp <- ggplot(fit_data61) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_network_complexity_decomp

## (Figure 5I) end


fit81 <- lm(beta_diversity_fungi ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit81)
summary(fit81)
alias(fit81)
vif(lm(beta_diversity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit82 <- lm(beta_diversity_fungi ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit82)
summary(fit82)
fit83 <- dredge(fit82)
subset(fit83, delta < 2)
fit84 <- lm(beta_diversity_fungi ~ SWC + TN + TK + C.N, data = diversity_multifun1)
check_normality(fit84)
summary(fit84)
AICc(fit81, fit82, fit84)
r2(fit81); r2(fit82); r2(fit84)
tab_model(fit84)
fit_data71 <- summary(fit84)$coefficients[, 1] %>% as.data.frame() %>% 
  mutate(lower = confint(fit84, full = TRUE)[, 1], upper = confint(fit84, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "TK", "TN", "C.N")),
         type = case_when(variable == "SWC" ~ "Soil structure",
                          variable %in% c("TN", "TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_beta_diversity_fungi <- ggplot(fit_data71) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_beta_diversity_fungi

## (Figure 5D) end


fit91 <- lm(beta_diversity_predator ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit91)
summary(fit91)
alias(fit91)
vif(lm(beta_diversity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit92 <- lm(beta_diversity_predator ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit92)
summary(fit92)
fit93 <- dredge(fit92)
subset(fit93, delta < 2)
best93 <- model.avg(fit93, subset = delta < 2)
summary(best93)
sw(best93)
fit94 <- lm(beta_diversity_predator ~ TK + SWC + Por + C.N + TN, data = diversity_multifun1)
check_normality(fit94)
summary(fit94)
AICc(fit91, fit92, fit94)
r2(fit91); r2(fit92); r2(fit94)
tab_model(fit94)
fit_data81 <- summary(best93)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best93, full = TRUE)[, 1], upper = confint(best93, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("SWC", "Por", "TK", "TN", "C.N")),
         type = case_when(variable %in% c("SWC", "Por") ~ "Soil structure",
                          variable %in% c("TN", "TK") ~ "Nutrient cycling",
                          variable == "C.N" ~ "Nutrient balance"))
f_beta_diversity_predator <- ggplot(fit_data81) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_beta_diversity_predator

## (Figure 5E) end


fit101 <- lm(beta_diversity_decomp ~ pH + Por + BD + SWC + NWC + SOM + TN + TP + TK + C.N + C.P + N.P, data = diversity_multifun1)
check_normality(fit101)
summary(fit101)
alias(fit101)
vif(lm(beta_diversity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1))
fit102 <- lm(beta_diversity_decomp ~ pH + Por + SWC + NWC + TN + TP + TK + C.N, data = diversity_multifun1)
check_normality(fit102)
summary(fit102)
fit103 <- dredge(fit102)
subset(fit103, delta < 2)
best103 <- model.avg(fit103, subset = delta < 2)
summary(best103)
sw(best103)
fit104 <- lm(beta_diversity_decomp ~ TN + SWC + NWC + pH + TP + Por, data = diversity_multifun1)
check_normality(fit104)
summary(fit104)
AICc(fit101, fit102, fit104)
r2(fit101); r2(fit102); r2(fit104)
tab_model(fit104)
fit_data91 <- summary(best103)$coefficients[1,] %>% as.data.frame() %>% 
  mutate(lower = confint(best103, full = TRUE)[, 1], upper = confint(best103, full = TRUE)[, 2]) %>% 
  slice(-1) %>% dplyr::rename(Estimate = ".") %>% mutate(variable = row.names(.)) %>% 
  mutate(variable = factor(variable, levels = c("Por", "NWC", "SWC", "TP", "pH", "TN")),
         type = case_when(variable %in% c("SWC", "Por", "NWC") ~ "Soil structure",
                          variable %in% c("TN", "TP", "pH") ~ "Nutrient cycling"))
f_beta_diversity_decomp <- ggplot(fit_data91) + 
  geom_hline(aes(yintercept = 0), linewidth = 0.5, colour = "grey50", linetype = 2) +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_bw_no_facet() +
  ylab("Parameter estimates"); f_beta_diversity_decomp

## (Figure 5F) end

P_MLR <- plot_grid(f_alpha_diversity_fungi, f_alpha_diversity_predator, f_alpha_diversity_decomp,
                   f_beta_diversity_fungi, f_beta_diversity_predator, f_beta_diversity_decomp,
                   f_network_complexity_fungi, f_network_complexity_predator, f_network_complexity_decomp,
                   ncol = 3, nrow = 3, align = "hv"); P_MLR
# ggsave(filename = "P_MLR .pdf", plot = P_MLR , device = "pdf", width = 15, height = 14, units = "cm", dpi = 300)






#### 6 Direct and indirect effects (Figure 6) ####
allVars0 <- qw(Lat, Alt)
envir_multifun0 <- cbind(envir, normalize(envir, allVars0))
data_psem <- data.frame(envir[,1:2], scale(envir[,3:17]), Geospatial_meanFunction = envir_multifun0$meanFunction, alpha_multifun[,-1], complexity_multifun[,-1], beta_multifun[,-1])

vif(lm(network_complexity_predator ~ pH + TK + C.N + SWC + NWC, data_psem))
compsite_soil_pre <- lm(network_complexity_predator ~ pH + TK + C.N + SWC + NWC, data_psem)
check_normality(compsite_soil_pre)
summary(compsite_soil_pre)
coefs(compsite_soil_pre, standardize = "scale")

Intercept2 <- summary(compsite_soil_pre)$coefficients[1, 1];Intercept2
pH2 <- summary(compsite_soil_pre)$coefficients[2, 1];pH2
TK2 <- summary(compsite_soil_pre)$coefficients[3, 1];TK2
C.N2 <- summary(compsite_soil_pre)$coefficients[4, 1];C.N2
SWC2 <- summary(compsite_soil_pre)$coefficients[5, 1];SWC2
NWC2 <- summary(compsite_soil_pre)$coefficients[6, 1];NWC2
compsite_soil_pre_2 <- pH2 * data_psem$pH + TK2 * data_psem$TK + C.N2 * data_psem$C.N + SWC2 * data_psem$SWC + NWC2 * data_psem$NWC + Intercept2

data_psem <- cbind(data_psem, compsite_soil_pre_2)
summary(lm(network_complexity_predator ~ compsite_soil_pre_2, data_psem))

vif(lm(alpha_diversity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2, data = data_psem))
vif(lm(beta_diversity_predator ~ Geospatial_meanFunction + alpha_diversity_predator + compsite_soil_pre_2, data = data_psem))
vif(lm(network_complexity_predator ~ Geospatial_meanFunction + alpha_diversity_predator + beta_diversity_predator + compsite_soil_pre_2, data = data_psem))

predator_psem <- psem(
  lm(compsite_soil_pre_2 ~ Geospatial_meanFunction, data = data_psem),
  lm(alpha_diversity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2, data = data_psem),
  lm(beta_diversity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2, data = data_psem),
  lm(network_complexity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2 + alpha_diversity_predator + beta_diversity_predator, data = data_psem),
  alpha_diversity_fungi %~~% beta_diversity_fungi, data = data_psem
)
summary(predator_psem)

predator_psem <- psem(
  lm(compsite_soil_pre_2 ~ Geospatial_meanFunction, data = data_psem),
  lm(alpha_diversity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2, data = data_psem),
  lm(beta_diversity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2 + alpha_diversity_predator, data = data_psem),
  lm(network_complexity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2 + alpha_diversity_predator + beta_diversity_predator, data = data_psem)
)
summary(predator_psem)

predator_psem1 <- update(predator_psem, network_complexity_predator ~ Geospatial_meanFunction + compsite_soil_pre_2 + beta_diversity_predator)
summary(predator_psem1)
AIC(predator_psem, predator_psem1)

predator_psem2 <- update(predator_psem1, network_complexity_predator ~ compsite_soil_pre_2 + beta_diversity_predator)
summary(predator_psem2)
AIC(predator_psem, predator_psem1, predator_psem2)

plot(predator_psem2)

## (Figure 6B) end


vif(lm(network_complexity_fungi ~ TK + SWC + NWC + TN + C.N + TP, data_psem))
compsite_soil_fu <- lm(network_complexity_fungi ~ TK + SWC + NWC + TN + C.N + TP, data_psem)
check_normality(compsite_soil_fu)
summary(compsite_soil_fu)
coefs(compsite_soil_fu, standardize = "scale")

Intercept1 <- summary(compsite_soil_fu)$coefficients[1, 1];Intercept1
TK1 <- summary(compsite_soil_fu)$coefficients[2, 1];TK1
SWC1 <- summary(compsite_soil_fu)$coefficients[3, 1];SWC1
NWC1 <- summary(compsite_soil_fu)$coefficients[4, 1];NWC1
TN1 <- summary(compsite_soil_fu)$coefficients[5, 1];TN1
C.N1 <- summary(compsite_soil_fu)$coefficients[6, 1];C.N1
TP1 <- summary(compsite_soil_fu)$coefficients[7, 1];TP1
compsite_soil_fu_1 <- TK1 * data_psem$TK + SWC1 * data_psem$SWC + TN1 * data_psem$TN + C.N1 * data_psem$C.N + NWC1 * data_psem$NWC + TP1 * data_psem$TP + Intercept1
data_psem <- cbind(data_psem, compsite_soil_fu_1)
names(data_psem)
summary(lm(network_complexity_fungi ~ compsite_soil_fu_1, data_psem))

vif(lm(alpha_diversity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1, data = data_psem))
vif(lm(beta_diversity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1, data = data_psem))
vif(lm(network_complexity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1 + alpha_diversity_fungi + beta_diversity_fungi, data = data_psem))

fungi_list <- list(
  lm(compsite_soil_fu_1 ~ Geospatial_meanFunction, data = data_psem),
  lm(alpha_diversity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1, data = data_psem),
  lm(beta_diversity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1, data = data_psem),
  lm(network_complexity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1 + alpha_diversity_fungi + beta_diversity_fungi, data = data_psem),
  alpha_diversity_fungi %~~% beta_diversity_fungi
)

fungi_psem <- as.psem(fungi_list)
summary(fungi_psem)

fungi_psem1 <- update(fungi_psem, network_complexity_fungi ~ Geospatial_meanFunction + compsite_soil_fu_1 + alpha_diversity_fungi)
summary(fungi_psem1)
AIC(fungi_psem, fungi_psem1)

fungi_psem2 <- update(fungi_psem1, network_complexity_fungi ~ compsite_soil_fu_1 + alpha_diversity_fungi)
summary(fungi_psem2)
AIC(fungi_psem, fungi_psem1, fungi_psem2)

plot(fungi_psem2)

## (Figure 6A) end


vif(lm(network_complexity_decomp ~ TK + C.N, data_psem))
compsite_soil_dec <- lm(network_complexity_decomp ~ TK + C.N, data_psem)
check_normality(compsite_soil_dec)
summary(compsite_soil_dec)
coefs(compsite_soil_dec, standardize = "scale")

Intercept3 <- summary(compsite_soil_dec)$coefficients[1, 1];print(Intercept3)
TK3 <- summary(compsite_soil_dec)$coefficients[2, 1];TK3
C.N3 <- summary(compsite_soil_dec)$coefficients[3, 1];C.N3
compsite_soil_dec_3 <- TK3 * data_psem$TK + C.N3 * data_psem$C.N + Intercept3

data_psem <- cbind(data_psem, compsite_soil_dec_3)
names(data_psem)
summary(lm(network_complexity_decomp ~ compsite_soil_dec_3, data_psem))

vif(lm(alpha_diversity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3, data = data_psem))
vif(lm(beta_diversity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3 + alpha_diversity_decomp, data = data_psem))
vif(lm(network_complexity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3 + beta_diversity_decomp + alpha_diversity_decomp, data = data_psem))

names(data_psem)
decomp_list <- list(
  lm(compsite_soil_dec_3 ~ Geospatial_meanFunction, data = data_psem),
  lm(alpha_diversity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3, data = data_psem),
  lm(beta_diversity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3, data = data_psem),
  lm(network_complexity_decomp ~ Geospatial_meanFunction + compsite_soil_dec_3 + beta_diversity_decomp + alpha_diversity_decomp, data = data_psem),
  alpha_diversity_decomp %~~% beta_diversity_decomp
)
decomp_psem <- as.psem(decomp_list)
summary(decomp_psem)

decomp_psem1 <- psem(
  #lm(compsite_soil_dec_3 ~ Geospatial_meanFunction, data = data_psem),
  lm(alpha_diversity_decomp ~ compsite_soil_dec_3, data = data_psem),
  lm(beta_diversity_decomp ~ compsite_soil_dec_3, data = data_psem),
  lm(network_complexity_decomp ~ compsite_soil_dec_3 + beta_diversity_decomp + alpha_diversity_decomp, data = data_psem),
  alpha_diversity_decomp %~~% beta_diversity_decomp,
  data = data_psem
)
summary(decomp_psem1)
AIC(decomp_psem, decomp_psem1)

decomp_psem2 <- update(decomp_psem1, network_complexity_decomp ~ beta_diversity_decomp + alpha_diversity_decomp)
summary(decomp_psem2)
AIC(decomp_psem, decomp_psem1, decomp_psem2)

decomp_psem3 <- update(decomp_psem2, network_complexity_decomp ~ alpha_diversity_decomp)
summary(decomp_psem3)
AIC(decomp_psem, decomp_psem1, decomp_psem2, decomp_psem3)

plot(decomp_psem3)

## (Figure 6C) end


build_graph_a <- function(coefs_df) {
  graph <- list()
  for (i in 1:nrow(coefs_df)) {
    from <- coefs_df$Predictor[i]
    to <- coefs_df$Response[i]
    effect <- coefs_df$effect[i]
    graph[[from]] <- if (!is.null(graph[[from]])) {
      rbind(graph[[from]], data.frame(to = to, effect = effect))
    } else {
      data.frame(to = to, effect = effect)
    }
  }
  return(graph)
}

find_all_paths_a <- function(graph, start, end, visited = character()) {
  if (start == end) return(list(c(end)))
  if (!start %in% names(graph)) return(list())
  visited <- c(visited, start)
  paths <- list()
  for (i in 1:nrow(graph[[start]])) {
    next_node <- graph[[start]]$to[i]
    if (!(next_node %in% visited)) {
      sub_paths <- find_all_paths_a(graph, next_node, end, visited)
      for (sp in sub_paths) {
        paths <- c(paths, list(c(start, sp)))
      }
    }
  }
  return(paths)
}

path_effect_a <- function(path, graph) {
  eff <- 1
  for (i in 1:(length(path)-1)) {
    from <- path[i]
    to <- path[i+1]
    edge <- graph[[from]]
    eff <- eff * edge$effect[edge$to == to]
  }
  return(eff)
}

compute_all_effects_a <- function(graph, variables) {
  library(dplyr)
  results <- tibble()
  for (pred in variables) {
    for (outc in variables) {
      if (pred != outc) {
        paths <- find_all_paths_a(graph, pred, outc)
        if (length(paths) > 0) {
          direct_effect <- sum(sapply(paths, function(p) if(length(p) == 2) path_effect_a(p, graph) else 0))
          indirect_effect <- sum(sapply(paths, function(p) if(length(p) > 2) path_effect_a(p, graph) else 0))
          results <- bind_rows(results, tibble(
            predictor = pred,
            outcome = outc,
            direct_effect = direct_effect,
            indirect_effect = indirect_effect,
            total_effect = direct_effect + indirect_effect
          ))
        }
      }
    }
  }
  return(results)
}

prepare_plot_data <- function(effects_df, desired_order, outcome_filter = "TRAD") {
  library(dplyr)
  library(reshape2)
  effects_df %>%
    melt(id.vars = c("predictor", "outcome"),
         measure.vars = c("direct_effect", "indirect_effect", "total_effect"),
         variable.name = "effect_type",
         value.name = "effect_size") %>%
    filter(predictor %in% desired_order, outcome == outcome_filter) %>%
    mutate(effect_type = factor(effect_type, levels = unique(effect_type)),
           predictor = factor(predictor, levels = desired_order))
}

coefs_df_fungi <- summary(fungi_psem2)$coefficients
coefs_df_fungi <- coefs_df_fungi[!grepl("~~", coefs_df_fungi$Predictor) & !grepl("~~", coefs_df_fungi$Response), ]
coefs_df_fungi$effect <- coefs_df_fungi$Std.Estimate
graph_fungi <- build_graph_a(coefs_df_fungi)
all_vars_fungi <- unique(c(coefs_df_fungi$Response, coefs_df_fungi$Predictor))
effects_df_fungi <- compute_all_effects_a(graph_fungi, all_vars_fungi)
effects_df_plot_fungi <- prepare_plot_data(
  effects_df_fungi,
  desired_order = c("Geospatial_meanFunction", "compsite_soil_fu_1", "alpha_diversity_fungi", "beta_diversity_fungi"),
  outcome_filter = "network_complexity_fungi"
)

effects_fungi_p <- ggplot(effects_df_plot_fungi, aes(x = effect_size, y = reorder(predictor, effect_size), fill = effect_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.5, color = "gray50") +
  geom_text(aes(label = round(effect_size, 2), group = effect_type), position = position_dodge(0.8), hjust = 0.4, vjust = 0.5, size = 3, color = "black") +
  scale_fill_manual(values = c("direct_effect" = "#95e1d3", "indirect_effect" = "#fce38a", "total_effect" = "#E41A1C")) +
  theme_bw_no_facet() +
  theme(axis.text.y = element_blank()); effects_fungi_p

## (Figure 6D) end


coefs_df_decomp <- summary(decomp_psem3)$coefficients
coefs_df_decomp <- coefs_df_decomp[!grepl("~~", coefs_df_decomp$Predictor) & !grepl("~~", coefs_df_decomp$Response), ]
coefs_df_decomp$effect <- coefs_df_decomp$Std.Estimate
graph_decomp <- build_graph_a(coefs_df_decomp)
all_vars_decomp <- unique(c(coefs_df_decomp$Response, coefs_df_decomp$Predictor))
effects_df_decomp <- compute_all_effects_a(graph_decomp, all_vars_decomp)
effects_df_plot_decomp <- prepare_plot_data(
  effects_df_decomp,
  desired_order = c("compsite_soil_dec_3", "alpha_diversity_decomp", "beta_diversity_decomp"),
  outcome_filter = "network_complexity_decomp"
)

effects_decomp_p <- ggplot(effects_df_plot_decomp, aes(x = effect_size, y = reorder(predictor, effect_size), fill = effect_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.5, color = "gray50") +
  geom_text(aes(label = round(effect_size, 2), group = effect_type), position = position_dodge(0.8), hjust = 0.4, vjust = 0.5, size = 3, color = "black") +
  scale_fill_manual(values = c("direct_effect" = "#95e1d3", "indirect_effect" = "#fce38a", "total_effect" = "#E41A1C")) +
  theme_bw_no_facet() +
  theme(axis.text.y = element_blank()); effects_decomp_p

## (Figure 6E) end


coefs_df_predator <- summary(predator_psem2)$coefficients
coefs_df_predator <- coefs_df_predator[!grepl("~~", coefs_df_predator$Predictor) & !grepl("~~", coefs_df_predator$Response), ]
coefs_df_predator$effect <- coefs_df_predator$Std.Estimate
graph_predator <- build_graph_a(coefs_df_predator)
all_vars_predator <- unique(c(coefs_df_predator$Response, coefs_df_predator$Predictor))
effects_df_predator <- compute_all_effects_a(graph_predator, all_vars_predator)
effects_df_plot_predator <- prepare_plot_data(
  effects_df_predator,
  desired_order = c("Geospatial_meanFunction", "compsite_soil_pre_2", "beta_diversity_predator"),
  outcome_filter = "network_complexity_predator"
)

effects_predator_p <- ggplot(effects_df_plot_predator, aes(x = effect_size, y = reorder(predictor, effect_size), fill = effect_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.5, color = "gray50") +
  geom_text(aes(label = round(effect_size, 2), group = effect_type), position = position_dodge(0.8), hjust = 0.4, vjust = 0.5, size = 3, color = "black") +
  scale_fill_manual(values = c("direct_effect" = "#95e1d3", "indirect_effect" = "#fce38a", "total_effect" = "#E41A1C")) +
  theme_bw_no_facet() +
  theme(axis.text.y = element_blank()); effects_predator_p

## (Figure 6F) end


P_sem_effe <- plot_grid(effects_fungi_p, effects_predator_p, effects_decomp_p, ncol = 3, nrow = 1, align = "hv"); P_sem_effe

# ggsave(filename = "P_sem_effe.pdf", plot = P_sem_effe, device = "pdf", width = 14, height = 4, units = "cm", dpi = 300)
## end

