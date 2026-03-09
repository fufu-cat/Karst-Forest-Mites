## These networks had already been constructed in the previous code and represent mite co-occurrence networks across different trophic guilds and forest habitat types.
network_list <- list(
  Total      = g_m_total,
  Fungi      = g_m_fungi,
  Predator   = g_m_predator,
  Decomposer = g_m_decomp,
  RBS        = g_m_RBS,
  DBF        = g_m_DBF,
  EBF        = g_m_EBF,
  CBF        = g_m_CBF
  )



## 1 R code used for Zi–Pi analysis, keystone-node removal simulation, and node vulnerability analysis.

analyze_mite_network <- function(g,
                                 net_name = "network",
                                 nrep = 1000,
                                 seed = 123) {
  
  # ----------------------------
  # 0. helper functions
  # ----------------------------
  largest_component_size <- function(graph) {
    if (igraph::vcount(graph) == 0) return(0)
    max(igraph::components(graph)$csize)
  }
  
  clean_taxon <- function(x) {
    stringr::str_replace(x, "^X\\.[0-9]+\\.\\.", "")
  }
  
  # ----------------------------
  # 1. Zi-Pi analysis
  # ----------------------------
  res <- ggClusterNet::ZiPiPlot(igraph = g, method = "cluster_fast_greedy")
  
  taxa.roles <- res[[2]] %>%
    tibble::rownames_to_column(var = "taxon") %>%
    dplyr::mutate(
      genus = clean_taxon(taxon),
      label_plot = ifelse(is.na(label) | label == "", genus, clean_taxon(label))
    ) %>%
    dplyr::filter(is.finite(z), is.finite(p))
  
  keystone <- taxa.roles %>%
    dplyr::filter(z > 2.5 | p > 0.62)
  
  keystone_nodes <- keystone$taxon
  
  # ----------------------------
  # 2. Keystone removal vs random removal
  # ----------------------------
  if (length(keystone_nodes) > 0) {
    
    g_keystone_removed <- igraph::delete_vertices(g, keystone_nodes)
    
    keystone_lcc <- largest_component_size(g_keystone_removed) / igraph::vcount(g)
    keystone_ncomp <- igraph::components(g_keystone_removed)$no
    
    set.seed(seed)
    random_lcc <- numeric(nrep)
    random_ncomp <- numeric(nrep)
    
    for (i in seq_len(nrep)) {
      random_nodes <- sample(igraph::V(g)$name, length(keystone_nodes))
      g_random_removed <- igraph::delete_vertices(g, random_nodes)
      
      random_lcc[i] <- largest_component_size(g_random_removed) / igraph::vcount(g)
      random_ncomp[i] <- igraph::components(g_random_removed)$no
    }
    
    empirical_prop_lcc <- mean(random_lcc <= keystone_lcc)
    empirical_prop_ncomp <- mean(random_ncomp >= keystone_ncomp)
    
  } else {
    
    keystone_lcc <- NA_real_
    keystone_ncomp <- NA_real_
    random_lcc <- numeric(0)
    random_ncomp <- numeric(0)
    empirical_prop_lcc <- NA_real_
    empirical_prop_ncomp <- NA_real_
  }
  
  # ----------------------------
  # 3. Node vulnerability
  # ----------------------------
  original_lcc <- largest_component_size(g)
  
  vul_table <- tibble::tibble(
    node = igraph::V(g)$name,
    vulnerability = sapply(igraph::V(g)$name, function(v) {
      g_removed <- igraph::delete_vertices(g, v)
      removed_lcc <- largest_component_size(g_removed)
      (original_lcc - removed_lcc) / original_lcc
    })
  ) %>%
    dplyr::arrange(dplyr::desc(vulnerability)) %>%
    dplyr::mutate(genus = clean_taxon(node))
  
  keystone_vulnerability <- vul_table %>%
    dplyr::filter(node %in% keystone_nodes)
  
  # ----------------------------
  # 4. Zi-Pi plot (your final selected version)
  # ----------------------------
  rect_df <- data.frame(
    xmin = c(0, 0.62, 0, 0.62),
    xmax = c(0.62, 1, 0.62, 1),
    ymin = c(-Inf, 2.5, 2.5, -Inf),
    ymax = c(2.5, Inf, Inf, 2.5),
    lab  = c("Peripherals", "Network hubs", "Module hubs", "Connectors")
  )
  
  p_zipi <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rect_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lab),
      inherit.aes = FALSE,
      alpha = 0.18,
      color = NA
    ) +
    ggplot2::geom_vline(
      xintercept = 0.62,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = 2.5,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = taxa.roles,
      ggplot2::aes(x = p, y = z),
      color = "grey45",
      size = 2,
      alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = keystone,
      ggplot2::aes(x = p, y = z),
      shape = 21,
      fill = "#D95F02",
      color = "black",
      size = 4,
      stroke = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = keystone,
      ggplot2::aes(x = p, y = z, label = label_plot),
      size = 3.4,
      color = "black",
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey40",
      segment.size = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Participation coefficient (Pi)",
      y = "Within-module connectivity (Zi)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::guides(fill = "none", color = "none")
  
  # ----------------------------
  # 5. Robustness plots (bar plot, hollow bars, black outline)
  # ----------------------------
  if (length(keystone_nodes) > 0) {
    
    robust_df <- tibble::tibble(
      strategy = factor(c("Keystone removal", "Random removal"),
                        levels = c("Keystone removal", "Random removal")),
      lcc = c(keystone_lcc, mean(random_lcc)),
      lcc_sd = c(NA_real_, stats::sd(random_lcc)),
      ncomp = c(keystone_ncomp, mean(random_ncomp)),
      ncomp_sd = c(NA_real_, stats::sd(random_ncomp))
    )
    
    robust_df <- tibble::tibble(
      
      strategy = factor(
        c("Keystone removal", "Random removal"),
        levels = c("Keystone removal", "Random removal")
      ),
      
      lcc = c(keystone_lcc, mean(random_lcc)),
      
      lcc_sd = c(NA, sd(random_lcc)),
      
      ncomp = c(keystone_ncomp, mean(random_ncomp)),
      
      ncomp_sd = c(NA, sd(random_ncomp))
      
    )
    
    cols <- c(
      "Keystone removal" = "#D95F02",
      "Random removal" = "grey65"
    )
    
    p_lcc <- ggplot2::ggplot(
      robust_df,
      ggplot2::aes(strategy, lcc, fill = strategy)
    ) +
      
      ggplot2::geom_col(
        width = 0.6,
        color = "black",
        linewidth = 0.5
      ) +
      
      ggplot2::geom_errorbar(
        data = robust_df %>% dplyr::filter(strategy == "Random removal"),
        ggplot2::aes(ymin = lcc - lcc_sd, ymax = lcc + lcc_sd),
        width = 0.15,
        linewidth = 0.6
      ) +
      
      ggplot2::scale_fill_manual(values = cols) +
      
      ggplot2::theme_bw(base_size = 12) +
      
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      ) +
      
      ggplot2::labs(
        x = NULL,
        y = "Relative size of largest connected component"
      )
    
    
    p_ncomp <- ggplot2::ggplot(
      robust_df,
      ggplot2::aes(strategy, ncomp, fill = strategy)
    ) +
      
      ggplot2::geom_col(
        width = 0.6,
        color = "black",
        linewidth = 0.5
      ) +
      
      ggplot2::geom_errorbar(
        data = robust_df %>% dplyr::filter(strategy == "Random removal"),
        ggplot2::aes(ymin = ncomp - ncomp_sd, ymax = ncomp + ncomp_sd),
        width = 0.15,
        linewidth = 0.6
      ) +
      
      ggplot2::scale_fill_manual(values = cols) +
      
      ggplot2::theme_bw(base_size = 12) +
      
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      ) +
      
      ggplot2::labs(
        x = NULL,
        y = "Number of connected components"
      )
    
  } else {
    
    p_lcc <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "No key nodes detected", size = 5) +
      ggplot2::theme_void()
    
    p_ncomp <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "No key nodes detected", size = 5) +
      ggplot2::theme_void()
  }
  
  # ----------------------------
  # 6. Vulnerability plot (clean style)
  # ----------------------------
  top_vul <- vul_table %>%
    dplyr::slice_max(vulnerability, n = 10) %>%
    dplyr::mutate(
      key_flag = ifelse(node %in% keystone_nodes, "Key node", "Other")
    )
  
  p_vul <- ggplot2::ggplot(
    top_vul,
    ggplot2::aes(
      x = stats::reorder(genus, vulnerability),
      y = vulnerability,
      fill = key_flag
    )
  ) +
    ggplot2::geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.7
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Key node" = "#D95F02", "Other" = "grey65")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Node vulnerability"
    )
  
  # ----------------------------
  # 7. Summary table
  # ----------------------------
  summary_table <- tibble::tibble(
    network = net_name,
    n_nodes = igraph::vcount(g),
    n_edges = igraph::ecount(g),
    n_connectors = sum(taxa.roles$roles == "Connectors", na.rm = TRUE),
    n_module_hubs = sum(taxa.roles$roles == "Module hubs", na.rm = TRUE),
    n_network_hubs = sum(taxa.roles$roles == "Network hubs", na.rm = TRUE),
    n_peripherals = sum(taxa.roles$roles == "Peripherals", na.rm = TRUE),
    n_keystone = nrow(keystone),
    keystone_lcc = keystone_lcc,
    random_lcc_mean = ifelse(length(random_lcc) > 0, mean(random_lcc), NA_real_),
    random_lcc_sd = ifelse(length(random_lcc) > 0, stats::sd(random_lcc), NA_real_),
    empirical_prop_lcc = empirical_prop_lcc,
    keystone_ncomp = keystone_ncomp,
    random_ncomp_mean = ifelse(length(random_ncomp) > 0, mean(random_ncomp), NA_real_),
    random_ncomp_sd = ifelse(length(random_ncomp) > 0, stats::sd(random_ncomp), NA_real_),
    empirical_prop_ncomp = empirical_prop_ncomp
  )
  
  # ----------------------------
  # 8. Return only
  # ----------------------------
  return(list(
    summary = summary_table,
    taxa_roles = taxa.roles,
    keystone = keystone,
    vulnerability = vul_table,
    keystone_vulnerability = keystone_vulnerability,
    random_lcc = random_lcc,
    random_ncomp = random_ncomp,
    zipi_plot = p_zipi,
    robustness_lcc_plot = p_lcc,
    robustness_ncomp_plot = p_ncomp,
    vulnerability_plot = p_vul
  ))
}



## Example
res_CBF <- analyze_mite_network(
  g = g_m_CBF,
  net_name = "CBF",
  nrep = 1000
)


print(res_CBF$summary, width = Inf)
res_CBF$keystone
res_CBF$keystone_vulnerability

res_CBF$zipi_plot
res_CBF$robustness_lcc_plot
res_CBF$robustness_ncomp_plot
res_CBF$vulnerability_plot





## The R workflow used for Zi–Pi analysis, keystone-node removal simulation, and batch processing of the eight networks is provided below.

run_all_mite_networks <- function(network_list, nrep = 1000) {
  
  all_res <- lapply(names(network_list), function(nm) {
    analyze_mite_network(
      g = network_list[[nm]],
      net_name = nm,
      nrep = nrep
    )
  })
  names(all_res) <- names(network_list)
  
  summary_all <- dplyr::bind_rows(lapply(all_res, function(x) x$summary))
  
  keystone_all <- dplyr::bind_rows(
    lapply(names(all_res), function(nm) {
      all_res[[nm]]$keystone %>%
        dplyr::mutate(network = nm, .before = 1)
    })
  )
  
  keystone_vulnerability_all <- dplyr::bind_rows(
    lapply(names(all_res), function(nm) {
      all_res[[nm]]$keystone_vulnerability %>%
        dplyr::mutate(network = nm, .before = 1)
    })
  )
  
  vulnerability_all <- dplyr::bind_rows(
    lapply(names(all_res), function(nm) {
      all_res[[nm]]$vulnerability %>%
        dplyr::mutate(network = nm, .before = 1)
    })
  )
  
   zipi_plots <- Map(function(p, nm) {
    p + ggplot2::ggtitle(nm) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
      )
  }, lapply(all_res, function(x) x$zipi_plot), names(all_res))
  
  lcc_plots <- Map(function(p, nm) {
    p + ggplot2::ggtitle(nm) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
      )
  }, lapply(all_res, function(x) x$robustness_lcc_plot), names(all_res))
  
  ncomp_plots <- Map(function(p, nm) {
    p + ggplot2::ggtitle(nm) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
      )
  }, lapply(all_res, function(x) x$robustness_ncomp_plot), names(all_res))
  
  vul_plots <- Map(function(p, nm) {
    p + ggplot2::ggtitle(nm) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
      )
  }, lapply(all_res, function(x) x$vulnerability_plot), names(all_res))
  

  zipi_panel <- patchwork::wrap_plots(zipi_plots, ncol = 4)
  lcc_panel <- patchwork::wrap_plots(lcc_plots, ncol = 4)
  ncomp_panel <- patchwork::wrap_plots(ncomp_plots, ncol = 4)
  vul_panel <- patchwork::wrap_plots(vul_plots, ncol = 4)
  
  return(list(
    all_res = all_res,
    summary_all = summary_all,
    keystone_all = keystone_all,
    keystone_vulnerability_all = keystone_vulnerability_all,
    vulnerability_all = vulnerability_all,
    zipi_panel = zipi_panel,
    lcc_panel = lcc_panel,
    ncomp_panel = ncomp_panel,
    vul_panel = vul_panel
  ))
}




all_output <- run_all_mite_networks(network_list, nrep = 1000)

print(all_output$summary_all, width = Inf)
as.data.frame(all_output$keystone_all)
print(all_output$keystone_vulnerability_all, width = Inf)
print(all_output$vulnerability_all, width = Inf)

all_output$zipi_panel
all_output$lcc_panel
all_output$ncomp_panel
all_output$vul_panel



# ggsave("zipi_panel.pdf", all_output$zipi_panel, width = 16, height = 8)
# ggsave("lcc_panel.pdf", all_output$lcc_panel, width = 16, height = 8)
# ggsave("ncomp_panel.pdf", all_output$ncomp_panel, width = 16, height = 8)
# ggsave("vulnerability_panel.pdf", all_output$vul_panel, width = 16, height = 8)

# write.csv(all_output$summary_all, "mite_network_summary.csv", row.names = FALSE)
# write.csv(as.data.frame(all_output$keystone_all), "mite_keystone_species.csv", row.names = FALSE)
# write.csv(all_output$keystone_vulnerability_all, "mite_keystone_vulnerability.csv", row.names = FALSE)
# write.csv(all_output$vulnerability_all, "mite_all_node_vulnerability.csv", row.names = FALSE)

