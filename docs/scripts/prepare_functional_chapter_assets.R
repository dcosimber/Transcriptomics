#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

project_dir <- normalizePath(getwd(), mustWork = TRUE)
default_source_dir <- "/mnt/lustre/scratch/nlsas/home/otras/pia/dci/holofish_rnaseq/analysis/results/05_functional_enrichment"
source_dir <- Sys.getenv("FUNCTIONAL_SOURCE_DIR", unset = default_source_dir)
if (!dir.exists(source_dir)) {
  stop("Functional enrichment source directory not found: ", source_dir, call. = FALSE)
}

contrast_order <- c(
  "N_OD_vs_N_BL",
  "P_OD_vs_P_OL",
  "P_BD_vs_P_BL",
  "P_OL_vs_P_BL",
  "P_BD_vs_P_OD",
  "P_surface_x_pigmentation",
  "P_OL_vs_N_OD",
  "P_BD_vs_N_BL"
)

short_label <- c(
  N_OD_vs_N_BL = "N_OD vs N_BL",
  P_OD_vs_P_OL = "P_OD vs P_OL",
  P_BD_vs_P_BL = "P_BD vs P_BL",
  P_OL_vs_P_BL = "P_OL vs P_BL",
  P_BD_vs_P_OD = "P_BD vs P_OD",
  P_surface_x_pigmentation = "Surface x pigment",
  P_OL_vs_N_OD = "P_OL vs N_OD",
  P_BD_vs_N_BL = "P_BD vs N_BL"
)

contrast_question <- c(
  N_OD_vs_N_BL = "Referencia fisiologica de pigmentacion y superficie.",
  P_OD_vs_P_OL = "Cambios funcionales dentro de la superficie ocular pseudoalbina.",
  P_BD_vs_P_BL = "Pigmentacion ectopica dentro de la superficie blind pseudoalbina.",
  P_OL_vs_P_BL = "Efecto de superficie bajo color claro pseudoalbino.",
  P_BD_vs_P_OD = "Efecto de superficie bajo color oscuro pseudoalbino.",
  P_surface_x_pigmentation = "Interaccion superficie por pigmentacion.",
  P_OL_vs_N_OD = "Ocular claro anomalo frente al ocular oscuro normal.",
  P_BD_vs_N_BL = "Blind oscuro anomalo frente al blind claro normal."
)

tables_dir <- file.path(project_dir, "tables", "06_functional")
figures_dir <- file.path(project_dir, "figures", "06_functional")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(tables_dir, "by_contrast"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "summary"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "contrasts"), showWarnings = FALSE, recursive = TRUE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

as_num <- function(x) suppressWarnings(as.numeric(x))

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", format(round(as_num(x), digits), nsmall = digits, trim = TRUE))
}

fmt_sci <- function(x) {
  ifelse(is.na(x), "", formatC(as_num(x), format = "e", digits = 2))
}

wrap_text <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

strip_go_prefix <- function(x) {
  x <- gsub("^[A-Z]{2}: ", "", x)
  x
}

category_patterns <- list(
  pigmentacion_melanocito = "melano|pigment|tyrosinase|dopamine biosynthetic|tetrahydrobiopterin",
  matriz_adhesion = "extracellular matrix|collagen|integrin|adhesion|focal adhesion|basement membrane|heparin binding",
  epidermis_epitelio = "epiderm|epithel|keratin|desmosome|cell-cell junction|intermediate filament",
  musculo_sarcomero = "muscle|sarcomere|myofibril|myofilament|myosin|actin filament|z disc|m band|contract|troponin|tropomyosin",
  neural_desarrollo = "neural crest|neuron|axon|retina|ganglion|brain|nerve",
  wnt_retinoides_senalizacion = "wnt|retinoic|retinoid|nodal|activin|tgf|bmp|serotonin|signaling receptor",
  lipidos = "lipid|fatty|acyl|triacylglycerol|lipid droplet|lipid storage|ketone",
  mitocondria_metabolismo = "mitochond|respiratory chain|oxidative phosphorylation|atp synthesis|nadh|oxidoreductase",
  inmunidad_inflamacion = "immune|immun|inflamm|interleukin|chemokine|antigen|lysosome|leukocyte|myeloid|nf-kappa",
  traduccion_ribosoma = "ribosome|translation|ribosomal|rrna",
  desarrollo_morfogenesis = "development|morphogenesis|differentiation|pattern formation|regeneration|organism development"
)

category_label <- c(
  pigmentacion_melanocito = "Pigmentacion/melanocito",
  matriz_adhesion = "Matriz extracelular/adhesion",
  epidermis_epitelio = "Epidermis/epitelio",
  musculo_sarcomero = "Musculo/sarcomero",
  neural_desarrollo = "Neural/axon/cresta",
  wnt_retinoides_senalizacion = "Wnt/retinoides/senalizacion",
  lipidos = "Lipidos",
  mitocondria_metabolismo = "Mitocondria/metabolismo",
  inmunidad_inflamacion = "Inmunidad/inflamacion",
  traduccion_ribosoma = "Traduccion/ribosoma",
  desarrollo_morfogenesis = "Desarrollo/morfogenesis",
  other = "Otros"
)

assign_primary_category <- function(term) {
  term <- tolower(term)
  for (nm in names(category_patterns)) {
    if (grepl(category_patterns[[nm]], term, perl = TRUE)) return(nm)
  }
  "other"
}

expand_categories <- function(df) {
  out <- list()
  term <- tolower(paste(df$term_label, df$term_name, df$Description))
  idx <- 1
  for (cat in names(category_patterns)) {
    keep <- grepl(category_patterns[[cat]], term, perl = TRUE)
    if (any(keep)) {
      tmp <- df[keep, , drop = FALSE]
      tmp$category <- cat
      tmp$category_label <- unname(category_label[cat])
      out[[idx]] <- tmp
      idx <- idx + 1
    }
  }
  if (!length(out)) return(df[0, , drop = FALSE])
  do.call(rbind, out)
}

select_cols <- function(df, cols) {
  if (is.null(df)) {
    out <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
    colnames(out) <- cols
    return(out)
  }
  if (!is.data.frame(df)) df <- as.data.frame(df)
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) df[missing] <- NA
  df[, cols, drop = FALSE]
}

top_rows <- function(df, n = 8, increasing_effect = FALSE) {
  if (!nrow(df)) return(df)
  ord_effect <- if (increasing_effect) df$effect else -abs(df$effect)
  df[order(df$p_adjust, ord_effect), , drop = FALSE][seq_len(min(n, nrow(df))), , drop = FALSE]
}

top_one <- function(df, positive = TRUE) {
  if (!nrow(df)) return(rep(NA_character_, 4))
  if (positive) {
    df <- df[df$effect > 0, , drop = FALSE]
    df <- df[order(df$p_adjust, -df$effect), , drop = FALSE]
  } else {
    df <- df[df$effect < 0, , drop = FALSE]
    df <- df[order(df$p_adjust, df$effect), , drop = FALSE]
  }
  if (!nrow(df)) return(rep(NA_character_, 4))
  c(df$term_label[1], fmt_num(df$effect[1], 2), fmt_sci(df$p_adjust[1]), df$category_label[1])
}

gsea <- read_tsv(file.path(source_dir, "go_gsea/summary/tables/go_gsea_results_only.tsv"))
ora <- read_tsv(file.path(source_dir, "go_ora/summary/tables/go_ora_results_only.tsv"))
focused <- read_tsv(file.path(source_dir, "development_neural_crest/tables/development_neural_crest_go_terms.tsv"))
status <- read_tsv(file.path(source_dir, "summary/tables/functional_enrichment_module_status.tsv"))

gsea$p_adjust <- as_num(gsea[["p.adjust"]])
gsea$effect <- as_num(gsea$NES)
gsea$count <- as_num(gsea$setSize)
gsea$method <- "GSEA"
gsea$effect_direction <- ifelse(gsea$effect >= 0, "NES_positive", "NES_negative")
gsea$term_label <- ifelse(is.na(gsea$term_label) | gsea$term_label == "", gsea$Description, gsea$term_label)
gsea$term_clean <- strip_go_prefix(gsea$term_label)
gsea$primary_category <- vapply(paste(gsea$term_label, gsea$term_name, gsea$Description), assign_primary_category, character(1))
gsea$category_label <- unname(category_label[gsea$primary_category])

ora$p_adjust <- as_num(ora[["p.adjust"]])
ora$effect <- as_num(ora$FoldEnrichment)
ora$signed_effect <- ifelse(
  ora$direction == "down",
  -ora$effect,
  ifelse(ora$direction == "up", ora$effect, as_num(ora$zScore))
)
ora$count <- as_num(ora$Count)
ora$method <- "ORA"
ora$effect_direction <- ifelse(ora$direction == "all", "all_DEGs", ora$direction)
ora$term_label <- ifelse(is.na(ora$term_label) | ora$term_label == "", ora$Description, ora$term_label)
ora$term_clean <- strip_go_prefix(ora$term_label)
ora$primary_category <- vapply(paste(ora$term_label, ora$term_name, ora$Description), assign_primary_category, character(1))
ora$category_label <- unname(category_label[ora$primary_category])

gsea_sig <- gsea[gsea$status == "ok" & !is.na(gsea$p_adjust) & gsea$p_adjust < 0.05, , drop = FALSE]
ora_sig <- ora[ora$status == "ok" & !is.na(ora$p_adjust) & ora$p_adjust < 0.05, , drop = FALSE]
gsea_sig$contrast <- factor(gsea_sig$contrast, levels = contrast_order)
ora_sig$contrast <- factor(ora_sig$contrast, levels = contrast_order)

sig_cols <- c(
  "method", "contrast", "effect_direction", "ontology", "GO_ID", "term_label",
  "term_clean", "pvalue", "p_adjust", "qvalue", "count", "effect",
  "signed_effect", "primary_category", "category_label"
)

gsea_export <- select_cols(gsea_sig, sig_cols)
ora_export <- select_cols(ora_sig, sig_cols)
gsea_export$contrast <- as.character(gsea_export$contrast)
ora_export$contrast <- as.character(ora_export$contrast)
write_tsv(gsea_export, file.path(tables_dir, "functional_gsea_significant.tsv"))
write_tsv(ora_export, file.path(tables_dir, "functional_ora_significant.tsv"))
write_tsv(rbind(gsea_export, ora_export), file.path(tables_dir, "functional_go_significant_combined.tsv"))
write_tsv(status, file.path(tables_dir, "functional_enrichment_module_status.tsv"))

gsea_cat <- expand_categories(gsea_sig)
gsea_cat$contrast <- factor(gsea_cat$contrast, levels = contrast_order)

cat_rows <- list()
idx <- 1
for (contrast in contrast_order) {
  for (cat in names(category_label)[names(category_label) != "other"]) {
    cx <- gsea_cat[as.character(gsea_cat$contrast) == contrast & gsea_cat$category == cat, , drop = FALSE]
    if (!nrow(cx)) next
    pos <- cx[cx$effect > 0, , drop = FALSE]
    neg <- cx[cx$effect < 0, , drop = FALSE]
    dom <- cx[order(-abs(cx$effect), cx$p_adjust), , drop = FALSE][1, , drop = FALSE]
    top_pos <- top_one(cx, positive = TRUE)
    top_neg <- top_one(cx, positive = FALSE)
    cat_rows[[idx]] <- data.frame(
      contrast = contrast,
      short_label = unname(short_label[contrast]),
      category = cat,
      category_label = unname(category_label[cat]),
      n_terms = nrow(unique(cx[, c("contrast", "GO_ID", "effect_direction")])),
      n_positive = nrow(unique(pos[, c("contrast", "GO_ID", "effect_direction")])),
      n_negative = nrow(unique(neg[, c("contrast", "GO_ID", "effect_direction")])),
      dominant_effect = dom$effect[1],
      dominant_term = dom$term_label[1],
      dominant_padj = dom$p_adjust[1],
      top_positive_term = top_pos[1],
      top_positive_effect = top_pos[2],
      top_positive_padj = top_pos[3],
      top_negative_term = top_neg[1],
      top_negative_effect = top_neg[2],
      top_negative_padj = top_neg[3],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}
category_summary <- if (length(cat_rows)) do.call(rbind, cat_rows) else data.frame()
category_summary$dominant_effect_fmt <- fmt_num(category_summary$dominant_effect, 2)
category_summary$dominant_padj_fmt <- fmt_sci(category_summary$dominant_padj)
write_tsv(category_summary, file.path(tables_dir, "functional_category_summary_gsea.tsv"))

contrast_rows <- list()
idx <- 1
for (contrast in contrast_order) {
  gx <- gsea_sig[as.character(gsea_sig$contrast) == contrast, , drop = FALSE]
  ox <- ora_sig[as.character(ora_sig$contrast) == contrast, , drop = FALSE]
  gp <- top_one(gx, positive = TRUE)
  gn <- top_one(gx, positive = FALSE)
  op <- ox[ox$effect_direction == "up", , drop = FALSE]
  od <- ox[ox$effect_direction == "down", , drop = FALSE]
  oa <- ox[ox$effect_direction == "all_DEGs", , drop = FALSE]
  op <- if (nrow(op)) op[order(op$p_adjust, -op$effect), , drop = FALSE][1, ] else NULL
  od <- if (nrow(od)) od[order(od$p_adjust, -od$effect), , drop = FALSE][1, ] else NULL
  oa <- if (nrow(oa)) oa[order(oa$p_adjust, -oa$effect), , drop = FALSE][1, ] else NULL
  contrast_rows[[idx]] <- data.frame(
    contrast = contrast,
    short_label = unname(short_label[contrast]),
    question = unname(contrast_question[contrast]),
    n_gsea_fdr_0_05 = nrow(gx),
    n_gsea_positive = sum(gx$effect > 0, na.rm = TRUE),
    n_gsea_negative = sum(gx$effect < 0, na.rm = TRUE),
    n_ora_fdr_0_05 = nrow(ox),
    n_ora_all = sum(ox$effect_direction == "all_DEGs", na.rm = TRUE),
    n_ora_up = sum(ox$effect_direction == "up", na.rm = TRUE),
    n_ora_down = sum(ox$effect_direction == "down", na.rm = TRUE),
    top_gsea_positive = gp[1],
    top_gsea_positive_effect = gp[2],
    top_gsea_positive_padj = gp[3],
    top_gsea_negative = gn[1],
    top_gsea_negative_effect = gn[2],
    top_gsea_negative_padj = gn[3],
    top_ora_up = if (is.null(op)) NA else op$term_label,
    top_ora_up_fold = if (is.null(op)) NA else fmt_num(op$effect, 2),
    top_ora_up_padj = if (is.null(op)) NA else fmt_sci(op$p_adjust),
    top_ora_down = if (is.null(od)) NA else od$term_label,
    top_ora_down_fold = if (is.null(od)) NA else fmt_num(od$effect, 2),
    top_ora_down_padj = if (is.null(od)) NA else fmt_sci(od$p_adjust),
    top_ora_all = if (is.null(oa)) NA else oa$term_label,
    top_ora_all_fold = if (is.null(oa)) NA else fmt_num(oa$effect, 2),
    top_ora_all_padj = if (is.null(oa)) NA else fmt_sci(oa$p_adjust),
    stringsAsFactors = FALSE
  )
  idx <- idx + 1
}
contrast_summary <- do.call(rbind, contrast_rows)
write_tsv(contrast_summary, file.path(tables_dir, "functional_contrast_summary.tsv"))

top_all <- list()
idx <- 1
for (contrast in contrast_order) {
  gx <- gsea_sig[as.character(gsea_sig$contrast) == contrast, , drop = FALSE]
  ox <- ora_sig[as.character(ora_sig$contrast) == contrast, , drop = FALSE]
  g_top <- rbind(
    top_rows(gx[gx$effect > 0, , drop = FALSE], n = 8),
    top_rows(gx[gx$effect < 0, , drop = FALSE], n = 8)
  )
  o_top <- do.call(rbind, lapply(c("all_DEGs", "up", "down"), function(direction) {
    yy <- ox[ox$effect_direction == direction, , drop = FALSE]
    if (!nrow(yy)) return(NULL)
    yy <- yy[order(yy$p_adjust, -yy$effect), , drop = FALSE]
    head(yy, 5)
  }))
  g_out <- select_cols(g_top, sig_cols)
  o_out <- select_cols(o_top, sig_cols)
  out <- rbind(g_out, o_out)
  out$contrast <- as.character(out$contrast)
  out$effect_fmt <- fmt_num(ifelse(out$method == "ORA", out$signed_effect, out$effect), 2)
  out$padj_fmt <- fmt_sci(out$p_adjust)
  write_tsv(out, file.path(tables_dir, "by_contrast", paste0(contrast, "_functional_top_terms.tsv")))
  write_tsv(select_cols(gx, sig_cols), file.path(tables_dir, "by_contrast", paste0(contrast, "_gsea_significant.tsv")))
  write_tsv(select_cols(ox, sig_cols), file.path(tables_dir, "by_contrast", paste0(contrast, "_ora_significant.tsv")))
  top_all[[idx]] <- out
  idx <- idx + 1
}
write_tsv(do.call(rbind, top_all), file.path(tables_dir, "functional_top_terms_by_contrast.tsv"))

focused_sig <- focused[focused$status == "ok" & as_num(focused[["p.adjust"]]) < 0.10, , drop = FALSE]
focused_sig$p_adjust <- as_num(focused_sig[["p.adjust"]])
focused_sig$effect <- ifelse(focused_sig$method == "GSEA", as_num(focused_sig$NES), as_num(focused_sig$zScore))
focused_sig$effect_fmt <- fmt_num(focused_sig$effect, 2)
focused_sig$padj_fmt <- fmt_sci(focused_sig$p_adjust)
focused_export <- select_cols(focused_sig, c(
  "contrast", "category", "method", "direction", "ontology", "GO_ID",
  "term_label", "p_adjust", "effect", "effect_fmt", "padj_fmt"
))
write_tsv(focused_export, file.path(tables_dir, "functional_focused_development_terms_fdr_0_10.tsv"))

theme_functional <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title.position = "plot",
      legend.position = "right"
    )
}

counts <- rbind(
  data.frame(
    contrast = as.character(gsea_sig$contrast),
    method = "GSEA",
    direction = ifelse(gsea_sig$effect > 0, "NES positivo", "NES negativo"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    contrast = as.character(ora_sig$contrast),
    method = "ORA",
    direction = ifelse(ora_sig$effect_direction == "all_DEGs", "todos los DEGs", ora_sig$effect_direction),
    stringsAsFactors = FALSE
  )
)
counts$n <- 1
counts <- aggregate(n ~ contrast + method + direction, counts, sum)
counts$short_label <- factor(unname(short_label[counts$contrast]), levels = unname(short_label[contrast_order]))

p_counts <- ggplot(counts, aes(x = short_label, y = n, fill = direction)) +
  geom_col(width = 0.75, color = "grey25", linewidth = 0.15) +
  facet_wrap(~method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "NES positivo" = "#B6473B",
    "NES negativo" = "#4267A5",
    "all" = "#777777",
    "todos los DEGs" = "#777777",
    "up" = "#C66C2E",
    "down" = "#4C7A6B"
  )) +
  labs(
    title = "Terminos GO significativos por contraste",
    x = NULL,
    y = "Terminos con FDR < 0.05",
    fill = "Direccion"
  ) +
  theme_functional(10)
ggsave(file.path(figures_dir, "summary", "functional_enrichment_counts.png"), p_counts, width = 10, height = 7, dpi = 200)

grid <- expand.grid(
  contrast = contrast_order,
  category = names(category_label)[names(category_label) != "other"],
  stringsAsFactors = FALSE
)
heat <- merge(grid, category_summary, by = c("contrast", "category"), all.x = TRUE)
heat$short_label <- factor(unname(short_label[heat$contrast]), levels = unname(short_label[contrast_order]))
heat$category_label <- factor(unname(category_label[heat$category]), levels = rev(unname(category_label[names(category_label) != "other"])))
heat$n_terms[is.na(heat$n_terms)] <- 0

p_heat <- ggplot(heat, aes(x = short_label, y = category_label)) +
  geom_point(aes(size = n_terms, fill = dominant_effect), shape = 21, color = "grey30", stroke = 0.2) +
  scale_fill_gradient2(low = "#4267A5", mid = "white", high = "#B6473B", midpoint = 0, na.value = "grey92") +
  scale_size_continuous(range = c(1.5, 8), breaks = c(1, 5, 15, 40, 80)) +
  labs(
    title = "Mapa funcional dirigido basado en GSEA",
    x = NULL,
    y = NULL,
    fill = "NES dominante",
    size = "n terminos"
  ) +
  theme_functional(10) +
  theme(axis.text.y = element_text(size = 9))
ggsave(file.path(figures_dir, "summary", "functional_category_heatmap.png"), p_heat, width = 11, height = 7.5, dpi = 220)

focused_summary <- read_tsv(file.path(project_dir, "tables", "05_de", "development_neural_crest_summary.tsv"))
focused_summary$contrast <- factor(focused_summary$contrast, levels = contrast_order)
focused_summary$category <- factor(focused_summary$category)
focused_summary$short_label <- factor(unname(short_label[as.character(focused_summary$contrast)]), levels = unname(short_label[contrast_order]))
focused_summary$top_NES <- as_num(focused_summary$top_NES)
focused_summary$n_terms_fdr_0_10 <- as_num(focused_summary$n_terms_fdr_0_10)
p_focused <- ggplot(focused_summary, aes(x = short_label, y = category, fill = top_NES)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(n_terms_fdr_0_10 > 0, n_terms_fdr_0_10, "")), size = 3) +
  scale_fill_gradient2(low = "#4267A5", mid = "white", high = "#B6473B", midpoint = 0, na.value = "grey92") +
  labs(
    title = "Categorias dirigidas de desarrollo, migracion, pigmentacion y Wnt/retinoides",
    x = NULL,
    y = NULL,
    fill = "NES termino lider"
  ) +
  theme_functional(9)
ggsave(file.path(figures_dir, "summary", "functional_focused_category_tilemap.png"), p_focused, width = 11, height = 5.8, dpi = 220)

save_empty_plot <- function(path, title) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Sin terminos significativos para esta vista", size = 5) +
    labs(title = title) +
    theme_void()
  ggsave(path, p, width = 8, height = 4, dpi = 180)
}

for (contrast in contrast_order) {
  contrast_fig_dir <- file.path(figures_dir, "contrasts", contrast)
  dir.create(contrast_fig_dir, showWarnings = FALSE, recursive = TRUE)

  gx <- gsea_sig[as.character(gsea_sig$contrast) == contrast, , drop = FALSE]
  g_top <- rbind(
    top_rows(gx[gx$effect > 0, , drop = FALSE], n = 10),
    top_rows(gx[gx$effect < 0, , drop = FALSE], n = 10)
  )
  if (nrow(g_top)) {
    g_top$term_wrapped <- wrap_text(g_top$term_clean, 42)
    g_top$term_wrapped <- factor(g_top$term_wrapped, levels = g_top$term_wrapped[order(g_top$effect)])
    g_top$ontology <- factor(g_top$ontology, levels = c("BP", "CC", "MF"))
    p_g <- ggplot(g_top, aes(x = effect, y = term_wrapped, color = ontology, size = -log10(p_adjust))) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
      geom_point(alpha = 0.92) +
      scale_color_manual(values = c(BP = "#4C7A6B", CC = "#B6473B", MF = "#5E5AA8"), drop = FALSE) +
      labs(
        title = paste0("GSEA por termino GO: ", short_label[contrast]),
        x = "NES",
        y = NULL,
        color = "Ontologia",
        size = "-log10(FDR)"
      ) +
      theme_functional(9) +
      theme(axis.text.x = element_text(angle = 0))
    ggsave(file.path(contrast_fig_dir, paste0("functional_gsea_top_terms_", contrast, ".png")), p_g, width = 9.5, height = 7, dpi = 220)
  } else {
    save_empty_plot(file.path(contrast_fig_dir, paste0("functional_gsea_top_terms_", contrast, ".png")), paste0("GSEA: ", short_label[contrast]))
  }

  ox <- ora_sig[as.character(ora_sig$contrast) == contrast, , drop = FALSE]
  o_top <- do.call(rbind, lapply(c("all_DEGs", "up", "down"), function(direction) {
    yy <- ox[ox$effect_direction == direction, , drop = FALSE]
    if (!nrow(yy)) return(NULL)
    yy <- yy[order(yy$p_adjust, -yy$effect), , drop = FALSE]
    head(yy, 6)
  }))
  if (!is.null(o_top) && nrow(o_top)) {
    o_top$term_wrapped <- wrap_text(o_top$term_clean, 38)
    o_top$direction_label <- factor(
      o_top$effect_direction,
      levels = c("all_DEGs", "up", "down"),
      labels = c("Todos los DEGs", "Up", "Down")
    )
    o_top$term_wrapped <- factor(o_top$term_wrapped, levels = unique(o_top$term_wrapped[order(o_top$signed_effect)]))
    p_o <- ggplot(o_top, aes(x = signed_effect, y = term_wrapped, color = ontology, size = -log10(p_adjust))) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
      geom_point(alpha = 0.92) +
      facet_wrap(~direction_label, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c(BP = "#4C7A6B", CC = "#B6473B", MF = "#5E5AA8"), drop = FALSE) +
      labs(
        title = paste0("ORA por termino GO: ", short_label[contrast]),
        x = "Fold enrichment con signo para Up/Down; zScore para todos",
        y = NULL,
        color = "Ontologia",
        size = "-log10(FDR)"
      ) +
      theme_functional(8.5) +
      theme(axis.text.x = element_text(angle = 0))
    ggsave(file.path(contrast_fig_dir, paste0("functional_ora_top_terms_", contrast, ".png")), p_o, width = 9.5, height = 8, dpi = 220)
  } else {
    save_empty_plot(file.path(contrast_fig_dir, paste0("functional_ora_top_terms_", contrast, ".png")), paste0("ORA: ", short_label[contrast]))
  }
}

message("Generated functional chapter tables in ", tables_dir)
message("Generated functional chapter figures in ", figures_dir)
