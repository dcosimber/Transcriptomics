#!/usr/bin/env Rscript

project_dir <- normalizePath(getwd(), mustWork = TRUE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

format_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
}

format_sci <- function(x) {
  ifelse(is.na(x), "", formatC(as.numeric(x), format = "e", digits = 2))
}

tables_dir <- file.path(project_dir, "tables", "05_de")
chapter_dir <- file.path(tables_dir, "chapter_summaries")
by_contrast_dir <- file.path(tables_dir, "by_contrast")

status <- read_tsv(file.path(tables_dir, "deg_contrast_status.tsv"))
sig <- read_tsv(file.path(tables_dir, "deg_significant_all_contrasts_annotated.tsv"))
top50 <- read_tsv(file.path(tables_dir, "deg_top50_by_contrast_annotated.tsv"))
go <- read_tsv(file.path(tables_dir, "go_enrichment_combined_summary.tsv"))

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

contrast_label <- c(
  N_OD_vs_N_BL = "Normal: superficie ocular oscura vs superficie ciega clara",
  P_OD_vs_P_OL = "Pseudoalbino: superficie ocular oscura vs superficie ocular clara",
  P_BD_vs_P_BL = "Pseudoalbino: superficie ciega oscura vs superficie ciega clara",
  P_OL_vs_P_BL = "Pseudoalbino: superficie ocular clara vs superficie ciega clara",
  P_BD_vs_P_OD = "Pseudoalbino: superficie ciega oscura vs superficie ocular oscura",
  P_surface_x_pigmentation = "Pseudoalbino: interaccion superficie ocular/ciega por pigmentacion clara/oscura",
  P_OL_vs_N_OD = "Pseudoalbino: superficie ocular clara vs Normal: superficie ocular oscura",
  P_BD_vs_N_BL = "Pseudoalbino: superficie ciega oscura vs Normal: superficie ciega clara"
)

short_label <- c(
  N_OD_vs_N_BL = "Normal: ocular oscura vs ciega clara",
  P_OD_vs_P_OL = "Pseudoalbino: ocular oscura vs ocular clara",
  P_BD_vs_P_BL = "Pseudoalbino: ciega oscura vs ciega clara",
  P_OL_vs_P_BL = "Pseudoalbino: ocular clara vs ciega clara",
  P_BD_vs_P_OD = "Pseudoalbino: ciega oscura vs ocular oscura",
  P_surface_x_pigmentation = "Pseudoalbino: superficie x pigmentacion",
  P_OL_vs_N_OD = "Pseudoalbino: ocular clara vs Normal: ocular oscura",
  P_BD_vs_N_BL = "Pseudoalbino: ciega oscura vs Normal: ciega clara"
)

status <- status[match(contrast_order, status$contrast), ]
status$label <- unname(contrast_label[status$contrast])
status$short_label <- unname(short_label[status$contrast])

sig$label <- unname(contrast_label[sig$contrast])
sig$short_label <- unname(short_label[sig$contrast])
top50$label <- unname(contrast_label[top50$contrast])
top50$short_label <- unname(short_label[top50$contrast])

quick_read <- c(
  N_OD_vs_N_BL = "Referencia fisiologica normal entre superficie ocular oscura y superficie ciega clara.",
  P_OD_vs_P_OL = "Pseudoalbinismo ocular; mezcla senal de pigmentacion con matriz extracelular y remodelado tisular.",
  P_BD_vs_P_BL = "Pigmentacion ectopica en superficie ciega pseudoalbina; no reproduce limpiamente el eje melanogenico normal.",
  P_OL_vs_P_BL = "Color claro constante; demuestra que las zonas claras pseudoalbinas retienen identidad de superficie.",
  P_BD_vs_P_OD = "Color oscuro constante; separa superficie ciega oscura de superficie ocular oscura con sesgo completo hacia la superficie ciega.",
  P_surface_x_pigmentation = "Interaccion estricta superficie por pigmentacion; lista corta y biologicamente focalizada.",
  P_OL_vs_N_OD = "Superficie ocular clara pseudoalbina frente a superficie ocular oscura normal; perdida parcial del programa normal.",
  P_BD_vs_N_BL = "Superficie ciega oscura pseudoalbina frente a superficie ciega clara normal; contraste mas intenso y mas divergente."
)

summary_table <- data.frame(
  contrast = status$contrast,
  short_label = status$short_label,
  priority = status$priority,
  contrast_class = status$contrast_class,
  n_samples = status$n_samples,
  n_genes_tested = status$n_genes_tested,
  n_degs = status$n_fdr_0_05,
  n_degs_lfc1 = status$n_fdr_0_05_lfc_1_up + status$n_fdr_0_05_lfc_1_down,
  up_or_interaction_positive = status$n_fdr_0_05_lfc_1_up,
  down_or_interaction_negative = status$n_fdr_0_05_lfc_1_down,
  up_fraction = ifelse(status$n_fdr_0_05 > 0, status$n_fdr_0_05_lfc_1_up / status$n_fdr_0_05, NA),
  quick_read = unname(quick_read[status$contrast])
)
summary_table$dominant_direction <- ifelse(
  summary_table$up_or_interaction_positive > summary_table$down_or_interaction_negative,
  "primer grupo / interaccion positiva",
  ifelse(
    summary_table$up_or_interaction_positive < summary_table$down_or_interaction_negative,
    "segundo grupo / interaccion negativa",
    "equilibrado"
  )
)
summary_table$up_fraction_fmt <- format_num(summary_table$up_fraction * 100, 1)
write_tsv(summary_table, file.path(chapter_dir, "chapter_contrast_summary.tsv"))

sig$contrast <- factor(sig$contrast, levels = contrast_order)
sig <- sig[order(sig$contrast, sig$padj, -abs(sig$log2FoldChange)), ]

for (contrast in contrast_order) {
  contrast_sig <- sig[as.character(sig$contrast) == contrast, ]
  contrast_sig$contrast <- as.character(contrast_sig$contrast)
  write_tsv(
    contrast_sig,
    file.path(by_contrast_dir, paste0(contrast, "_significant_annotated.tsv"))
  )
}

select_top <- function(x, n = 6) {
  x <- x[order(x$padj, -abs(x$log2FoldChange)), ]
  head(x, n)
}

top_by_direction <- do.call(
  rbind,
  lapply(contrast_order, function(contrast) {
    cx <- sig[as.character(sig$contrast) == contrast, ]
    directions <- intersect(c("up", "down", "interaction_positive", "interaction_negative"), unique(cx$direction))
    do.call(
      rbind,
      lapply(directions, function(direction) {
        dx <- select_top(cx[cx$direction == direction, ], n = 6)
        if (nrow(dx) == 0) return(NULL)
        dx$chapter_direction <- direction
        dx
      })
    )
  })
)

top_by_direction <- top_by_direction[, c(
  "contrast", "short_label", "chapter_direction", "display_label", "gene_id",
  "baseMean", "log2FoldChange", "padj", "best_description", "annotation_confidence"
)]
top_by_direction$contrast <- as.character(top_by_direction$contrast)
top_by_direction$log2FoldChange_fmt <- format_num(top_by_direction$log2FoldChange, 2)
top_by_direction$padj_fmt <- format_sci(top_by_direction$padj)
top_by_direction$baseMean_fmt <- format_num(top_by_direction$baseMean, 1)
write_tsv(top_by_direction, file.path(chapter_dir, "chapter_top_genes_by_contrast_direction.tsv"))

top_compact <- do.call(
  rbind,
  lapply(contrast_order, function(contrast) {
    select_top(sig[as.character(sig$contrast) == contrast, ], n = 8)
  })
)
top_compact <- top_compact[, c(
  "contrast", "short_label", "direction", "display_label", "gene_id",
  "baseMean", "log2FoldChange", "padj", "best_description"
)]
top_compact$contrast <- as.character(top_compact$contrast)
top_compact$log2FoldChange_fmt <- format_num(top_compact$log2FoldChange, 2)
top_compact$padj_fmt <- format_sci(top_compact$padj)
write_tsv(top_compact, file.path(chapter_dir, "chapter_top_genes_compact.tsv"))

go$contrast <- factor(go$contrast, levels = contrast_order)
go <- go[order(go$contrast, go$method, go$direction, go$`p.adjust`), ]
go_top <- do.call(
  rbind,
  lapply(split(go, list(go$contrast, go$method, go$direction), drop = TRUE), function(x) {
    head(x[order(x$`p.adjust`), ], 5)
  })
)
go_top <- go_top[!is.na(go_top$contrast), ]
go_top$contrast <- as.character(go_top$contrast)
go_top$p_adjust_fmt <- format_sci(go_top$`p.adjust`)
go_top$signed_effect_fmt <- format_num(go_top$signed_effect, 2)
write_tsv(go_top, file.path(chapter_dir, "chapter_go_terms_by_contrast.tsv"))

inventory <- data.frame(
  contrast = contrast_order,
  significant_table = file.path("tables/05_de/by_contrast", paste0(contrast_order, "_significant_annotated.tsv")),
  volcano = file.path("figures/05_de/contrasts", contrast_order, paste0("volcano_", contrast_order, ".png")),
  ma_plot = file.path("figures/05_de/contrasts", contrast_order, paste0("ma_", contrast_order, ".png")),
  heatmap = file.path("figures/05_de/contrasts", contrast_order, paste0("heatmap_", contrast_order, ".png"))
)
write_tsv(inventory, file.path(chapter_dir, "chapter_contrast_inventory.tsv"))

message("Generated DE chapter tables in ", chapter_dir)
