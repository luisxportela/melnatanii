pacotes <- c("terra", "readxl", "dplyr", "ggplot2", "ggrepel", "patchwork",
             "ecospat", "usdm", "MASS", "openxlsx")
faltando <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(faltando) > 0) install.packages(faltando)
invisible(lapply(pacotes, library, character.only = TRUE))

set.seed(123)
dir_saida <- "resultados_nicho"
dir.create(dir_saida, showWarnings = FALSE)
arq <- function(x) file.path(dir_saida, x)

# Variables retained by the VIF selection (Supplementary material 2)
usar_vars_publicadas <- TRUE
vars_publicadas <- c("bio_02", "bio_03", "bio_04", "bio_08", "bio_09", "bio_12",
                     "bio_13", "bio_14", "bio_15", "bio_18", "bio_19", "elev")

relacionadas <- c("eichleri", "frigidus", "lhotzkyanus", "salzmannianus")
n_background <- 10000
n_rep_equiv  <- 500

# 1. Climatic layers (select any .tif; all bio and elev layers in the folder are read)
pasta <- dirname(file.choose())
arquivos <- list.files(pasta, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
arquivos <- arquivos[grepl("bio|elev", basename(arquivos), ignore.case = TRUE)]
if (length(arquivos) != 20) stop("Expected 20 layers (bio1-bio19 + elevation); found ", length(arquivos))

rasters <- lapply(arquivos, rast)
ref <- rasters[[1]]
rasters <- lapply(rasters, function(r) {
  if (!compareGeom(r, ref, stopOnError = FALSE)) r <- resample(r, ref, method = "bilinear")
  r
})
raster_stack <- rast(rasters)
nomes <- tolower(basename(arquivos))
names(raster_stack) <- ifelse(grepl("elev", nomes), "elev",
                              sprintf("bio_%02d", suppressWarnings(as.integer(
                                gsub(".*bio[_]?([0-9]{1,2}).*", "\\1", nomes)))))

# 2. Occurrence records (Supplementary material 1)
oco_df <- as.data.frame(read_excel(file.choose(), sheet = 1,
                                   range = cell_limits(c(4, 1), c(NA, 7))))
oco_df <- oco_df[!is.na(oco_df$Record) & !is.na(oco_df$Species), ]
oco_df$especie <- sub("^Mitracarpus ", "", oco_df$Species)

# Study area: extent of all records plus a 1-degree buffer
e0 <- ext(vect(oco_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326"))
e  <- ext(e0[1] - 1, e0[2] + 1, e0[3] - 1, e0[4] + 1)
raster_stack <- crop(raster_stack, e)
gc()

# 3. VIF selection (threshold = 80, Spearman, 1,000 points)
viftest <- vifstep(raster_stack, th = 80, size = 1000, method = "spearman")
write.csv(viftest@results, arq("VIF_results.csv"), row.names = FALSE)
selected_vars <- if (usar_vars_publicadas) vars_publicadas else viftest@results$Variables
raster_stack <- raster_stack[[selected_vars]]

# 4. Extraction and cleaning
casas_decimais <- function(x) {
  s <- sub("0+$", "", formatC(abs(x), format = "f", digits = 8))
  nchar(sub("^[0-9]*\\.?", "", s))
}

pts   <- vect(oco_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
clima <- extract(raster_stack, pts, ID = FALSE)
oco_df <- cbind(oco_df, clima)
oco_df$casas  <- pmax(casas_decimais(oco_df$Longitude), casas_decimais(oco_df$Latitude))
oco_df$celula <- cellFromXY(raster_stack, as.matrix(oco_df[, c("Longitude", "Latitude")]))

oco_df$motivo <- NA_character_
oco_df$motivo[!complete.cases(oco_df[, selected_vars])] <- "No climatic data in the 30 arc-sec cell"
oco_df$motivo[is.na(oco_df$motivo) & oco_df$casas < 2] <-
  "Coordinates given with fewer than two decimal places"
cand <- which(is.na(oco_df$motivo))
dup  <- cand[duplicated(paste(oco_df$especie[cand], oco_df$celula[cand]))]
oco_df$motivo[dup] <- "Same 30 arc-sec cell as another record of the same species"
oco_df$usado <- ifelse(is.na(oco_df$motivo), "Yes", "No")
write.csv(oco_df[, c("Record", "Species", "Longitude", "Latitude", "usado", "motivo")],
          arq("records_used.csv"), row.names = FALSE)

dados <- oco_df[oco_df$usado == "Yes", ]
dados_cong <- dados[dados$especie %in% relacionadas, ]
dados_cong$especie <- factor(dados_cong$especie, levels = relacionadas)
dados_novo <- dados[dados$especie == "elnatanii", ]

# 5. PCA-env calibrated on background points; records projected (Broennimann et al. 2012)
bg <- spatSample(raster_stack, size = n_background, method = "random",
                 na.rm = TRUE, as.df = TRUE)
bg <- bg[complete.cases(bg), selected_vars]

pca_env <- prcomp(bg, center = TRUE, scale. = TRUE)
autoval <- pca_env$sdev^2
var_pct <- 100 * autoval / sum(autoval)

scores_bg <- as.data.frame(predict(pca_env, bg)[, 1:2])
scores    <- as.data.frame(predict(pca_env, dados_cong[, selected_vars])[, 1:2])
scores$especie <- dados_cong$especie
pc_novo <- predict(pca_env, dados_novo[, selected_vars])[1, 1:2]

loadings <- as.data.frame(pca_env$rotation[, 1:2])
loadings$variable <- sub("bio_0?", "bio", rownames(loadings))
contrib  <- data.frame(variable = loadings$variable,
                       PC1 = 100 * loadings$PC1^2, PC2 = 100 * loadings$PC2^2)

# 6. MANOVA and LDA (leave-one-out cross-validation, proportional priors)
dados_manova <- dados_cong[, c(selected_vars, "especie")]
form_manova  <- as.formula(paste0("cbind(", paste(selected_vars, collapse = ", "), ") ~ especie"))
res_manova   <- summary(manova(form_manova, data = dados_manova), test = "Pillai")
st <- res_manova$stats
manova_tab <- data.frame(Effect = "Species", Pillai = st[1, "Pillai"], approx_F = st[1, "approx F"],
                         df_num = st[1, "num Df"], df_den = st[1, "den Df"], P = st[1, "Pr(>F)"])

lda_loo <- lda(especie ~ ., data = dados_manova, CV = TRUE)
conf <- table(Observed = dados_manova$especie, Predicted = lda_loo$class)
acc  <- sum(diag(conf)) / sum(conf)
p_obs  <- rowSums(conf) / sum(conf)
p_pred <- colSums(conf) / sum(conf)
acaso  <- sum(p_obs^2)
pe     <- sum(p_obs * p_pred)
kappa  <- (acc - pe) / (1 - pe)
acerto_sp <- diag(conf) / rowSums(conf)

lda_final  <- lda(especie ~ ., data = dados_manova)
prop_traco <- lda_final$svd^2 / sum(lda_final$svd^2)
coef_ld1   <- lda_final$scaling[, 1]

# Standardized and structure coefficients of LD1 (pooled within-species)
X <- as.matrix(dados_manova[, selected_vars])
grupos <- dados_manova$especie
Sw <- Reduce(`+`, lapply(split(as.data.frame(X), grupos),
                         function(s) cov(s) * (nrow(s) - 1))) / (nrow(X) - nlevels(grupos))
sd_w <- sqrt(diag(Sw))
coef_pad <- coef_ld1 * sd_w
Rw <- cov2cor(Sw)
coef_estrutura <- as.numeric(Rw %*% coef_pad) / sqrt(as.numeric(t(coef_pad) %*% Rw %*% coef_pad))
medias <- aggregate(dados_manova[, selected_vars], by = list(especie = dados_manova$especie), mean)

# 7. Schoener's D and niche equivalency (ecospat), Holm-adjusted P
glob <- rbind(scores_bg[, c("PC1", "PC2")], scores[, c("PC1", "PC2")])
grades <- lapply(relacionadas, function(sp)
  ecospat.grid.clim.dyn(glob = glob, glob1 = glob,
                        sp = scores[scores$especie == sp, c("PC1", "PC2")], R = 100))
names(grades) <- relacionadas

pares <- combn(relacionadas, 2, simplify = FALSE)
resultado_D <- do.call(rbind, lapply(pares, function(par) {
  z1 <- grades[[par[1]]]; z2 <- grades[[par[2]]]
  D  <- ecospat.niche.overlap(z1, z2, cor = TRUE)$D
  set.seed(123)
  eq <- ecospat.niche.equivalency.test(z1, z2, rep = n_rep_equiv, overlap.alternative = "lower")
  data.frame(sp1 = par[1], sp2 = par[2], D = D, P = eq$p.D)
}))
resultado_D$P_Holm <- p.adjust(resultado_D$P, method = "holm")

tab2 <- matrix("–", 4, 4, dimnames = list(paste("M.", relacionadas), paste("M.", relacionadas)))
for (i in seq_len(nrow(resultado_D))) {
  a <- match(resultado_D$sp1[i], relacionadas); b <- match(resultado_D$sp2[i], relacionadas)
  tab2[a, b] <- sprintf("%.3f", resultado_D$D[i])
  tab2[b, a] <- ifelse(resultado_D$P_Holm[i] < 0.001, "< 0.001",
                       sprintf("%.3f", resultado_D$P_Holm[i]))
}
write.csv(resultado_D, arq("Table2_D_P.csv"), row.names = FALSE)
write.csv(tab2, arq("Table2.csv"))

# 8. Position of M. elnatanii (95% ellipses and closest records)
dentro_elipse <- sapply(relacionadas, function(sp) {
  s <- scores[scores$especie == sp, c("PC1", "PC2")]
  mahalanobis(pc_novo, colMeans(s), cov(s)) <= qchisq(0.95, df = 2)
})
dist_reg <- sqrt((scores$PC1 - pc_novo[1])^2 + (scores$PC2 - pc_novo[2])^2)
mais_proximos <- head(data.frame(Record = dados_cong$Record, especie = scores$especie,
                                 distancia = dist_reg)[order(dist_reg), ], 10)
write.csv(mais_proximos, arq("closest_records_elnatanii.csv"), row.names = FALSE)

# 9. Figure 5
cores  <- c(eichleri = "#E1A62F", frigidus = "#1B9E77", lhotzkyanus = "#8C6BB1", salzmannianus = "#D6604D")
formas <- c(eichleri = 16, frigidus = 15, lhotzkyanus = 17, salzmannianus = 18)
rotulos <- setNames(paste("M.", relacionadas), relacionadas)

elipse_95 <- function(x, y, n = 100) {
  if (length(x) < 3) return(NULL)
  centro <- c(mean(x), mean(y)); eig <- eigen(cov(cbind(x, y)))
  ang <- seq(0, 2 * pi, length.out = n)
  pts <- cbind(cos(ang), sin(ang)) %*% diag(sqrt(eig$values * qchisq(0.95, df = 2))) %*% t(eig$vectors)
  data.frame(x = pts[, 1] + centro[1], y = pts[, 2] + centro[2])
}
elipses <- bind_rows(lapply(relacionadas, function(sp) {
  s <- scores[scores$especie == sp, ]
  el <- elipse_95(s$PC1, s$PC2)
  if (!is.null(el)) el$especie <- sp
  el
}))

lim_x <- max(abs(c(scores$PC1, pc_novo[1]))) * 1.08
lim_y <- max(abs(c(scores$PC2, pc_novo[2]))) * 1.08
mult  <- 0.85 * min(lim_x, lim_y) / max(sqrt(loadings$PC1^2 + loadings$PC2^2))

estrela <- function(cx, cy, rx, ry, n = 5) {
  ang <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n * 2 + 1)[-(n * 2 + 1)]
  f <- rep(c(1, 0.4), n)
  data.frame(x = cx + rx * f * cos(ang), y = cy + ry * f * sin(ang))
}
rx <- 0.045 * lim_x; ry <- 0.045 * lim_y
estrela_novo <- estrela(pc_novo[1], pc_novo[2], rx, ry)

tema <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        plot.tag = element_text(face = "bold", size = 16))

p_main <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_polygon(data = elipses, aes(x, y, fill = especie, color = especie), alpha = 0.15, linewidth = 0.5) +
  geom_point(data = scores, aes(PC1, PC2, color = especie, shape = especie), size = 2, alpha = 0.75) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.45) +
  geom_text_repel(data = loadings, aes(PC1 * mult * 1.1, PC2 * mult * 1.1, label = variable),
                  fontface = "italic", size = 3, seed = 123, min.segment.length = Inf,
                  box.padding = 0.2) +
  geom_polygon(data = estrela_novo, aes(x, y), fill = "#C44E52", color = "black", linewidth = 0.5) +
  annotate("text", x = pc_novo[1], y = pc_novo[2] - 2.2 * ry, label = "M. elnatanii",
           fontface = "italic", size = 3.2) +
  scale_color_manual(values = cores, labels = rotulos, name = "Species") +
  scale_fill_manual(values = cores, labels = rotulos, name = "Species") +
  scale_shape_manual(values = formas, labels = rotulos, name = "Species") +
  coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]), tag = "A") +
  tema

p_top <- ggplot(scores, aes(PC1, color = especie, fill = especie)) +
  geom_density(alpha = 0.15, adjust = 1.8) +
  scale_color_manual(values = cores, guide = "none") + scale_fill_manual(values = cores, guide = "none") +
  coord_cartesian(xlim = c(-lim_x, lim_x)) + theme_void()

p_right <- ggplot(scores, aes(PC2, color = especie, fill = especie)) +
  geom_density(alpha = 0.15, adjust = 1.8) +
  scale_color_manual(values = cores, guide = "none") + scale_fill_manual(values = cores, guide = "none") +
  coord_flip(xlim = c(-lim_y, lim_y)) + theme_void()

plot_contrib <- function(eixo, rotulo, tag) {
  cc <- contrib[order(contrib[[eixo]]), ]
  cc$variable <- factor(cc$variable, levels = cc$variable)
  ggplot(cc, aes(.data[[eixo]], variable)) +
    geom_col(fill = "grey30", width = 0.6) +
    geom_vline(xintercept = 100 / nrow(cc), color = "#C44E52", linetype = "dashed") +
    labs(x = rotulo, y = NULL, tag = tag) + tema
}
p_c1 <- plot_contrib("PC1", "Contribution to PC1 (%)", "B")
p_c2 <- plot_contrib("PC2", "Contribution to PC2 (%)", "C")

desenho <- "
AAAAA#
BBBBBC
BBBBBC
BBBBBC
DDDEEE
DDDEEE
"
fig5 <- wrap_plots(A = p_top, B = p_main, C = p_right, D = p_c1, E = p_c2, design = desenho) +
  plot_layout(guides = "collect")
ggsave(arq("Figure5.png"), fig5, width = 8, height = 10, dpi = 600, bg = "white")
ggsave(arq("Figure5.pdf"), fig5, width = 8, height = 10, bg = "white")

# 10. Supplementary materials 3-5
escreve_supl <- function(arquivo, titulo, legenda, tabelas) {
  wb <- createWorkbook()
  negrito <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         border = "TopBottom", wrapText = TRUE)
  for (nome in names(tabelas)) {
    addWorksheet(wb, nome)
    writeData(wb, nome, titulo[[nome]], startRow = 1)
    writeData(wb, nome, legenda[[nome]], startRow = 2)
    mergeCells(wb, nome, cols = 1:ncol(tabelas[[nome]]), rows = 1)
    mergeCells(wb, nome, cols = 1:ncol(tabelas[[nome]]), rows = 2)
    addStyle(wb, nome, createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold"), rows = 1, cols = 1)
    addStyle(wb, nome, createStyle(fontName = "Arial", fontSize = 10, wrapText = TRUE, valign = "top"),
             rows = 2, cols = 1)
    setRowHeights(wb, nome, rows = 2, heights = 75)
    writeData(wb, nome, tabelas[[nome]], startRow = 4, headerStyle = negrito)
    setColWidths(wb, nome, cols = 1:ncol(tabelas[[nome]]), widths = c(24, rep(14, ncol(tabelas[[nome]]) - 1)))
  }
  saveWorkbook(wb, arq(arquivo), overwrite = TRUE)
}

n_total <- nrow(dados_cong)
supl3 <- data.frame(Variable = c(sub("bio_0?", "bio", rownames(pca_env$rotation)) |> sub(pattern = "^elev$", replacement = "Elevation"),
                                 "Variance explained (%)", "Cumulative variance (%)", "Eigenvalue"),
                    rbind(pca_env$rotation, var_pct, cumsum(var_pct), autoval), check.names = FALSE)
escreve_supl("Suppl_material_3_PCA_all_axes.xlsx",
  list(`PCA loadings` = "Supplementary material 3. Loadings of the environmental variables on all axes of the principal component analysis (PCA-env)."),
  list(`PCA loadings` = sprintf(paste(
    "PCA on the correlation matrix of the 12 retained variables, calibrated on %s random background points drawn within the study area",
    "(extent of the occurrence records plus a 1° buffer). The %d records of M. eichleri, M. frigidus, M. lhotzkyanus and M. salzmannianus",
    "and the record of M. elnatanii were projected onto the resulting axes. The sign of each axis is arbitrary."),
    format(nrow(bg), big.mark = ","), n_total)),
  list(`PCA loadings` = supl3))

supl4 <- data.frame(Variable = sub("^elev$", "Elevation", sub("bio_0?", "bio", selected_vars)),
                    `Raw coefficient (LD1)` = coef_ld1, `Standardized coefficient (LD1)` = coef_pad,
                    `Structure coefficient (LD1)` = coef_estrutura, check.names = FALSE)
for (sp in relacionadas) supl4[[paste0("Mean M. ", sp)]] <- as.numeric(unlist(medias[medias$especie == sp, selected_vars]))
manova_out <- data.frame(Effect = "Species", `Pillai's trace` = manova_tab$Pillai,
                         `Approximate F` = manova_tab$approx_F, `df (numerator)` = manova_tab$df_num,
                         `df (denominator)` = manova_tab$df_den,
                         P = ifelse(manova_tab$P < 0.001, "< 0.001", sprintf("%.3f", manova_tab$P)),
                         check.names = FALSE)
escreve_supl("Suppl_material_4_LDA_MANOVA.xlsx",
  list(LDA = "Supplementary material 4. Linear discriminant analysis (LDA) and multivariate analysis of variance (MANOVA).",
       MANOVA = "MANOVA (Pillai's trace)"),
  list(LDA = sprintf(paste(
    "LDA and one-way MANOVA on the 12 retained variables for the %d records of the four related species (one record per species in each",
    "30 arc-second cell; Supplementary material 1). The three discriminant functions account for %.1f%%, %.1f%% and %.1f%% of the",
    "between-species variance (proportion of trace). Raw coefficients depend on the measurement scale of each variable and are not comparable",
    "among variables; standardized coefficients are the raw coefficients multiplied by the pooled within-species standard deviation; structure",
    "coefficients are the pooled within-species correlations between each variable and LD1. The sign of LD1 is arbitrary. Species means are",
    "given in the units of the climatic layers; isothermality (bio3) is expressed as the ratio bio2/bio7."),
    n_total, 100 * prop_traco[1], 100 * prop_traco[2], 100 * prop_traco[3]),
    MANOVA = "One-way MANOVA testing the difference among the four related species in the 12 retained variables."),
  list(LDA = supl4, MANOVA = manova_out))

conf_out <- as.data.frame.matrix(conf)
names(conf_out) <- paste("M.", names(conf_out))
conf_out <- cbind(`Observed \\ Predicted` = paste("M.", rownames(conf_out)), conf_out,
                  Total = rowSums(conf), `Correctly classified (%)` = round(100 * acerto_sp, 1))
conf_out <- rbind(conf_out, data.frame(`Observed \\ Predicted` = "Total", t(colSums(conf)),
                                       Total = sum(conf), `Correctly classified (%)` = round(100 * acc, 1),
                                       check.names = FALSE) |> setNames(names(conf_out)))
escreve_supl("Suppl_material_5_LDA_confusion_matrix.xlsx",
  list(`Confusion matrix` = "Supplementary material 5. Classification of records by linear discriminant analysis with leave-one-out cross-validation."),
  list(`Confusion matrix` = sprintf(paste(
    "Rows: species to which each record was identified; columns: species predicted by the LDA when that record was left out of the training set.",
    "Values on the diagonal are correct classifications. Prior probabilities were proportional to sample sizes; the proportion of records expected",
    "to be correctly classified by chance is %.1f%%, and Cohen's kappa is %.2f."), 100 * acaso, kappa)),
  list(`Confusion matrix` = conf_out))

# 11. Main results
print(resultado_D)
print(manova_tab)
cat(sprintf("LDA accuracy = %.1f%%; kappa = %.2f; chance = %.1f%%\n", 100 * acc, kappa, 100 * acaso))
print(dentro_elipse)
