### HC: Clostridium, Prevotella 
### MMC: Escherichia coli 
### IHD: Bacteroides vulgatus 
# (hba1c, bmi, LVEF,status(mutliclass:HC reference))
### references: 
# [1] Atarashi, Koji, et al. "Treg induction by a rationally selected mixture of Clostridia strains from the human microbiota." Nature 500.7461 (2013): 232-236.
# [2] Kovatcheva-Datchary, Petia, et al. "Dietary fiber-induced improvement in glucose metabolism is associated with increased abundance of Prevotella." Cell metabolism 22.6 (2015): 971-982.
# [3] Cani, Patrice D., et al. "Metabolic endotoxemia initiates obesity and insulin resistance." Diabetes 56.7 (2007): 1761-1772.
# [4] Jie, Zhuye, et al. "The gut microbiome in atherosclerotic cardiovascular disease." Nature communications 8.1 (2017): 845.
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(car)
library(quantreg)
library(dplyr)
library(ggsignif)
library(hillR)
library(WasserSpectrum)

load("./Real Data/cardiovascular_wgs.RData")
meta <- as(sample_data(ps),"data.frame")
meta$SampleID <- rownames(meta)
meta <- meta %>%
  filter(grepl("^MMC", Status) | grepl("^HC", Status) | grepl("^IHD", Status)) %>%
  mutate(Status = ifelse(grepl("^MMC", Status), "MMC",
                         ifelse(grepl("^HC", Status), "HC", 
                                ifelse(grepl("^IHD", Status), "IHD", NA))))

meta <- meta %>%
  mutate(HbA1c.... = ifelse(HbA1c.... == "NA", NA, HbA1c....)) %>% 
  filter(!is.na(HbA1c....)) %>%
  mutate(HbA1c = as.numeric(HbA1c....))
meta <- meta %>%
  filter(Gender %in% c("Female", "Male"))

ps_filt <- prune_samples(meta$SampleID, ps)
sample_data(ps_filt) <- sample_data(meta)
ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)
meta$Depth <- colSums(otu_table(ps_filt))
ps_filt@sam_data$Age = ps_filt@sam_data$Age..years.
meta$BMI..kg.m.. <- as.numeric(meta$BMI..kg.m..)
meta <- meta[!is.na(meta$BMI..kg.m..), ]
meta$Left.ventricular.ejection.fraction.... <- as.numeric(meta$Left.ventricular.ejection.fraction....)
meta <- meta[!is.na(meta$Left.ventricular.ejection.fraction....), ]
ps_filt <- prune_samples(meta$SampleID, ps)
sample_data(ps_filt) <- sample_data(meta)
# Genus-level aggregation (for HC taxa: Clostridium, Prevotella)
ps_genus <- aggregate_taxa(ps_filt, "Genus")

# Species-level aggregation (for MMC and IHD taxa: Escherichia coli, Bacteroides vulgatus)
ps_species <- aggregate_taxa(ps_filt, "Species")

# Convert to matrix with samples as rows
genus_abund <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) genus_abund <- t(genus_abund)

species_abund <- as(otu_table(ps_species), "matrix")
if (taxa_are_rows(ps_species)) species_abund <- t(species_abund)

# Taxonomy tables for label matching
genus_tax <- tax_table(ps_genus)
species_tax <- tax_table(ps_species)


genus_targets <- c("Clostridium", "Prevotella")
species_targets <- c("Escherichia coli", "Bacteroides vulgatus")

# Genus level
genus_features <- rownames(genus_tax)[genus_tax[, "Genus"] %in% genus_targets]

# Species level
species_features <- rownames(species_tax)[species_tax[, "Species"] %in% species_targets]

# Genus
test_df = data.frame(
  "Clostridium" = as.numeric(otu_table(ps_genus)["Clostridium",]),
  "Prevotella" = as.numeric(otu_table(ps_genus)["Prevotella",]),
  "BMI" = sample_data(ps_genus)$BMI..kg.m..,
  "Status" = sample_data(ps_genus)$Status,
  "Age" = sample_data(ps_genus)$Age..years.,
  "Gender" = sample_data(ps_genus)$Gender,
  "HbA1c" = sample_data(ps_genus)$HbA1c,
  "Depth" = sample_data(ps_genus)$Depth,
  "LVEF" = sample_data(ps_genus)$Left.ventricular.ejection.fraction....
)

## Colstridium
res_Clostridium_BMI <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Clostridium",
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_BMI$circular_plot
ggsave(
  filename = "figures/Colstridium_BMI.eps",
  plot = res_Clostridium_BMI$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)


res_Clostridium_HbA1c <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Clostridium",
  outcome_col = "HbA1c",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_HbA1c$circular_plot
ggsave(
  filename = "figures/Colstridium_HbA1c.eps",
  plot = res_Clostridium_HbA1c$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)


res_Clostridium_LVEF <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Clostridium",
  outcome_col = "LVEF",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_LVEF$circular_plot
ggsave(
  filename = "figures/Colstridium_LVEF.eps",
  plot = res_Clostridium_LVEF$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Clostridium_Status <- wasserstein_spectrum_multiclass(
  df = test_df,
  feature_col = "Clostridium",
  outcome_col = "Status",
  reference_level = "HC",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)

res_Clostridium_Status$plot
ggsave(
  filename = "figures/Colstridium_Status.eps",
  plot = res_Clostridium_Status$plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



## Prevotella

res_Prevotella_BMI <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Prevotella",
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_BMI$circular_plot
ggsave(
  filename = "figures/Prevotella_BMI.eps",
  plot = res_Prevotella_BMI$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Prevotella_HbA1c <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Prevotella",
  outcome_col = "HbA1c",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_HbA1c$circular_plot
ggsave(
  filename = "figures/Prevotella_HbA1c.eps",
  plot = res_Prevotella_HbA1c$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)


res_Prevotella_LVEF <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Prevotella",
  outcome_col = "LVEF",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_LVEF$circular_plot
ggsave(
  filename = "figures/Prevotella_LVEF.eps",
  plot = res_Prevotella_LVEF$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Prevotella_Status <- wasserstein_spectrum_multiclass(
  df = test_df,
  feature_col = "Prevotella",
  outcome_col = "Status",
  reference_level = "HC",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)

res_Prevotella_Status$plot
ggsave(
  filename = "figures/Prevotella_Status.eps",
  plot = res_Prevotella_Status$plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



# Species "Escherichia coli", "Bacteroides vulgatus"
test_df = data.frame(
  "Escherichia_coli" = as.numeric(otu_table(ps_species)["Escherichia coli",]),
  "Bacteroides_vulgatus" = as.numeric(otu_table(ps_species)["Bacteroides vulgatus",]),
  "BMI" = sample_data(ps_species)$BMI..kg.m..,
  "Status" = sample_data(ps_species)$Status,
  "Age" = sample_data(ps_species)$Age..years.,
  "Gender" = sample_data(ps_species)$Gender,
  "HbA1c" = sample_data(ps_species)$HbA1c,
  "Depth" = sample_data(ps_species)$Depth,
  "LVEF" = sample_data(ps_species)$Left.ventricular.ejection.fraction....
)

## Escherichia coli
res_Clostridium_BMI <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Escherichia_coli",
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_BMI$circular_plot
ggsave(
  filename = "figures/Escherichia_coli_BMI.eps",
  plot = res_Clostridium_BMI$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Clostridium_HbA1c <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Escherichia_coli",
  outcome_col = "HbA1c",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_HbA1c$circular_plot
ggsave(
  filename = "figures/Escherichia_coli_HbA1c.eps",
  plot = res_Clostridium_HbA1c$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Clostridium_LVEF <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Escherichia_coli",
  outcome_col = "LVEF",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Clostridium_LVEF$circular_plot
ggsave(
  filename = "figures/Escherichia_coli_LVEF.eps",
  plot = res_Clostridium_LVEF$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Clostridium_Status <- wasserstein_spectrum_multiclass(
  df = test_df,
  feature_col = "Escherichia_coli",
  outcome_col = "Status",
  reference_level = "HC",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)

res_Clostridium_Status$plot
ggsave(
  filename = "figures/Escherichia_coli_Status.eps",
  plot = res_Clostridium_Status$plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)

## Bacteroides vulgatus

res_Prevotella_BMI <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Bacteroides_vulgatus",
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_BMI$circular_plot
ggsave(
  filename = "figures/Bacteroides_vulgatus_BMI.eps",
  plot = res_Prevotella_BMI$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)




res_Prevotella_HbA1c <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Bacteroides_vulgatus",
  outcome_col = "HbA1c",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_HbA1c$circular_plot
ggsave(
  filename = "figures/Bacteroides_vulgatus_HbA1c.eps",
  plot = res_Prevotella_HbA1c$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)


res_Prevotella_LVEF <- wasserstein_spectrum(
  df = test_df,
  feature_col = "Bacteroides_vulgatus",
  outcome_col = "LVEF",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)
res_Prevotella_LVEF$circular_plot
ggsave(
  filename = "figures/Bacteroides_vulgatus_LVEF.eps",
  plot = res_Prevotella_LVEF$circular_plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)



res_Prevotella_Status <- wasserstein_spectrum_multiclass(
  df = test_df,
  feature_col = "Bacteroides_vulgatus",
  outcome_col = "Status",
  reference_level = "HC",
  confounder_cols = c("Age", "Gender", "Depth"),
  plot = TRUE
)

res_Prevotella_Status$plot
ggsave(
  filename = "figures/Bacteroides_vulgatus_Status.eps",
  plot = res_Prevotella_Status$plot,
  device = cairo_ps,
  width = 7,
  height = 5,
  fallback_resolution = 600
)

