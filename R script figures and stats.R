#Packages necessary
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ggrepel)
library(lubridate)
library(readr)
library(RColorBrewer)
library(ggsci)
library(FactoMineR)
library(stringr)
library(emmeans)
library(forcats)
library(lme4)
library(tibble)
library(broom)
library(ggh4x)

Sys.setlocale("LC_TIME", "C")


###Figure 1###
#A - Soil Occupation 3km around the hives

soiloccup <- read.table("Soil.occupation.3km.csv",sep=";",header=TRUE,dec=",",na.strings="NA")

head(soiloccup)
str(soiloccup)

library(data.table)
long <- melt(setDT(soiloccup), id.vars = "Location", variable.name = "SoilOccupation")

long$Location <- factor(long$Location,  levels = c("Forest","Suburban","Urban","Vineyard"))

p <- long %>% 
  ggplot(aes(x = Location, y = value)) +
  geom_col(aes(fill = SoilOccupation), width = 0.8, col="black") +
  scale_fill_manual(
    values = c("#848586","#4C9F70","#136F63","#ACACDE","#ABDAFC","#ECDD7B"),
    labels = c(
      "Urban fabric", 
      "Grassland", 
      "Forest", 
      "Vineyards", 
      "Water", 
      "Annual crops"
    )
  ) +
  guides(fill = guide_legend(title = "")) +
  ylab("Soil occupation 3 km radius around hives (%)") +
  xlab("Location") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=14),
    axis.text.y = element_text(size=15),
    axis.title.y = element_text(size=15, margin=margin(r=5)),
    axis.title.x = element_text(size=15, vjust=6),
    legend.text = element_text(size=12),
    strip.text.x = element_text(size=15),
    legend.position = 'bottom',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black")
  )
p

#B - PCA
df_metals <- read.table("metals_bees.csv", header = TRUE, sep = ";", dec = ",", na.strings = "NA")
df_soil <- read.table("soil.occupation.3km.csv", header = TRUE, sep = ";", dec = ",", na.strings = "NA")

# Cleaning values <LOD
replace_less_than <- function(x) {
  x <- gsub(",", ".", x)
  ifelse(grepl("^<", x), as.numeric(sub("<\\s*", "", x)) / 3, as.numeric(x))
}
metals_cols <- grep("mg.kg|g.kg", names(df_metals), value = TRUE)

df_metals_clean <- df_metals %>%
  mutate(across(all_of(metals_cols), replace_less_than)) %>%
  mutate(across(all_of(metals_cols), as.numeric))

# Class
metal_classes <- data.frame(
  Colname = metals_cols,
  Classe = c("Trace elements", "Trace elements", "Trace elements", "Macroelements", "Trace elements",
             "Trace elements", "Trace elements", "Trace elements", "Macroelements", "Trace elements",
             "Trace elements", "Macroelements", "Trace elements", "Macroelements",
             "Macroelements", "Macroelements", "Trace elements")
)

# Select only trace elements
trace_cols <- metal_classes %>% filter(Classe == "Trace elements") %>% pull(Colname)
# or select all
df_traces <- df_metals_clean %>%
  select(Location, all_of(trace_cols))

# Join with soil occupation
df_full <- df_traces %>% left_join(df_soil, by = "Location")
# Prepare for PCA
df_acp <- df_full %>% select(-Location)

# avec FactoMineR
res.pca <- PCA(df_acp, scale.unit = TRUE, graph = FALSE)

# Coordinates indiv.
ind_coords <- as.data.frame(res.pca$ind$coord)
ind_coords$Location <- df_full$Location
ind_coords$Location <- factor(ind_coords$Location, levels = c("Urban", "Suburban", "Forest", "Vineyard"))

# Coordinates
var_coords <- as.data.frame(res.pca$var$coord)
var_coords$Variable <- rownames(var_coords)
var_coords$Variable <- var_coords$Variable %>% str_replace_all("\\..*", "")  # remove units

# Calculate centers
offset <- 1.2
centers <- ind_coords %>%
  group_by(Location) %>%
  summarise(mean_Dim1 = mean(Dim.1), mean_Dim2 = mean(Dim.2)) %>%
  rowwise() %>%
  mutate(
    dec_Dim1 = ifelse(mean_Dim1 >= 0, mean_Dim1 + offset, mean_Dim1 - offset),
    dec_Dim2 = ifelse(mean_Dim2 >= 0, mean_Dim2 + offset, mean_Dim2 - offset)
  )

# Préparation des abréviations éléments trace à afficher
abbr_table <- c(
  "Cadmium" = "Cd", "Cobalt" = "Co", "Copper" = "Cu", "Iron" = "Fe",
  "Manganese" = "Mn", "Lead" = "Pb", "Zinc" = "Zn", 
  "Arsenic" = "As", "Chromium" = "Cr", "Aluminium" = "Al", "Calcium" = "Ca", "Magnesium" = "Mg", "Molybdenum" = "Mo")

var_coords$Variable <- var_coords$Variable %>% str_replace_all("\\..*", "")
var_coords$Abbr <- ifelse(var_coords$Variable %in% names(abbr_table),
                          abbr_table[var_coords$Variable],
                          var_coords$Variable)
var_coords$Abbr_lab <- paste0('bold("', var_coords$Abbr, '")')

# Classif. variables
var_coords$Type <- ifelse(var_coords$Variable %in% str_replace_all(trace_cols, "\\..*", ""),
                          "Trace elements", "Soil occupation")

# Longer arrows
arrow_length_factor <- 2.5
var_coords$Dim.1 <- var_coords$Dim.1 * arrow_length_factor
var_coords$Dim.2 <- var_coords$Dim.2 * arrow_length_factor

# Show or not variables soil occupation
show_soil_occupation <- FALSE  # TRUE show all, FALSE only trace elements
if (!show_soil_occupation) {
  var_plot <- var_coords %>% filter(Type == "Trace elements")
} else {
  var_plot <- var_coords
}

# Colors
location_colors <- c(
  "Urban" = "#d95f02",
  "Suburban" = "#e7298a",
  "Forest" = "#1b9e77",
  "Vineyard" = "#7570b3"
)
variable_colors <- c(
  "Trace elements" = "black",
  "Soil occupation" = "forestgreen"
)
location_shapes <- c(
  "Urban"    = 16,  # circle
  "Suburban" = 17,  # triangle
  "Forest"   = 15,  # square
  "Vineyard" = 18   # diamond
)

#calculation % PCA
percent_var <- res.pca$eig[1:2, 2]

#PLOT
ggplot() +
  geom_point(data = ind_coords,
             aes(x = Dim.1, y = Dim.2, color = Location, shape = Location),
             size = 3) +
  stat_ellipse(data = ind_coords,
               aes(x = Dim.1, y = Dim.2, fill = Location, color = Location),
               type = "norm", geom = "polygon", alpha = 0.12, linewidth = 0.8,
               show.legend = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_segment(data = var_plot,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.35, "cm"), type = "open"),
               color = "black", linewidth = 0.7,
               show.legend = FALSE) +
  geom_text_repel(data = var_plot,
                  aes(x = Dim.1, y = Dim.2, label = Abbr_lab),
                  size = 5, fontface = "bold", parse = TRUE,
                  nudge_x = 0.1 * sign(var_plot$Dim.1),
                  nudge_y = 0.1 * sign(var_plot$Dim.2),
                  color = "black",
                  show.legend = FALSE) +
  geom_text_repel(data = centers,
                  aes(x = dec_Dim1, y = dec_Dim2, label = Location, color = Location),
                  size = 5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(
    values = location_colors,
    breaks = names(location_colors),
    labels = names(location_colors),
    name = "Location"
  ) +
  scale_fill_manual(
    values = location_colors,
    breaks = names(location_colors),
    labels = names(location_colors),
    name = "Location"
  ) +
  scale_shape_manual(
    values = location_shapes,
    breaks = names(location_shapes),
    labels = names(location_shapes),
    name = "Location",
  ) +
  guides(
    shape = "none",          
    fill  = "none"           
  ) +
  labs(
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
    color = "Location",
    fill = "Location",
    shape = "Location"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


################################################################################
###Figure 2###

df_metals <- read.table("metals_bees.csv", header = TRUE, sep = ";", dec = ",", na.strings = "NA")
df_metals$Date <- as.Date(df_metals$Date, format = "%d/%m/%Y")
df_metals$Date <- factor(df_metals$Date)

# Cleaning values <LOD
replace_less_than <- function(x) {
  x <- gsub(",", ".", x)
  x <- ifelse(grepl("^<", x),
              as.numeric(sub("<\\s*", "", x)) / 3,
              as.numeric(x))
  return(x)
}
cols_with_less_than <- c("Aluminium.mg.kg", "Arsenic.mg.kg", "Cadmium.mg.kg",
                         "Chromium.mg.kg", "Cobalt.mg.kg", "Lead.mg.kg")
df_metals_clean <- df_metals %>%
  mutate(across(all_of(cols_with_less_than), replace_less_than)) %>%
  mutate(across(-c(Hive, Location, Date), ~ as.numeric(gsub(",", ".", .))))

# Metal names
format_metal_label <- function(metal_name) {
  metal_name %>%
    str_remove("\\.g\\.kg|\\.mg\\.kg") %>%
    str_replace_all("\\.", " ") %>%
    str_to_title()
}
df_metals_long <- df_metals_clean %>%
  pivot_longer(cols = -c(Date, Hive, Location),
               names_to = "Metal_raw", values_to = "Concentration") %>%
  mutate(Metal = format_metal_label(Metal_raw))
df_metals_long$Location <- as.factor(df_metals_long$Location)
levels(df_metals_long$Location) <- gsub("[()]", "", levels(df_metals_long$Location))
levels(df_metals_long$Location) <- gsub(" ", "_", levels(df_metals_long$Location))

df_metals_long <- df_metals_long %>%
  mutate(LogConcentration = log(Concentration + 1e-6))

#LMER mixed model
modele_mixte <- lmer(LogConcentration ~ Metal * Location + (1 | Date), data = df_metals_long)

# Contrasts Vineyard vs others
contrast_vin <- list(
  "Vineyard vs Forest"    = c(-1, 0, 0, 1),
  "Vineyard vs Sub_urban" = c(0, -1, 0, 1),
  "Vineyard vs Urban"     = c(0, 0, -1, 1)
)
# Contrasts Urban/Suburban/Forest
contrast_other <- list(
  "Urban vs Sub_urban" = c(0, -1, 1, 0),
  "Urban vs Forest" = c(-1, 0, 1, 0),
  "Sub_urban vs Forest" = c(-1, 1, 0, 0)
)
metaux <- unique(df_metals_long$Metal)

results_cmp <- lapply(metaux, function(m) {
  emm_metal <- emmeans(modele_mixte, ~ Location, at = list(Metal = m))
  contr <- contrast(emm_metal, method = contrast_vin, adjust = "none")
  dfc <- as.data.frame(contr)
  dfc$Metal <- m
  return(dfc)
}) %>% bind_rows() %>%
  mutate(Significant = p.value < 0.05,
         Metal = factor(Metal, levels = rev(sort(unique(Metal))))
  )

results_extra_cmp <- lapply(metaux, function(m){
  emm_metal <- emmeans(modele_mixte, ~ Location, at = list(Metal = m))
  contr <- contrast(emm_metal, method = contrast_other, adjust = "none")
  dfc <- as.data.frame(contr)
  dfc$Metal <- m
  return(dfc)
}) %>% bind_rows() %>%
  mutate(Significant = p.value < 0.05,
         Metal = factor(Metal, levels = rev(sort(unique(Metal))))
  )

# All in a data.frame for facet
results_all <- bind_rows(
  results_cmp,
  results_extra_cmp
) %>%
  mutate(
    Contrast_clean = factor(
      contrast,
      levels = c(
        "Vineyard vs Forest",
        "Vineyard vs Sub_urban",
        "Vineyard vs Urban",
        "Urban vs Sub_urban",
        "Urban vs Forest",
        "Sub_urban vs Forest"
      )
    )
  )

# Plot facet 2 lines 3 columns
ggplot(results_all, aes(x = estimate, y = Metal)) +
  geom_point(aes(color = Significant), size = 3) +
  geom_errorbarh(aes(xmin = estimate - 1.96 * SE, xmax = estimate + 1.96 * SE),
                 height = 0.15, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c("TRUE" = "#B22222", "FALSE" = "black"),  # black signif., firebrick non signif.
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
    name = "Significance"
  ) +
  facet_wrap(~ Contrast_clean, ncol = 3, scales = "free_x") +
  labs(x = "Estimates (log-difference)", y = "Element") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(face = "italic", size = 12),
    strip.background = element_blank(),  # supprimer le fond
    strip.text = element_text(face = "bold", size = 12, color = "black", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    panel.border = element_blank(),   # supprime toute bordure
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90")
  )

################################################################################
###Figure 3###

df_metals <- read.table("metals_bees.csv", header = TRUE, sep = ";", dec = ",", na.strings = "NA")

df_metals$Date <- as.Date(df_metals$Date, format = "%d/%m/%Y")

# Fonction pour gérer les valeurs "<"
replace_less_than <- function(x) {
  x <- gsub(",", ".", x)
  x <- ifelse(grepl("^<", x),
              as.numeric(sub("<\\s*", "", x)) / 3,
              as.numeric(x))
  return(x)
}

# Clean values <LOD
cols_with_less_than <- c("Aluminium.mg.kg", "Arsenic.mg.kg", "Cadmium.mg.kg",
                         "Chromium.mg.kg", "Cobalt.mg.kg", "Lead.mg.kg")
# Clean
df_metals_clean <- df_metals %>%
  mutate(across(all_of(cols_with_less_than), replace_less_than)) %>%
  mutate(across(-c(Hive, Location, Date), ~ as.numeric(gsub(",", ".", .))))
# Names
format_metal_label <- function(metal_name) {
  metal_name %>%
    str_remove("\\.g\\.kg|\\.mg\\.kg") %>%
    str_replace_all("\\.", " ") %>%
    str_to_title()
}

# Pivot longer
df_metals_long <- df_metals_clean %>%
  pivot_longer(cols = -c(Date, Hive, Location),
               names_to = "Metal_raw", values_to = "Concentration") %>%
  mutate(Metal = format_metal_label(Metal_raw))

# Classes
metal_classes <- data.frame(
  Metal = c("Aluminium", "Arsenic", "Cadmium", "Calcium", "Chromium",
            "Cobalt", "Copper", "Iron", "Magnesium", "Manganese",
            "Molybdenum", "Phosphorous", "Lead", "Potassium",
            "Sodium", "Sulfur", "Zinc"),
  Classe = c("Trace elements", "Trace elements", "Trace elements", "Macroelements", "Trace elements",
             "Trace elements", "Trace elements", "Trace elements", "Macroelements", "Trace elements",
             "Trace elements", "Macroelements", "Trace elements", "Macroelements",
             "Macroelements", "Macroelements", "Trace elements")
)


df_metals_long_classed <- df_metals_long %>%
  left_join(metal_classes, by = "Metal")

# Mean of two hives per Date, Location, Metal
df_all_avg <- df_metals_long_classed %>%
  group_by(Date, Location, Metal) %>%
  summarise(
    Concentration = mean(Concentration, na.rm = TRUE),  
    n = n(),  
    .groups = "drop"
  ) %>%
  mutate(
    SD = df_metals_long_classed %>% 
      group_by(Date, Location, Metal) %>% 
      summarise(SD = sd(Concentration, na.rm = TRUE)) %>%
      pull(SD), 
    SE = ifelse(n > 1, SD / sqrt(n), NA_real_)  # SE if more than 1 observation, otherwise NA
  )

#Only Vineyards
df_vin <- df_all_avg %>% 
  filter(Location == "Vineyard", Metal %in% c("Copper","Sulfur","Aluminium", "Arsenic","Cobalt","Chromium")) %>%
  mutate(
    Date  = as.Date(Date),
    Month = floor_date(Date, "month"),
    Metal = factor(Metal,
                   levels = c("Copper","Sulfur","Aluminium","Arsenic","Cobalt","Chromium"),
                   labels = c("Copper (mg/kg)","Sulfur (g/kg)","Aluminium (mg/kg)","Arsenic (mg/kg)","Cobalt (mg/kg)","Chromium (mg/kg)")),
    Month_num = interval(start = min(Date), end = Date) %/% months(1)
  )

# solpes and p-values
slopes_tbl <- df_vin %>%
  group_by(Metal) %>%
  do(tidy(lm(Concentration ~ Month_num, data = .))) %>%
  filter(term == "Month_num") %>%
  select(Metal, estimate, p.value) %>%
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    label = sprintf("Slope = %.2f\np = %.3g %s", estimate, p.value, stars)
  )

# y_max_tbl
y_max_tbl <- df_vin %>%
  group_by(Metal) %>%
  summarise(y_max = max(Concentration, na.rm = TRUE), .groups = "drop")

# Phyto treaments applied on nearest plot
df_events_all <- tribble(
  ~md,     ~Metal,      ~Dose,
  "04-20", "Copper",     0.75,
  "05-20", "Copper",     2,
  "05-28", "Copper",     1,
  "06-10", "Copper",     1,
  "06-21", "Copper",     1.2,
  "07-05", "Copper",     1.2,
  "07-15", "Copper",     1,
  "07-29", "Copper",     2.5,
  "04-20", "Sulfur",     4,
  "05-12", "Sulfur",     3,
  "05-28", "Sulfur",     3,
  "06-10", "Sulfur",     2,
  "06-21", "Sulfur",     3,
  "07-05", "Sulfur",     4,
  "07-15", "Sulfur",     4,
  "07-29", "Sulfur",     3,
  "05-28", "Aluminium",  2,
  "07-29", "Aluminium",  2.5
) %>%
  mutate(
    Metal   = factor(Metal,
                     levels = c("Copper","Sulfur","Aluminium","Arsenic","Cobalt","Chromium"),
                     labels = c("Copper (mg/kg)","Sulfur (g/kg)","Aluminium (mg/kg)","Arsenic (mg/kg)","Cobalt (mg/kg)","Chromium (mg/kg)")),
    Month   = as.Date(paste0("2021-", md)),
    label_x = Month + days(2)
  ) %>%
  left_join(y_max_tbl, by = "Metal") %>%     
  mutate(
    label_y = case_when(
      Metal == "Sulfur (g/kg)"     ~ y_max * 1.35,
      Metal == "Copper (mg/kg)"     ~ y_max * 1.30,
      Metal == "Aluminium (mg/kg)"     ~ y_max * 1.30)
  )

# Annotation p-value
annot_pvalues <- slopes_tbl %>%
  left_join(y_max_tbl, by = "Metal") %>%
  mutate(
    x = as.Date("2021-09-15"),
    y = y_max * 1.4        
  )

# Plot
ggplot(df_vin, aes(x = Date, y = Concentration, color = Metal)) +
  geom_pointrange(aes(ymin = Concentration - SE, ymax = Concentration + SE),
                  position = position_dodge(width = 2), fatten = 2, size = 0.8) +
  geom_smooth(aes(x = Month, fill = Metal),
              method    = "lm", se = TRUE, fullrange = TRUE,
              size = 1, alpha = 0.1) +
  geom_vline(data = df_events_all,
             aes(xintercept = as.numeric(Month), color = Metal),
             linetype = "dashed") +
  geom_text(data = df_events_all,
            aes(x = label_x, y = label_y, label = Dose),
            angle = 90, hjust = 0.5, vjust = 0.5, size = 3) +
  geom_text(data = annot_pvalues,
            aes(x = x, y = y, label = label, color = Metal),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~Metal, scales = "free_y", ncol = 2, dir = "v") +
  scale_x_date(
    date_labels = "%b",
    date_breaks = "1 month",
    limits      = c(as.Date("2021-04-01"), max(df_vin$Month) + months(1)),
    expand      = expansion(add = c(0, 0))
  ) +
  scale_color_manual(values = c(
    "Copper (mg/kg)"   = "#B77729",
    "Sulfur (g/kg)"     = "#fdcc36",
    "Aluminium (mg/kg)" = "#888B8D",
    "Arsenic (mg/kg)" = "#80b692",
    "Cobalt (mg/kg)" = "#0050B5",
    "Chromium (mg/kg)" = "#8685ab"
  )) +
  scale_fill_manual(values = c(
    "Copper (mg/kg)"   = "#B77729",
    "Sulfur (g/kg)"     = "#fdcc36",
    "Aluminium (mg/kg)" = "#888B8D",
    "Arsenic (mg/kg)" = "#80b692",
    "Cobalt (mg/kg)" = "#0050B5",
    "Chromium (mg/kg)" = "#8685ab"
  )) +
  #coord_cartesian(ylim = c(0, NA)) +
  labs(
    #title    = "Monthly evolution of concentrations and treatments (Vineyard)",
    #subtitle = "Points = mean ± SE; dashed lines = treatment dose (kg/ha)",
    x        = "Month",
    y        = "Concentration",
    color    = "Compound",
    fill     = "Compound"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    #panel.grid.major.y = element_line(color = "grey80", size = 0.5),
    panel.grid.minor.y = element_blank(), 
    panel.background = element_blank(),   
    axis.line = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks.x = element_line(color = "grey30"),
    axis.title.y = element_text(size=14, margin=margin(r=5)),
    axis.title.x = element_text(size=14),
    strip.text.x = element_text(size=13),
    legend.position   = "none")


################################################################################
###Figure 4###

df <- read.table(file = "pesticide_bees.csv", dec = ",", sep = ";", header = TRUE, na.strings = "NA")
df$Date <- as.Date(df$Date, format = "%d/%m/%Y")

# Change format and handling names and NA
df_long <- df %>%
  pivot_longer(
    cols = starts_with("Pest_"),
    names_to = "Pesticide",
    values_to = "Concentration"
  ) %>%
  mutate(
    Hive_Location = paste(Location, Hive, sep = "_"),
    Concentration = ifelse(is.na(Concentration), 0, Concentration),
    Pesticide = sub("^Pest_", "", Pesticide)
  )

# Concentration max Pesticide/Date/Location
df_max <- df_long %>%
  group_by(Location, Date, Pesticide) %>%
  summarise(Max_Concentration = max(Concentration, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Date_formatted = factor(format(Date, "%b %d"), levels = format(sort(unique(Date)), "%b %d")),
    Conc_log = log10(Max_Concentration + 1))

# Classes
pesticide_classes <- data.frame(
  Pesticide = c(
    "DMF", "ametoctradin", "cyflufenamid", "cymiazole", "difenoconazole",
    "dimethomorph", "fenbuconazole", "fludioxonil", "fluopicolide",
    "fluopyram", "mandipropamid", "metrafenone",
    "pyrimethanil", "pyraclostrobin", "tebuconazole", "trifloxystrobin", "zoxamide"
  ),
  Classe = c(
    "Acaricides", "Fungicides", "Fungicides", "Acaricides", "Fungicides",
    "Fungicides", "Fungicides", "Fungicides", "Fungicides",
    "Fungicides", "Fungicides", "Fungicides",
    "Fungicides", "Fungicides", "Fungicides", "Fungicides", "Fungicides"
  )
)
# Join to table
df_max <- df_max %>%
  left_join(pesticide_classes, by = "Pesticide") %>%
  filter(!is.na(Classe)) %>%
  group_by(Pesticide) %>%
  filter(any(Max_Concentration > 0)) %>%
  ungroup() %>%
  mutate(
    Classe = factor(Classe, levels = c("Fungicides", "Acaricides")),
    Pesticide = fct_reorder(Pesticide, as.numeric(Classe))
  )

# Labels
class_labels <- c(
  "Fungicides" = "Fungicides",
  "Acaricides" = "Acaricides"
)

#Plot
ggplot(df_max, aes(x = Date_formatted, y = Pesticide, fill = Conc_log)) +
  geom_tile(color = "grey90") +
  facet_nested(
    rows = vars(Classe),
    cols = vars(Location),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(Classe = class_labels),
    switch = "y",
    strip = strip_themed(
      background_y = list(
        "Fungicides" = element_rect(fill = NA, color = NA),
        "Acaricides" = element_rect(fill = NA, color = NA)
      ),
      text_y = list(
        "Fungicides" = element_text(color = "black", size = 14, face = 'bold'),
        "Acaricides" = element_text(color = "black", size = 14, face='bold')
      ),
      background_x = element_blank(),
      text_x = element_text(face = "bold", size=18)
    )
  ) +
  scale_fill_gradientn(
    colours = c("grey90", "yellow", "orange", "red", "darkred", "purple4"),
    values = scales::rescale(c(0, 0.1, 0.5, 1, 1.5, 2.5)),
    name = "log10(Conc + 1)",
    na.value = "grey80"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size=13),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=13, margin=margin(r=5)),
    axis.title.x = element_text(size=13, vjust=6),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    strip.text = element_text(face = "bold", size = 14, color = "black", hjust = 0.5),
    panel.grid = element_blank(),
    strip.placement = "outside",
    legend.position = "right"
  ) +
  labs(x = "", y = "")


################################################################################
###Figure 5###

df <- read.csv2(file = "pesticide_bees.csv",
                 dec = ",", sep = ";", header = TRUE, na.strings = "NA")
df$Date <- as.Date(df$Date, format = "%d/%m/%Y")

df <- df %>%
  # keep only Vineyard
  filter(Location == "Vineyard") %>%
  # drop DMF and cymiazole columns because acaricides do not come from environment
  select(-Pest_DMF, -Pest_cymiazole) %>%
  # remove "Pest_" prefix from remaining pesticide columns
  rename_with(~ sub("^Pest_", "", .x), starts_with("Pest_"))


molecules <- names(df)[5:ncol(df)]
for(mol in molecules) df[[mol]] <- as.numeric(gsub(",", ".", df[[mol]]))
#df$Date <- dmy(df$Date)

df_long <- df %>%
  select(Date, Hive, all_of(molecules)) %>%
  pivot_longer(cols = all_of(molecules), names_to = "Molecule", values_to = "Concentration_ng_g")

df_max <- df_long %>%
  group_by(Date, Molecule) %>%
  summarise(Max_Conc_ng_g = if(all(is.na(Concentration_ng_g))) NA_real_ else max(Concentration_ng_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(Max_Conc_ng_g = ifelse(is.na(Max_Conc_ng_g), 0, Max_Conc_ng_g))

traitements <- read_csv2("treatments_pesticides.csv")
traitements$Date <- dmy(traitements$Date)
traitements <- traitements %>% filter(!is.na(Date))

common_mols <- intersect(unique(df_max$Molecule), unique(traitements$Molecule))
only_df <- setdiff(unique(df_max$Molecule), unique(traitements$Molecule))
only_trait <- setdiff(unique(traitements$Molecule), unique(df_max$Molecule))

fixed_colors <- brewer.pal(8, "Set3")
fixed_colors <- setNames(fixed_colors[seq_along(common_mols)], common_mols)

set.seed(42)
random_colors_df <- setNames(grDevices::rainbow(length(only_df)), only_df)
random_colors_trait <- setNames(grDevices::rainbow(length(only_trait)), only_trait)

mycolors <- c(fixed_colors, random_colors_df, random_colors_trait)
all_levels <- names(mycolors)

df_max$Molecule <- factor(df_max$Molecule, levels = all_levels)
traitements$Molecule <- factor(traitements$Molecule, levels = all_levels)

df_max <- df_max %>% mutate(is_treatment = Molecule %in% traitements$Molecule)
df_max <- df_max %>% mutate(is_not_treatment = !(Molecule %in% traitements$Molecule))

delta_days <- 0.5  
traitements <- traitements %>%
  group_by(Date) %>%
  arrange(Molecule) %>%
  mutate(pos = row_number(),
         offset = -((pos - 1) * delta_days),
         x = as.numeric(Date) + offset,
         x_date = as.Date(x, origin = "1970-01-01"))

ymax <- max(df_max$Max_Conc_ng_g, na.rm = TRUE) * 1.1

df_max <- df_max %>%
  group_by(Molecule) %>%
  mutate(pos = as.numeric(factor(Molecule))) %>%
  ungroup() %>%
  mutate(Date_num = as.numeric(Date),
         Date_shifted = as.Date(Date_num + (pos - mean(pos)) * 0.5, origin = "1970-01-01"))

ymax <- df_max %>%
  group_by(Date) %>%
  summarise(sum_conc = sum(Max_Conc_ng_g, na.rm=TRUE)) %>%
  summarise(ymax = max(sum_conc)) %>%
  pull(ymax)
ymax <- ymax * 1.1  # rajouter marge

#Plot
ggplot(df_max, aes(x = Date_shifted, y = Max_Conc_ng_g, fill = Molecule, pattern = is_not_treatment)) +
  geom_col_pattern(aes(pattern = is_not_treatment),
                   position = "stack", color = "black", width = 6,
                   pattern_fill = "black", pattern_angle = 45,
                   pattern_density = 0.3, pattern_spacing = 0.01,
                   pattern_size = 0.2, pattern_linetype = 1) +
  geom_segment(data = traitements,
               aes(x = x_date, xend = x_date, y = 0, yend = ymax, color = Molecule),
               inherit.aes = FALSE,
               linetype = "dashed", linewidth = 1.3) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  scale_pattern_manual(
    name = "Substance use",
    values = c("TRUE" = "stripe", "FALSE" = "none"),
    labels = c("TRUE" = "Treatment", "FALSE" = "No treatment")
  ) +
  labs(pattern = NULL) +
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill = c("white", "white"),     
        color = "black",                  
        pattern_fill = "black",         
        pattern = c("stripe", "none")
      ),
      title = "Treatment indicator"
    ),
    fill = guide_legend(override.aes = list(pattern = "none"))
  ) +
  scale_x_date(
    breaks = unique(df_max$Date),
    labels = scales::date_format("%b %d"),
    expand = expansion(mult = c(0.02, 0.02))) +
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill = c("white", "white"),     
        color = "black",                  
        pattern_fill = "black",         
        pattern = c("stripe", "none")
      ),
      title = "Treatment indicator",
      ncol = 1,
      title.position = "top"
    ),
    fill = guide_legend(
      override.aes = list(pattern = "none"),
      ncol = 4,
      title.position = "top"
    ),
    color = guide_legend(
      ncol = 2,
      title.position = "top")
      )+
  labs(title = "",
       x = "Sampling date", y = "Concentration (ng/g)",
       fill = "Pesticide residues in bees", color = "Pesticide treatments") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks.x = element_line(color = "grey30"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "lines")
  )+
  coord_cartesian(ylim = c(0, ymax))



