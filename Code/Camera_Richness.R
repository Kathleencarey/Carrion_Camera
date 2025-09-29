# Required libraries
library(lme4)
library(glmmTMB)
library(dplyr)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(ggpubr)

# Load data
Camera_Richness <- read.csv("Camera_Richness.csv")

# Run GLMM
Richness_GLMM <- glmmTMB(
  Richness ~ Design + Month + (1 | Site/Random_Effect),
  data = Camera_Richness,
  family = poisson
)

# View results
summary(Richness_GLMM)


#Tukeys Pairwise Comparison

# Compare by month
Month_Pairwise <- emmeans(Richness_GLMM, ~ Month)

pairs(Month_Pairwise, adjust = "tukey")

# Check model fit
res <- simulateResiduals(fittedModel = Richness_GLMM)
plot(res) 
testUniformity(res) 
testDispersion(res) 
testOutliers(res) 
testZeroInflation(res) 
residual_deviance <- sum(residuals(Richness_GLMM, type = "pearson")^2)
df_residuals <- df.residual(Richness_GLMM)
dispersion_ratio <- residual_deviance / df_residuals
dispersion_ratio 

# Create supplementary species richness table

Richness_Table <- Camera_Richness %>%
  group_by(Month, Design) %>%
  summarise(
    Mean = mean(Richness, na.rm = TRUE),
    Standard_Deviation = sd(Richness, na.rm = TRUE),
    Standard_Error = Standard_Deviation / sqrt(n()),
    Minimum = min(Richness, na.rm = TRUE),
    Maximum = max(Richness, na.rm = TRUE),
    # Unique_Cameras = n_distinct(Camera),  # Removed - add later Count unique cameras
    .groups = 'drop'
  )

# View table

print(Richness_Table)

# Save table
write.csv(Richness_Table, "Table_S1.csv", row.names = FALSE)

# Split Violin Plot 

#Create geom split plot function

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}



# Split-violin Plot by month and design

# Custom colors
custom_colors <- c("H" = "#bf2986ff", "V" = "#6a0dad")
custom_fills <- c("H" = "#eeb4d7ff", "V" = "#c5a3d5")

# Set months to the correct order
Camera_Richness$Month <- factor(Camera_Richness$Month, levels = c("August", "September", "November"))


Split_Violin_Richness <- ggplot() + 
  geom_boxplot(data = Camera_Richness %>% filter(Month == "August" & Design == "V"),
               aes(x = Month, y = Richness, fill = Design, color = Design),
               width = 0.2, alpha = 0.6,
               position = position_nudge(x = 0.1)) +
  stat_summary(data = Camera_Richness %>% filter(Month == "August" & Design == "V"),
               aes(x = Month, y = Richness),
               fun = mean, geom = "point", shape = 18, size = 3,
               position = position_nudge(x = 0.1), color = "black") +
  geom_split_violin(data = Camera_Richness %>% filter(Month == "August" & Design == "H"),
                    aes(x = Month, y = Richness, fill = Design, color = Design),
                    trim = FALSE, alpha = 0.4,
                    position = position_nudge(x = -0.05)) +
  geom_boxplot(data = Camera_Richness %>% filter(Month == "August" & Design == "H"),
               aes(x = Month, y = Richness, fill = Design, color = Design),
               width = 0.2, alpha = 0.6,
               position = position_nudge(x = -0.1)) +
  stat_summary(data = Camera_Richness %>% filter(Month == "August" & Design == "H"),
               aes(x = Month, y = Richness),
               fun = mean, geom = "point", shape = 18, size = 3,
               position = position_nudge(x = -0.10), color = "black") +

  geom_split_violin(data = Camera_Richness %>% 
                      filter(Month %in% c("September", "November")),
                    aes(x = factor(Month), y = Richness, fill = Design, color = Design),
                    trim = FALSE, alpha = 0.4,
                    position = position_dodge(0.25)) +
  geom_boxplot(data = Camera_Richness%>% 
                 filter(Month %in% c("September", "November")),
               aes(x = Month, y = Richness, fill = Design, color = Design,
                   group = interaction(Month, Design)),
               width = 0.2, alpha = 0.6,
               position = position_dodge(0.25)) +
  stat_summary(data = Camera_Richness %>% 
                 filter(Month %in% c("September", "November")),
               aes(x = Month, y = Richness,
                   group = interaction(Month, Design)),
               fun = mean, geom = "point", shape = 18, size = 3,
               position = position_dodge(0.25), color = "black",
               show.legend = TRUE) +

  scale_fill_manual(values = custom_fills,
                    labels = c("H" = "Horizontal", "V" = "Vertical")) +
  scale_color_manual(values = custom_colors,
                     labels = c("H" = "Horizontal", "V" = "Vertical")) +

  labs(x = "", y = "Species Richness", fill = "Design") + 
  theme_pubr() + 
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.60, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_blank())


# View plot

print(Split_Violin_Richness)


# Save plot

ggsave("Figure_2.jpg", plot = Split_Violin_Richness, dpi = 600, height = 6, width = 10)

