# Required libraries

library(dplyr)
library(tidyr)
library(vegan) 
library(reshape2)
library(ggplot2)

# Load data
Camera_Detections <- read.csv("Camera_Detections.csv")



# Create community matrix
Community_Matrix <- dcast(Camera_Detections, Site + Plot + Month + Design ~ Species, 
                          value.var ='Number_of_Individuals', fun.aggregate = sum, na.rm = FALSE)


# Create indices
indices <- Community_Matrix[, 1:4]
abund <- Community_Matrix[, 5:23]
pa_matrix <- decostand(abund, method = "pa")

#Calculate distance 
Community_Dist <- vegdist(pa_matrix, "jaccard")

# Set seed and run NMDS ordination 
set.seed(123)
nmdsCarrion <- metaMDS(Community_Dist, distance = "jaccard", k = 2) 
nmdsCarrion

# Plot NMDS ordination fit
stressplot(nmdsCarrion)

# Set strata to the interaction between site and month
strata_perm <- interaction(indices$Site, indices$Month)

#Run perMANOVA test
set.seed(123)
permCarrion <- adonis2(Community_Dist ~ indices$Design, method = "jaccard",
                       data = indices, 
                       permutations = 999, 
                       strata = strata_perm)

# View results
print(permCarrion)


# Plot NMDS 

# Extract NMDS scores (coordinates) 
nmds_scores <- as.data.frame(scores(nmdsCarrion)) 

# Add the indices to the NMDS scores data frame 
nmds_scores <- cbind(nmds_scores, indices) 

# Custom colors
custom_colors <- c("H" = "#bf2986ff", "V" = "#6a0dad")
custom_fills <- c("H" = "#eeb4d7ff", "V" = "#c5a3d5")


# Extract the stress value
stress_value <- round(nmdsCarrion$stress, 2)  

# Jitter points for visualization of overlapping points
nmds_scores_jitter <- nmds_scores
nmds_scores_jitter$NMDS1 <- jitter(nmds_scores$NMDS1, amount = 0.05)
nmds_scores_jitter$NMDS2 <- jitter(nmds_scores_jitter$NMDS2, amount = 0.05)

# Switch months to chronological
nmds_scores$Month <- factor(nmds_scores$Month, levels = month.name)

# Plot data
Carrion_NMDS <- ggplot(nmds_scores_jitter, aes(x = NMDS1, y = NMDS2, col = Design, fill = Design)) +
  geom_point(data = nmds_scores_jitter, 
             aes(x = NMDS1, y = NMDS2, color = Design), 
             size = 3, alpha = 1) +
  stat_ellipse(aes(group = Design, fill = Design, color = Design), 
               geom = "polygon",    
               alpha = 0.3) +       
  theme_bw() +
  xlim(-2.5, 3.5) +
  ylim(-2.0, 3.0) +
  labs(
    title = NULL,  
    x = "NMDS 1",  
    y = "NMDS 2",
    color = NULL,   
    fill = NULL
  ) +
  scale_color_manual(
    values = custom_colors,   
    labels = c("H" = "Horizontal", "V" = "Vertical")
  ) +
  scale_fill_manual(
    values = custom_fills, 
    labels = c("H" = "Horizontal", "V" = "Vertical")
  ) +
  theme(
    axis.text = element_text(size = 16), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16), 
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 16),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", linewidth = 1),  
    legend.position = c(0.15, 0.95),  # Adjust legend position
    legend.justification = c(1, 1),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  ) + 
  geom_text(data = subset(nmds_scores, Month == "August"),
            aes(x = -2, y = 2, label = paste("Stress =", stress_value)),
            size = 5, color = "black", hjust = 0) +
  facet_wrap(~Month)


# View plot

print(Carrion_NMDS)


# Save plot

ggsave("Figure_3.jpg", plot = Carrion_NMDS, dpi = 600, height = 6, width = 12)

