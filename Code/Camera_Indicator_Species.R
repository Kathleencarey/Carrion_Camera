# Required libraries

library(indicspecies)
library(reshape2)

# Read in data
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

# Run indicator species test 

# Set seed
set.seed(123)

# Run indicator species analysis
indicators <- multipatt(pa_matrix, indices$Design, control = how(nperm = 999))

# View the significant results (p < 0.05)
summary(indicators)

# View the IndVal for every species
summary(indicators, indvalcomp = TRUE, alpha = 1)

