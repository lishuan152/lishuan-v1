rm(list = ls())
install.packages("vegan")
install.packages("geosphere")
library(vegan)
library(geosphere)
soil_data <- data.frame(
  Sample = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
  Longitude = c(109.292, 109.292, 110.290, 110.286, 111.391, 
                111.665, 109.273, 109.324, 107.111, 106.943),
  Latitude = c(21.5013, 21.5106, 25.0813, 25.0424, 24.4493, 
               24.4243, 24.3549, 24.3817, 23.3429, 23.6733),
  Se_Conc = c(0.29112, 0.81354, 1.00234, 1.54641, 0.32697, 
         0.54163, 0.91645, 1.25108, 0.31808, 0.18304)
)
coords <- cbind(soil_data$Longitude, soil_data$Latitude)
geo_dist <- distm(coords, fun = distHaversine)
rownames(geo_dist) <- soil_data$Sample
colnames(geo_dist) <- soil_data$Sample
geo_dist_km <- geo_dist / 1000
print(round(geo_dist_km, 2))
se_dist <- dist(soil_data$Se_Conc)
se_dist_matrix <- as.matrix(se_dist)
rownames(se_dist_matrix) <- soil_data$Sample
colnames(se_dist_matrix) <- soil_data$Sample
print(round(se_dist_matrix, 4))
mantel_result <- mantel(
  as.dist(geo_dist),
  as.dist(se_dist_matrix), 
  method = "pearson",
  permutations = 9999
)