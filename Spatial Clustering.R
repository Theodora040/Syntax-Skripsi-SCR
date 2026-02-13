## ANALISIS TINGKAT KEMISKINAN DI PROVINSI JAWA BARAT
## MENGGUNAKAN METODE SPATIAL CLUSTERING REGRESSION

library(maptools)
library(spdep)
library(sp)
library(RColorBrewer)
library(lattice)
library(gstat)
library(raster)
require(splancs)
library(dplyr)
library(tidyr)
library(extrafont) 
library(SpatialEpi)
library(rgdal)
library(maps)
library(maptools)
library(ggplot2)
library(grid)
library(rgdal)
library(colourvalues)
library(viridis) 
library(spatialreg)
library(tidyverse)

setwd("D:\\Kuliah\\SKRIPSI LANCAR\\SIDANG\\ANALISIS")
library(readxl)

# Analisis Deskriptif #
data <- read_excel("Data Penelitian.xlsx")
head(data)
summary(data)

# Spatial Clustering #
tingkat_miskin <- read_excel("Data Kemiskinan.xlsx")
head(tingkat_miskin)

# Matriks Bobot Spasial #
library(spdep)
library(sf)

jabar=readShapePoly("D:\\Kuliah\\SKRIPSI LANCAR\\SIDANG\\Jabar.shp")
jabar@data=tingkat_miskin
jabar
ID<-c(1:27)
jabar$ID<-c(1:27)

# Matriks Pembobot 
C <- poly2nb(jabar, row.names=jabar$ID, queen=TRUE) 
WC <- nb2mat(C, style='B', zero.policy =TRUE) #menyajikan dalam bentuk matrix kontiguitas = unstandardized (1,0)
WC[1:27,1:27]
WW <- nb2mat(C, style='W', zero.policy =TRUE) #menyajikan dalam bentuk matrix pembobot spasial = standardized
WW[1:27,1:27]
W <-nb2listw(C)
summary(W)

# Pengujian Autokorelasi Spasial #
moran = moran.test(tingkat_miskin$Tingkat_Miskin, W); moran

# Ward's Method Spatial Clustering #
library(ClustGeo)

# Matriks jarak Euclidean (D0) dan matriks pembobot spasial (D1)
D0 <- dist(tingkat_miskin$Tingkat_Miskin)
D1 <- as.dist(1 - WC)

# Penentuan Alpha 
range.alpha <- seq(0, 1, by=0.1)
cr <- choicealpha(D0, D1, range.alpha, K=2, graph = TRUE, scale = TRUE)  
round(cr$Qnorm, 4)

# Dendrogram 
alpha = 0.3
spatialclust <- hclustgeo(D0, D1, alpha = alpha, scale = TRUE)  
plot(spatialclust,labels = tingkat_miskin$Wilayah,
     main = "Dendrogram Spatial Clustering",
     xlab = "", sub = "")

# Nilai Silhouette 
library(cluster)
max_k <- 10
silhouette_avg_k <- numeric(max_k - 1)

for (k in 2:max_k) {
  # Memotong dendrogram untuk mendapatkan k klaster
  clusters <- cutree(spatialclust, k = k)
  
  # Menghitung silhouette (menggunakan D0 sebagai matriks jarak)
  sil <- silhouette(clusters, dist = D0) 
  
  # Menyimpan rata-rata silhouette coefficient
  silhouette_avg_k[k - 1] <- mean(sil[, 3])
}

# Membuat DataFrame hasil (untuk plotting)
silhouette_results <- data.frame(
  Jumlah_Klaster = 2:max_k,
  Average_Silhouette_Width = silhouette_avg_k)

optimal_k <- silhouette_results$Jumlah_Klaster[which.max(silhouette_results$Average_Silhouette_Width)]
max_silhouette <- max(silhouette_results$Average_Silhouette_Width)

# Plot Rata-rata Silhouette Coefficient
plot(silhouette_results$Jumlah_Klaster, 
     silhouette_results$Average_Silhouette_Width, 
     type = "b", 
     pch = 19,
     col = "black",
     xlab = "Jumlah Klaster", 
     ylab = "Average Silhouette Width", 
     main = "Graph Silhouette Clustering Dengan Efek Spasial")
abline(v = optimal_k, col = "red", lty = 2)

optimal_k; max_silhouette
cluster_membership <- cutree(spatialclust, k = optimal_k)
jabar$cluster <- factor(cluster_membership)

# Dendrogram Final + Peta Spatial Clustering 
plot(spatialclust, 
     labels = tingkat_miskin$Wilayah,
     main = "Dendrogram Spatial Clustering",
     xlab = "", sub = "")
rect.hclust(spatialclust, k = optimal_k, border = c("red", "blue"))

library(tmap)
tmap_options(check.and.fix = TRUE)
tm_shape(jabar) +
  tm_fill("cluster", style = "cat", palette = c("cyan", "pink"), title = "KLASTER") +
  tm_borders() +
  tm_layout(frame = FALSE,            
            legend.position = c(1, 0.3),
            legend.title.fontface = "bold",
            legend.text.size = 1,   
            legend.title.size = 1)

# Rata-Rata Tingkat Kemiskinan di Tiap Klaster 
tingkat_miskin$Klaster <- cluster_membership

cluster_summary <- aggregate(tingkat_miskin$Tingkat_Miskin ~ Klaster, 
                             data = tingkat_miskin, 
                             FUN = mean)
names(cluster_summary)[2] <- "Rata_rata_Tingkat_Kemiskinan"
round(cluster_summary, 2)
