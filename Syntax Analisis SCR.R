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

# Analisis Regresi Linier Multigrup #
datafinal <- read_excel("Data Regresi.xlsx")
head(datafinal)

library(lmtest)
model_ols <- lm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + D + X1D + X2D
                + X3D + X4D + X5D + X6D + X7D + X8D, data = datafinal)
summary(model_ols)

# Uji Asumsi Klasik #
# Uji Normalitas
ks.test(residuals(model_ols), "pnorm",
        mean = mean(residuals(model_ols)),
        sd   = sd(residuals(model_ols)))

# Uji Homoskedastisitas
bptest(model_ols) 

# Uji Non-Autokorelasi
library(car)
durbinWatsonTest(model_ols)

# Uji Multikolinearitas
vif(model_ols) #Tidak Terpenuhi (VIF > 10)

# Ridge Regression # 
predictors <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", 
                "D", "X1D", "X2D", "X3D", "X4D", "X5D", "X6D", "X7D", "X8D")
X_matrix <- as.matrix(datafinal[, predictors])
Y_vector <- datafinal$Y

penalty <- c(rep(1, 8),  # X1–X8
             0,          # D (TIDAK dipenalti)
             rep(1, 8))  # X1D–X8D

library(glmnet)
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=3); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=4); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=5); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=6); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=7); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=8);cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=9); cv_ridge
set.seed(123)
cv_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, penalty.factor = penalty, nfolds=10); cv_ridge
best_lambda <- 0.6392
ridge_model <- glmnet(X_matrix, Y_vector, alpha = 0, lambda = best_lambda, penalty.factor = penalty)
coef(ridge_model)

# Identifikasi Pencilan #
cooks_dist <- cooks.distance(model_ols)

n <- nrow(datafinal)               # Jumlah observasi
p <- length(coef(model_ols))       # Jumlah parameter (termasuk intercept)
alpha <- 0.5

# Metode F-distribution (α=0.5)
df1 <- p + 1                       # Derajat kebebasan 1 = p + 1
df2 <- n - (p + 1)                 # Derajat kebebasan 2 = n - (p + 1)
cutoff_f <- qf(alpha, df1 = df1, df2 = df2)

# Metode Fixed threshold 1.0 
cutoff_fixed <- 1.0

# Identifikasi outliers untuk setiap metode
outliers_f <- which(cooks_dist > cutoff_f)
outliers_fixed <- which(cooks_dist > cutoff_fixed)

# Gabungkan semua outliers (unique)
all_outliers <- unique(c(outliers_f, outliers_fixed))

cat("Model Information:\n")
cat("Observations (n):", n, "\n")
cat("Parameters (p):", p, "\n")
cat("Degrees of freedom (df1 = p+1):", df1, "\n")
cat("Degrees of freedom (df2 = n-(p+1)):", df2, "\n\n")

cat("Cut-off Values:\n")
cat("F-distribution (α=0.5):", round(cutoff_f, 4), 
    " (F(", df1, ",", df2, "))\n")
cat("Fixed threshold:", cutoff_fixed, "\n\n")

cat("Outlier Detection Results:\n")
cat("Method F-distribution:", length(outliers_f), "outliers detected\n")
cat("Method Fixed 1.0     :", length(outliers_fixed), "outliers detected\n")
cat("Total unique outliers    :", length(all_outliers), "outliers\n\n")

if (length(all_outliers) > 0) {
  # Buat dataframe untuk outliers
  outliers_df <- datafinal[all_outliers, ]
  outliers_df$Observation <- all_outliers
  outliers_df$CooksD <- cooks_dist[all_outliers]
  outliers_df$Exceeds_F <- outliers_df$CooksD > cutoff_f
  outliers_df$Exceeds_Fixed <- outliers_df$CooksD > cutoff_fixed
  outliers_df$Methods_Detected <- rowSums(outliers_df[, c("Exceeds_F", "Exceeds_Fixed")])
  
  # Pilih kolom untuk ditampilkan
  cols_to_show <- c("Observation", "Y", predictors[1:min(4, length(predictors))], 
                    "CooksD", "Methods_Detected", "Exceeds_F", "Exceeds_Fixed")
  cols_to_show <- cols_to_show[cols_to_show %in% names(outliers_df)]
  
  # Ambil subset dataframe
  display_df <- outliers_df[, cols_to_show]
  
  # Format nama kolom
  names(display_df)[names(display_df) == "Exceeds_F"] <- "> F(α=0.5)"
  names(display_df)[names(display_df) == "Exceeds_Fixed"] <- "> 1.0"
  names(display_df)[names(display_df) == "Methods_Detected"] <- "# Methods"
  
  cat("=========================================================\n")
  cat("DETAILED OUTLIER INFORMATION\n")
  cat("=========================================================\n\n")
  
  # Tampilkan dalam format tabel sederhana
  print(display_df, digits = 4)
  
  # Tampilkan observasi dengan Cook's Distance tertinggi
  max_cooks <- which.max(cooks_dist)
  cat("\n=========================================================\n")
  cat("OBSERVASI DENGAN COOK'S DISTANCE TERTINGGI\n")
  cat("=========================================================\n")
  cat("Observation:", max_cooks, "\n")
  cat("Cook's Distance:", round(cooks_dist[max_cooks], 4), "\n")
  cat("Y value:", datafinal$Y[max_cooks], "\n")
  cat("Exceeds all thresholds:", 
      ifelse(cooks_dist[max_cooks] > cutoff_f &
               cooks_dist[max_cooks] > cutoff_fixed, "YES", "NO"), "\n")
} else {
  cat("No outliers detected by any method.\n")
}

# Robust Regression #
library(MASS)
robust_model <- rlm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + D + 
                      X1D + X2D + X3D + X4D + X5D + X6D + X7D + X8D, 
                    data = datafinal, method = "M", maxit = 110)
summary(robust_model)

# Tampilkan bobot dari robust regression
weights <- robust_model$w
min(weights)
max(weights)
mean(weights)
cat("Observations with weight < 0.5:", sum(weights < 0.5), "\n")
cat("Observations with weight < 0.3:", sum(weights < 0.3), "\n")

# Identifikasi observasi dengan bobot rendah (potensi outliers)
low_weight_obs <- which(weights < 0.5)
cat("\nObservations with low weights (< 0.5):", toString(low_weight_obs), "\n")

# Ridge Robust Regression #
# Dapatkan bobot dari robust regression
weights <- robust_model$w

# Cross-validation untuk Ridge Robust dengan Bobbot
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 3); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 4); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 5); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 6); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 7); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 8); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 9); cv_robust_ridge
set.seed(123)
cv_robust_ridge <- cv.glmnet(X_matrix, Y_vector, alpha = 0, weights = weights, penalty.factor = penalty, nfolds = 10); cv_robust_ridge

lambda_robust_min <- 0.508

# Estimasi Ridge Robust Regression dengan lambda optimal
robust_ridge_model <- glmnet(X_matrix, Y_vector, 
                             alpha = 0, 
                             lambda = lambda_robust_min,
                             weights = weights,
                             penalty.factor = penalty)

robust_ridge_coef <- as.matrix(coef(robust_ridge_model)); robust_ridge_coef

# Pemilihan Model Terbaik #
calculate_metrics <- function(y_true, y_pred, n_p = NULL) {
  # Hitung R²
  SSE <- sum((y_true - y_pred)^2)  # Sum of Squared Errors
  SST <- sum((y_true - mean(y_true))^2)  # Total Sum of Squares
  R2 <- 1 - (SSE / SST)
  
  # Hitung Adjusted R² jika n_p diberikan
  if (!is.null(n_p)) {
    n <- length(y_true)
    adj_R2 <- 1 - (1 - R2) * ((n - 1) / (n - n_p - 1))
  } else {
    adj_R2 <- NA
  }
  
  # Hitung RMSE
  RMSE <- sqrt(mean((y_true - y_pred)^2))
  
  # Hitung MAE
  MAE <- mean(abs(y_true - y_pred))
  
  return(list(R2 = R2, adj_R2 = adj_R2, RMSE = RMSE, MAE = MAE, 
              SSE = SSE, SST = SST))
}

# Prediksi OLS
y_pred_ols <- predict(model_ols)

# Prediksi Ridge
y_pred_ridge <- predict(ridge_model, newx = X_matrix)[,1]

# Prediksi Robust Ridge
y_pred_robust_ridge <- predict(robust_ridge_model, newx = X_matrix)[,1]

# Hitung Metrik untuk Setiap Model
# OLS Metrics
n_p_ols <- length(coef(model_ols))  # jumlah parameter
metrics_ols <- calculate_metrics(Y_vector, y_pred_ols, n_p_ols)

# Ridge Metrics
n_p_ridge <- sum(coef(ridge_model) != 0)  # jumlah parameter non-zero
metrics_ridge <- calculate_metrics(Y_vector, y_pred_ridge)

# Robust Ridge Metrics
n_p_robust_ridge <- sum(coef(robust_ridge_model) != 0)
metrics_robust_ridge <- calculate_metrics(Y_vector, y_pred_robust_ridge)

comparison_table <- data.frame(
  Model = c("OLS", "Ridge", "Robust Ridge"),
  R2 = c(round(metrics_ols$R2, 4), 
         round(metrics_ridge$R2, 4), 
         round(metrics_robust_ridge$R2, 4)),
  Adj_R2 = c(round(metrics_ols$adj_R2, 4), NA, NA),
  RMSE = c(round(metrics_ols$RMSE, 4), 
           round(metrics_ridge$RMSE, 4), 
           round(metrics_robust_ridge$RMSE, 4)),
  MAE = c(round(metrics_ols$MAE, 4), 
          round(metrics_ridge$MAE, 4), 
          round(metrics_robust_ridge$MAE, 4))); comparison_table