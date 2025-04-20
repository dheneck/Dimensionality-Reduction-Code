# Clear
rm(list = ls())

# Remove Scientific Notation 
options(scipen = 999)

# Set Working Directory 
setwd("C:\\Users\\dhene\\OneDrive\\Desktop\\Dimens Reduction")

# Packages
library(FactoMineR)
library(factoextra)
library(psych)
library(dplyr)
library(ggplot2)
library(reshape2)
library(sem)
library(stargazer)
library(lavaan)
library(knitr)

ess_data <- read.csv("C:\\Users\\dhene\\OneDrive\\Desktop\\Dimens Reduction\\ESS9e03_2.csv")

# Example variables from your data
trust_vars <- c("trstprl", "trstlgl", "trstplc", 
                "trstplt", "trstprt", "trstep", "trstun")
# Subset dataset
trust_data <- ess_data[trust_vars]

# High-trust: Scandinavian
scandi <- c("SE", "NO", "DK", "FI")

# Mediterranean 
medi <- c("ES", "PT", "IT", "CY", "FR")

# Pearson correlation
cor_pearson <- cor(trust_data, use = "pairwise.complete.obs", method = "pearson")

# Polychoric correlation 
poly_cor <- lavCor(trust_data, ordered = trust_vars)
pearson_cor <- cor(trust_data, use = "pairwise.complete.obs")

kable(round(pearson_cor, 2), caption = "Pearson Correlation Matrix")
kable(round(poly_cor, 2), caption = "Polychoric Correlation Matrix")

# Create dataframes for both subsets 
scandi_data <- subset(ess_data, cntry %in% scandi)
scandi_trust_data <- na.omit(scandi_data[trust_vars])

medi_data <- subset(ess_data, cntry %in% medi)
medi_trust_data <- na.omit(medi_data[trust_vars])

data_scaled_scandi <- scale(scandi_trust_data) 
data_scaled_medi <- scale(medi_trust_data)  
data_scaled <- scale(trust_data)

# Run PCA on Scandinavian data
data_scaled <- scale(scandi_trust_data)  # Standardize for PCA

pca_result <- prcomp(data_scaled_scandi, center = TRUE, scale. = TRUE)
res.pca <- PCA(data_scaled, graph = FALSE)

# Extract loadings, contributions and cos2 for the first 3 PCs
loadings <- res.pca$var$coord[, 1:2]
contribs <- res.pca$var$contrib[, 1:2]
cos2 <- res.pca$var$cos2[, 1:2]

# Round and combine into a table
pca_summary <- data.frame(
  Loading_PC1 = round(loadings[,1], 3),
  Contrib_PC1 = round(contribs[,1], 1),
  Cos2_PC1 = round(cos2[,1], 3),
  
  Loading_PC2 = round(loadings[,2], 3),
  Contrib_PC2 = round(contribs[,2], 1),
  Cos2_PC2 = round(cos2[,2], 3)
)
rownames(pca_summary) <- c("Trust in parliament", "Trust in the legal system", 
                           "Trust in the police", "Trust in politicians", 
                           "Trust in political parties", "Trust in the European Parliament", 
                           "Trust in the United Nations")
kable(pca_summary, caption = "Principal Component Analysis Results for EU Countries",
      col.names = c("Loading PC1", "Contrib PC1", "Cos² PC1", 
                    "Loading PC2", "Contrib PC2", "Cos² PC2"), align = "c")

res.pca <- PCA(data_scaled_scandi, graph = FALSE)

# Extract loadings, contributions and cos2 for the first 3 PCs
loadings <- res.pca$var$coord[, 1:2]
contribs <- res.pca$var$contrib[, 1:2]
cos2 <- res.pca$var$cos2[, 1:2]

# Round and combine into a table
pca_summary <- data.frame(
  Loading_PC1 = round(loadings[,1], 3),
  Contrib_PC1 = round(contribs[,1], 1),
  Cos2_PC1 = round(cos2[,1], 3),
  
  Loading_PC2 = round(loadings[,2], 3),
  Contrib_PC2 = round(contribs[,2], 1),
  Cos2_PC2 = round(cos2[,2], 3)
)
rownames(pca_summary) <- c("Trust in parliament", "Trust in the legal system", 
                           "Trust in the police", "Trust in politicians", 
                           "Trust in political parties", "Trust in the European Parliament", 
                           "Trust in the United Nations")
kable(pca_summary, caption = "Principal Component Analysis Results for Scandinavian Countries",
      col.names = c("Loading PC1", "Contrib PC1", "Cos² PC1", 
                    "Loading PC2", "Contrib PC2", "Cos² PC2"), align = "c")

res.pca <- PCA(data_scaled_medi, graph = FALSE)

# Extract loadings, contributions and cos2 for the first 2 PCs
loadings <- res.pca$var$coord[, 1:2]
contribs <- res.pca$var$contrib[, 1:2]
cos2 <- res.pca$var$cos2[, 1:2]

# Round and combine into a table
pca_summary <- data.frame(
  Loading_PC1 = round(loadings[,1], 3),
  Contrib_PC1 = round(contribs[,1], 1),
  Cos2_PC1 = round(cos2[,1], 3),
  
  Loading_PC2 = round(loadings[,2], 3),
  Contrib_PC2 = round(contribs[,2], 1),
  Cos2_PC2 = round(cos2[,2], 3)
)
rownames(pca_summary) <- c("Trust in parliament", "Trust in the legal system", 
                           "Trust in the police", "Trust in politicians", 
                           "Trust in political parties", "Trust in the European Parliament", 
                           "Trust in the United Nations")

kable(pca_summary, caption = "Principal Component Analysis Results for Mediterranean Countries",
      col.names = c("Loading PC1", "Contrib PC1", "Cos² PC1", 
                    "Loading PC2", "Contrib PC2", "Cos² PC2"), align = "c")

# PCA for high-trust countries
pca_scandi <- PCA(data_scaled_scandi, scale.unit = TRUE, ncp = 2, graph = FALSE)

# PCA for med countries
pca_medi <- PCA(medi_trust_data, scale.unit = TRUE, ncp = 2, graph = FALSE)

# PCA for EU
pca_eu <- PCA(trust_data, scale.unit = TRUE, ncp = 2, graph = FALSE)

fviz_pca_var(pca_eu, 
             repel = TRUE,
             label = "var",
             title = "EU (Full Sample)",
             xlab = "Dimension 1",
             ylab = "Dimension 2",
             col.var = "darkgreen", col.ind = "gray40")


fviz_pca_var(pca_scandi, 
             repel = TRUE,
             label = "var",
             xlab = "Dimension 1",
             ylab = "Dimension 2",
             title = "Scandinavian Countries",
             col.var = "darkgreen", col.ind = "gray40")
fviz_pca_var(pca_medi, 
             repel = TRUE,
             label = "var",
             xlab = "Dimension 1",
             ylab = "Dimension 2",
             title = "Mediterranean Countries",
             col.var = "darkgreen", col.ind = "gray40")

pca_result <- prcomp(trust_data, scale. = TRUE)
# Extract the PCA scores
pca_scores <- as.data.frame(pca_result$x[, 1:2])  # Only PC1 and PC2

# Rename columns for clarity
colnames(pca_scores) <- c("PC1", "PC2")

# Add the scores to your original data
ess_data_complete <- cbind(trust_data, pca_scores)


# Cleaning Data 
regression_vars <- c(
  "lrscale",      # Left-right ideology
  "agea",         # Age
  "eduyrs",       # Years of education
  "gndr",         # Gender
  "stfdem",       # Satisfaction with democracy
  "stflife",      # Satisfaction with life
  "PC1",          # First principal component (outcome)
  "PC2"           # Second principal component (outcome)
)

# Add additional regression variables from original data
extra_vars <- c("lrscale", "agea", "eduyrs", "gndr", "stfdem", "stflife", "polintr",
                "ppltrst", "pplfair", "pplhlp")

# Subset rows from ess_data that match trust_data
extra_data <- ess_data[rownames(trust_data), extra_vars]

# Combine with PC scores
regression_data <- cbind(extra_data, pca_scores)

# Remove rows with any NAs
regression_data <- na.omit(regression_data)

# Clean Gender
regression_data$gndr <- factor(regression_data$gndr, 
                               levels = c(1, 2), 
                               labels = c("Male", "Female"))

# Clean L-R Scale
regression_data <- regression_data[regression_data$lrscale >= 0 & regression_data$lrscale <= 10, ]

# Clean Education 
regression_data <- regression_data[regression_data$eduyrs < 90, ]

# Overall 
pca_all <- PCA(trust_data, scale.unit = TRUE, ncp = 5, graph = FALSE)
ess_data_clean <- ess_data[rownames(trust_data), ]  # Match full data rows

pca_scores_all <- get_pca_ind(pca_all)$coord
ess_data_clean$PC1 <- pca_scores_all[, 1]
ess_data_clean$PC2 <- pca_scores_all[, 2]

model_pc1 <- lm(PC1 ~ lrscale + agea + eduyrs + gndr + stfdem 
                + stflife + ppltrst + polintr, data = regression_data)

model_pc2 <- lm(PC2 ~ lrscale + agea + eduyrs + gndr + stfdem 
                + stflife + ppltrst + polintr, data = regression_data)

cat("\\renewcommand{\\arraystretch}{0.7}\n")
stargazer(model_pc1, model_pc2,
          type = "latex",
          title = "Regression Results: PC1 and PC2",
          column.labels = c("PC1", "PC2"),      
          dep.var.labels.include = FALSE,       
          model.names = FALSE,
          header = FALSE,
          covariate.labels = c("Left-Right Scale", "Age", "Education (Years)", 
                               "Female", "Sat. with Democracy", "Sat. with Life",
                               "General trust in People", "Political Interest"
          ),
          digits = 3
)

# Run PCA on eu data
data_scaled <- scale(trust_data)  # Standardize for PCA

pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)


# Eigenvalues (squared standard deviations)
eigenvalues <- pca_result$sdev^2

# Proportion of variance explained
variance <- eigenvalues / sum(eigenvalues) * 100

# Cumulative variance
cumvar <- cumsum(variance)

# Combine into a data frame
pca_table <- data.frame(
  Eigenvalue = round(eigenvalues, 3),
  Variance = round(variance, 1),
  `Cumulative Variance` = round(cumvar, 1)
)

# Optional: Rename rows as PC1, PC2, etc.
rownames(pca_table) <- paste0("PC", 1:nrow(pca_table))

# View the table
kable(pca_table, caption = "Principal Component Eigenvalues and Varaince (EU)")

# Eigenvalues (squared standard deviations)
eigenvalues <- pca_result$sdev^2

# Proportion of variance explained
variance <- eigenvalues / sum(eigenvalues) * 100

# Cumulative variance
cumvar <- cumsum(variance)

# Combine into a data frame
pca_table <- data.frame(
  Eigenvalue = round(eigenvalues, 3),
  Variance = round(variance, 1),
  `Cumulative Variance` = round(cumvar, 1)
)

# Optional: Rename rows as PC1, PC2, etc.
rownames(pca_table) <- paste0("PC", 1:nrow(pca_table))

# View the table
kable(pca_table, caption = "Principal Component Eigenvalues and Varaince (Scandinavia)")

pca_result <- prcomp(data_scaled_medi, center = TRUE, scale. = TRUE)


# Eigenvalues (squared standard deviations)
eigenvalues <- pca_result$sdev^2

# Proportion of variance explained
variance <- eigenvalues / sum(eigenvalues) * 100

# Cumulative variance
cumvar <- cumsum(variance)

# Combine into a data frame
pca_table <- data.frame(
  Eigenvalue = round(eigenvalues, 3),
  Variance = round(variance, 1),
  `Cumulative Variance` = round(cumvar, 1)
)

# Optional: Rename rows as PC1, PC2, etc.
rownames(pca_table) <- paste0("PC", 1:nrow(pca_table))

# View the table
kable(pca_table, caption = "Principal Component Eigenvalues and Variance (Mediterranean)")