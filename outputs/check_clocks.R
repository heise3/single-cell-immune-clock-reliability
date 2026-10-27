# Native public-model smoke tests. No study-level inference is made here.
.libPaths(c("work/Rlib", .libPaths()))
set.seed(20261025)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(purrr)
library(glmnet)
repo <- "work/models/scImmuAging-112f8fab10c72aec6e7a7d83762a40117a1e6e24"
source(file.path(repo,"R/main.R"))
model_set <- readRDS(file.path(repo,"data/all_model.RDS"))
feature_set <- readRDS(file.path(repo,"data/all_model_inputfeatures.RDS"))
input <- readRDS(file.path(repo,"data/demo_bcg_cd8t.RDS"))
cat("sc-ImmuAging demo class", class(input), "dim", dim(input), "\n")
pre <- PreProcess(input, "CD8T", model_set[["CD8T"]], feature_set[["CD8T"]])
pred <- AgingClockCalculator(pre, model_set[["CD8T"]], feature_set[["CD8T"]])
write.csv(Age_Donor(pred), "outputs/scImmuAging_native_demo.csv",row.names=FALSE)
cat("sc-ImmuAging native demo successful:", nrow(pred),"pseudocell predictions\n")
coeffs <- list()
for (celltype in names(model_set)) {
  co <- as.matrix(coef(model_set[[celltype]],s="lambda.min"))
  cat(celltype, "inputs", length(feature_set[[celltype]]),"nonzero",sum(co[,1]!=0),"\n")
  coeffs[[celltype]] <- data.frame(celltype=celltype, gene=rownames(co), coefficient=co[,1])
}
write.csv(bind_rows(coeffs),"outputs/scImmuAging_coefficients.csv",row.names=FALSE)
if (requireNamespace("randomForest",quietly=TRUE)) {
  library(randomForest)
  model <- readRDS("work/siAgeModel.rds")
  print(class(model));print(dim(model$importance))
  features <- rownames(model$importance)
  test <- as.data.frame(matrix(0,nrow=1,ncol=length(features)))
  colnames(test) <- features
  result <- predict(model,test)
  stopifnot(length(result)==1,is.finite(result))
  write.csv(data.frame(input="zero-vector software smoke test only",prediction=result),"outputs/siAge_software_smoke.csv",row.names=FALSE)
  writeLines(features,"outputs/siAge_required_features.txt")
  cat("siAge native model loaded and predicted finite output; annotation reference not retrieved\n")
}
writeLines(capture.output(sessionInfo()),"outputs/R_sessionInfo.txt")
