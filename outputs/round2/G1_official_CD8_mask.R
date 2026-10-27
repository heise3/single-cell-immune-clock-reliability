# Same pseudocell membership for complete-available and residual-gene-masked data.
.libPaths(c("work/Rlib", .libPaths()))
suppressPackageStartupMessages({library(SeuratObject); library(dplyr); library(tidyr); library(purrr); library(glmnet); library(jsonlite)})
set.seed(20261025)
repo <- "work/models/scImmuAging-112f8fab10c72aec6e7a7d83762a40117a1e6e24"
source(file.path(repo, "R/main.R"))
models <- readRDS(file.path(repo, "data/all_model.RDS"))
features <- readRDS(file.path(repo, "data/all_model_inputfeatures.RDS"))
input <- readRDS(file.path(repo, "data/demo_bcg_cd8t.RDS"))
mapping <- read.csv("outputs/round2/scImmuAging_gene_mapping.csv", na.strings=NULL)
co <- read.csv("outputs/scImmuAging_coefficients.csv")
required <- co$gene[co$celltype=="CD8T" & co$coefficient!=0 & co$gene!="(Intercept)"]
missing <- intersect(required, mapping$model_ensembl[mapping$mapped_geo_ids==""])
pre <- PreProcess(input, "CD8T", models$CD8T, features$CD8T)
full <- AgingClockCalculator(pre, models$CD8T, features$CD8T)
masked_pre <- pre
present <- intersect(missing, colnames(pre))
masked_pre[,present] <- 0
masked <- AgingClockCalculator(masked_pre, models$CD8T, features$CD8T)
a <- aggregate(Prediction~donor_id+age, full, mean)
b <- aggregate(Prediction~donor_id+age, masked, mean)
names(a)[3] <- "full"; names(b)[3] <- "masked"
out <- merge(a,b,by=c("donor_id","age"))
write.csv(out,"outputs/round2/G1_official_CD8_demo_mask_predictions.csv",row.names=FALSE)
r <- cor(out$full,out$masked)
summary <- list(seed=20261025,clock="scImmuAging",celltype="CD8T",validation_source="official demo_bcg_cd8t.RDS",
 n_donors=nrow(out),n_cells=ncol(input),Pearson_r=r,bias_years=mean(out$masked-out$full),
 full_MAE=mean(abs(out$full-out$age)),masked_MAE=mean(abs(out$masked-out$age)),
 MAE_change=mean(abs(out$masked-out$age))-mean(abs(out$full-out$age)),
 residual_mask_genes=missing,mask_genes_observed_in_demo=present,
 baseline_nonzero_model_features_in_demo=length(intersect(required,rownames(input))),baseline_nonzero_model_features=length(required),
 status=ifelse(r<0.9,"EXCLUDE_FROM_MAIN","PASS"),
 limitation="Five official-example donors, CD8 only; native full prediction already zero-fills 39 unavailable nonzero-weight genes. This is a perturbation check, not an external accuracy validation.")
write_json(summary,"outputs/round2/G1_official_CD8_demo_mask_summary.json",pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),"outputs/round2/G1_mask_R_sessionInfo.txt")
print(summary)
