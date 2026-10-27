# Verify the Python affine-score implementation against the original glmnet
# predictor on normalized real cells, using the exact published gene order.
library(glmnet)
repo <- "work/models/scImmuAging-112f8fab10c72aec6e7a7d83762a40117a1e6e24"
models <- readRDS(file.path(repo,"data/all_model.RDS"))
features <- readRDS(file.path(repo,"data/all_model_inputfeatures.RDS"))
input <- read.csv("outputs/K4_linear_validation_inputs.csv",check.names=FALSE)
out <- list()
for (ct in names(models)) {
  rows <- input[input$celltype==ct,,drop=FALSE]
  genes <- features[[ct]]
  x <- matrix(0,nrow=nrow(rows),ncol=length(genes),dimnames=list(NULL,genes))
  common <- intersect(genes,colnames(rows))
  x[,common] <- as.matrix(rows[,common,drop=FALSE])
  x[is.na(x)] <- 0
  native <- as.numeric(predict(models[[ct]],newx=x,s="lambda.min"))
  out[[ct]] <- data.frame(celltype=ct,python=rows$python_prediction,native=native,abs_error=abs(native-rows$python_prediction))
}
out <- do.call(rbind,out)
write.csv(out,"outputs/K4_native_glmnet_validation.csv",row.names=FALSE)
stopifnot(max(out$abs_error)<1e-8)
cat("Native glmnet equivalence PASS; maximum absolute error",max(out$abs_error),"years\n")
