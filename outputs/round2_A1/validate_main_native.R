suppressPackageStartupMessages(library(glmnet))
set.seed(20261025)
repo <- 'work/models/scImmuAging-112f8fab10c72aec6e7a7d83762a40117a1e6e24'
models <- readRDS(file.path(repo,'data/all_model.RDS'))
features <- readRDS(file.path(repo,'data/all_model_inputfeatures.RDS'))
input <- read.csv('outputs/round2_A1/main_native_prediction_input.csv',check.names=FALSE)
out <- list()
for (ct in names(models)) {
 z <- input[input$celltype==ct,,drop=FALSE]
 x <- matrix(0,nrow=nrow(z),ncol=length(features[[ct]]),dimnames=list(NULL,features[[ct]]))
 common <- intersect(colnames(z),colnames(x));x[,common] <- as.matrix(z[,common,drop=FALSE]);x[is.na(x)] <- 0
 pred <- as.numeric(predict(models[[ct]],newx=x,s='lambda.min'))
 out[[ct]] <- data.frame(celltype=ct,python=z$python_prediction,native=pred,absolute_error=abs(pred-z$python_prediction))
}
out <- do.call(rbind,out)
write.csv(out,'outputs/round2_A1/main_native_glmnet_validation.csv',row.names=FALSE)
stopifnot(max(out$absolute_error)<1e-8)
writeLines(capture.output(sessionInfo()),'outputs/round2_A1/main_R_sessionInfo.txt')
cat('PASS max absolute error:',max(out$absolute_error),'years\n')
