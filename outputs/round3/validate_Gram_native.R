library(glmnet)
root <- "work/round3/"
x <- as.matrix(read.table(paste0(root,"native_Gram_X.tsv"),sep="\t",header=FALSE))
yw <- read.delim(paste0(root,"native_Gram_yw.tsv"))
lam <- read.delim(paste0(root,"native_Gram_lambda.tsv"))$lambda
fit <- glmnet(x,yw$y,weights=yw$w,alpha=1,standardize=FALSE,intercept=TRUE,lambda=sort(lam,decreasing=TRUE),thresh=1e-12,maxit=1000000)
pred <- predict(fit,newx=x,s=lam)
write.table(pred,paste0(root,"native_Gram_R.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
capture.output(sessionInfo(),file="outputs/round3/R2_native_R_sessionInfo.txt")
