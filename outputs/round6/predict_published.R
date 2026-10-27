library(Matrix)
args<-commandArgs(trailingOnly=TRUE);ct<-args[[1]];source<-args[[2]];selection<-args[[3]];out<-args[[4]];work<-args[[5]]
stopifnot(file.exists(file.path(out,'all_count_inputs_ready.json')))
x<-readRDS(source);m<-read.csv(gzfile(selection),check.names=FALSE,stringsAsFactors=FALSE);ix<-match(m$barcode,colnames(x));stopifnot(!anyNA(ix));v<-read.csv(file.path(out,paste0('source_',ct,'_gene_annotation.csv')));stopifnot(identical(as.character(v$feature_name),rownames(x)))
co<-read.csv(file.path(out,paste0('published_',ct,'_weights.csv')));tr<-read.csv(file.path(out,paste0('published_',ct,'_source_mapping_triplets.csv')));b<-scan(file.path(out,paste0('published_',ct,'_intercept.txt')),quiet=TRUE);used<-sort(unique(tr$feature_index));M<-sparseMatrix(i=tr$model_index,j=match(tr$feature_index,used),x=1,dims=c(nrow(co),length(used)));pred<-numeric(length(ix))
for(lo in seq(1,length(ix),by=5000)) {
 jj<-lo:min(lo+4999,length(ix));a<-x[used,ix[jj],drop=FALSE];a@x<-expm1(a@x);y<-M%*%a;y@x<-log1p(y@x);pred[jj]<-as.numeric(crossprod(co$coefficient,y))+b
 if(lo==1) {
  t<-as.data.frame(as.matrix(y[,seq_len(min(3,ncol(y))),drop=FALSE]));names(t)<-paste0('cell_',seq_len(ncol(t)));t$gene<-co$gene;t$coefficient<-co$coefficient;write.csv(t,file.path(out,paste0('native_prediction_check_',ct,'.csv')),row.names=FALSE);write.csv(data.frame(prediction=pred[seq_len(min(3,length(pred)))]),file.path(out,paste0('native_prediction_check_',ct,'_scores.csv')),row.names=FALSE)
 }
 if(lo%%50000==1)cat('PREDICT',ct,lo,'/',length(ix),'\n')
}
stopifnot(all(is.finite(pred)));z<-data.frame(barcode=m$barcode,score=pred);target<-file.path(work,paste0('published_',ct,'_per_cell_scores.csv.gz'));part<-paste0(target,'.part');con<-gzfile(part,'wt');write.csv(z,con,row.names=FALSE);close(con);stopifnot(file.rename(part,target));cat('PREDICTION_COMPLETE',ct,length(pred),'\n')
