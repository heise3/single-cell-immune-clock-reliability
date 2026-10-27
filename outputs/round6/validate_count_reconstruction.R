library(Matrix)
args <- commandArgs(trailingOnly=TRUE);ct<-args[[1]];source<-args[[2]];selection<-args[[3]];out<-args[[4]]
x<-readRDS(source);stopifnot(inherits(x,'dgCMatrix'));m<-read.csv(gzfile(selection),check.names=FALSE,stringsAsFactors=FALSE);ix<-match(m$barcode,colnames(x));stopifnot(!anyNA(ix),!anyDuplicated(m$barcode));mxerr<-0;mxlib<-0;mxfeature<-0;nn<-0;rawdirectinteger<-0
for(lo in seq(1,length(ix),by=2000)) {
 jj<-lo:min(lo+1999,length(ix));a<-x[,ix[jj],drop=FALSE];stopifnot(all(is.finite(a@x)),all(a@x>=0));nn<-nn+length(a@x);rawdirectinteger<-rawdirectinteger+sum(abs(a@x-round(a@x))<1e-8);a@x<-expm1(a@x)*rep(m$nCount_RNA[jj],diff(a@p))/10000;err<-abs(a@x-round(a@x));mxerr<-max(mxerr,max(err));mxlib<-max(mxlib,max(abs(colSums(a)-m$nCount_RNA[jj])));mxfeature<-max(mxfeature,max(abs(diff(a@p)-m$nFeature_RNA[jj])));stopifnot(mxerr<1e-7,mxlib<1e-7,mxfeature==0)
 if(lo%%20000==1)cat('RAW_QA',ct,lo,'/',length(ix),'error',mxerr,'\n')
}
write.csv(data.frame(feature_index=seq_len(nrow(x)),feature_name=rownames(x)),file.path(out,paste0('source_',ct,'_gene_annotation.csv')),row.names=FALSE)
z<-data.frame(celltype=ct,source_genes=nrow(x),source_cells=ncol(x),analysis_cells=length(ix),analysis_nonzero_values=nn,released_fraction_integer=rawdirectinteger/nn,count_reconstruction='expm1(logCP10K)*released_nCount_RNA/10000',max_reconstruction_integer_error=mxerr,max_full_library_sum_error=mxlib,max_feature_count_vs_metadata_error=mxfeature,native_QC_verified=TRUE,clock_values_computed=FALSE)
write.csv(z,file.path(out,paste0('count_reconstruction_QA_',ct,'.csv')),row.names=FALSE);cat('COUNT_RECONSTRUCTION_COMPLETE',ct,'\n')
