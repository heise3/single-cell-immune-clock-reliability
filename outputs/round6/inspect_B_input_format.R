library(Matrix)
args <- commandArgs(trailingOnly=TRUE)
base <- args[[1]]
out <- args[[2]]
x <- readRDS(file.path(base, 'b_cells_rna.rds'))
m <- read.csv(file.path(base,'b_cells_metadata.csv'),check.names=FALSE)
rownames(m) <- m[[1]]
stopifnot(all(colnames(x) %in% rownames(m)))
m <- m[colnames(x),,drop=FALSE]
stopifnot(inherits(x,'dgCMatrix'))
nn <- length(x@x)
z <- expm1(x@x)*rep(m$nCount_RNA,diff(x@p))/10000
err <- abs(z-round(z))
rawsum <- as.numeric(rowsum(z,rep(seq_len(ncol(x)),diff(x@p)),reorder=FALSE))
q <- data.frame(property=c('matrix_class','n_genes','n_cells','nonzeros','min_value','max_value','fraction_integer_as_released','hypothesis_log1p_CP10K_with_metadata_library','reconstructed_fraction_integer_tol1e-6','reconstructed_max_rounding_error','reconstructed_median_rounding_error','max_absolute_library_sum_error','median_absolute_library_sum_error','max_abs_feature_count_vs_metadata_error','negative_reconstructed_values','feature_id_examples'),value=c(class(x)[1],nrow(x),ncol(x),nn,min(x@x),max(x@x),mean(abs(x@x-round(x@x))<1e-8),'format_validation_only_not_clock_prediction',mean(err<1e-6),max(err),median(err),max(abs(rawsum-m$nCount_RNA)),median(abs(rawsum-m$nCount_RNA)),max(abs(diff(x@p)-m$nFeature_RNA)),sum(z<0),paste(head(rownames(x)),collapse='|')))
write.csv(q,file.path(out,'K_a_B_matrix_format_QA.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_input_inspection_sessionInfo.txt'))
print(q,row.names=FALSE)
