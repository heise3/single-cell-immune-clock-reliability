#!/usr/bin/env python3
"""Recover completed CSR X after AnnData failed while writing obs index metadata."""
from pathlib import Path
import json,hashlib,re
import numpy as np,pandas as pd,h5py,anndata as ad
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';W=P/'work/round5'
def digest(x):
 h=hashlib.sha256()
 for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
 return h.hexdigest()
def main():
 target=W/'P4_additional_measured_counts.h5ad';receipt=O/'P4_count_write_recovery_QA.json'
 if receipt.exists():return
 log=(W/'P4_analysis_rowquery.log').read_text();assert 'P4_EXTRA_BLOCK 36805' in log and "DataFrame.index.name ('soma_joinid')" in log
 obs=pd.read_parquet(W/'P4_selected_cells.parquet').sort_values('soma_joinid');m=pd.read_csv(O/'P4_AIDA_published_gene_mapping.csv').fillna('').set_index('model_ensembl');co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');f=ad.read_h5ad(P/'work/round3/AIDA_evaluation_counts/AIDA_000.h5ad',backed='r');cached=set(f.var.feature_id);f.file.close();needed=sorted(set(v for e in co[co.gene.ne('(Intercept)')&co.coefficient.ne(0)].gene for v in m.loc[e,'mapped_feature_ids'].split('|') if v)-cached)
 with h5py.File(target,'r') as h:x=ad.io.read_elem(h['X'])
 assert x.shape==(len(obs),len(needed))==(36805,1145) and x.nnz==7040984 and np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all();before=digest(x);oi=pd.DataFrame({'soma_joinid':obs.soma_joinid.to_numpy()},index=pd.Index(obs.soma_joinid.astype(str),name=None));oi.index.name=None;vi=pd.DataFrame({'feature_id':needed},index=pd.Index(needed,name=None));a=ad.AnnData(x,obs=oi,var=vi);new=target.with_suffix('.recovered.h5ad');a.write_h5ad(new,compression='lzf');b=ad.read_h5ad(new);assert digest(b.X)==before and list(b.obs.soma_joinid)==list(obs.soma_joinid) and list(b.var.feature_id)==needed;partial=target.with_suffix('.partial_WRITE_ERROR.h5ad');target.rename(partial);new.rename(target);receipt.write_text(json.dumps({'status':'PASS','failure':'AnnData obs index name collided with integer soma_joinid column after completed X write','complete_input_query_log':'work/round5/P4_analysis_rowquery.log;all36805cells read before obs write error','recovery':'Read fully written CSR X;reconstruct obs/var labels from same frozen selected-cell sort and deterministic gene mapping; no clock outcome or data imputation','shape':[36805,1145],'nnz':int(x.nnz),'CSR_buffer_sha256_before_and_after':before,'integer_nonnegative_finite_verified':True,'partial_file_preserved':str(partial.relative_to(P)),'no_new_network_query':True},indent=2))
if __name__=='__main__':main()
