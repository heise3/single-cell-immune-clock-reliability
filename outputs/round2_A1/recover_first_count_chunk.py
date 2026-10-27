#!/usr/bin/env python3
"""Recover only a complete X written before a known obs-index-name write error.
Never recover an interrupted or partial X. No expression values are changed.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
P=Path(__file__).resolve().parents[2];ap=argparse.ArgumentParser();ap.add_argument('--sensitivity87',action='store_true');a=ap.parse_args();O=P/('outputs/round2_A1/sensitivity87' if a.sensitivity87 else 'outputs/round2_A1');W=P/('work/round2_A1/sensitivity87' if a.sensitivity87 else 'work/round2_A1');p=W/'main_Census_counts/HIHA_K0_000.part.h5ad';dest=p.with_name('HIHA_K0_000.h5ad')
log=(W/'fetch_main_Census_counts.pipeline.log').read_text();assert "DataFrame.index.name ('soma_joinid')" in log,'Only the documented post-X obs-index error authorizes this recovery'
obs=pd.read_parquet(W/'cache/HIHA_K0_clock_cells.parquet').sort_values('soma_joinid').iloc[:100000].copy();v=pd.read_csv(O/'HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv');m=pd.read_csv(O/'gene_mapping_HIHA_Census.csv',keep_default_na=False);co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');targets=set(co.loc[co.gene.ne('(Intercept)')&co.coefficient.ne(0),'gene']);genes={e for s in m.loc[m.model_ensembl.isin(targets),'mapped_feature_ids'] for e in s.split('|') if e};ram=pd.read_csv(O/'RAM_mapped_gene_panels.csv',keep_default_na=False);genes|={e for s in ram.mapped_feature_ids for e in s.split('|') if e};v=v[v['_index'].isin(genes)].sort_values('soma_joinid').copy()
with h5py.File(p,'r') as h:
 z=h['X'];shape=tuple(z.attrs['shape']);assert shape==(len(obs),len(v));data=z['data'][:];ind=z['indices'][:];ptr=z['indptr'][:];assert len(ptr)==len(obs)+1 and ptr[0]==0 and ptr[-1]==len(data)==len(ind) and np.diff(ptr).min()>=0;assert np.all(ind>=0) and np.all(ind<len(v));x=sparse.csr_matrix((data,ind,ptr),shape=shape);assert z.attrs['encoding-type']=='csr_matrix'
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data));assert np.all(np.asarray(x.sum(axis=1)).ravel()<=obs.raw_sum.to_numpy()+1e-5)
# The API's obs/var joinids are ordered; requested coordinates were ascending.
obs.index=obs.soma_joinid.astype(str);obs.index.name=None;v=v.rename(columns={'_index':'feature_id'});v.index=v.feature_id;v.index.name=None;out=ad.AnnData(x,obs=obs,var=v);out.write_h5ad(dest.with_suffix('.recovered.h5ad'),compression='gzip');dest.with_suffix('.recovered.h5ad').rename(dest)
report={'status':'RECOVERED_COMPLETE_X_AFTER_KNOWN_OBS_INDEX_WRITE_ERROR','original_partial_file':str(p.relative_to(P)),'output_file':str(dest.relative_to(P)),'expression_values_unchanged':True,'n_cells':len(obs),'n_features':len(v),'X_values_sha256':hashlib.sha256(x.data.tobytes()+x.indices.tobytes()+x.indptr.tobytes()).hexdigest(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'additional_verification_required':'Compare recovered var-ID order to the next native Census AnnData chunk; metadata keys are the exact prespecified first100000 target coordinates.'};(O/'first_count_chunk_recovery.json').write_text(json.dumps(report,indent=2));print(report)
