#!/usr/bin/env python3
"""Exact shared *storage* for two separately frozen populations; no pooled analysis.
Use official full-row blockwise reads to avoid large two-axis point projections.
Every full row sum is checked before retaining the fixed model gene projection.
"""
from pathlib import Path
import json,hashlib,time
import numpy as np,pandas as pd,anndata as ad,cellxgene_census as c
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';D=W/'union_count_cache';D.mkdir(exist_ok=True)
CFG={'soma.init_buffer_bytes':134217728,'py.init_buffer_bytes':134217728,'sm.mem.total_budget':17179869184,'sm.compute_concurrency_level':16,'sm.io_concurrency_level':32}

def main():
 qa=json.loads((O/'Census_blockwise_equivalence_QA.json').read_text());assert qa['status']=='PASS_EXACT_RAW_COUNT_EQUIVALENCE'
 cols=['soma_joinid','gsm','donor','raw_sum','AIFI_L2','AIFI_L3','dataset_id','clock_type'];a=pd.read_parquet(W/'cache/HIHA_K0_clock_cells.parquet',columns=cols);b=pd.read_parquet(W/'sensitivity87/cache/HIHA_K0_clock_cells.parquet',columns=cols);pa=set(a.gsm);pb=set(b.gsm)
 overlap=a.merge(b,on='soma_joinid',suffixes=('_primary','_sensitivity'),validate='one_to_one')
 for key in cols[1:]:assert overlap[key+'_primary'].eq(overlap[key+'_sensitivity']).all(),key
 obs=pd.concat([a,b],ignore_index=True).drop_duplicates('soma_joinid').sort_values('soma_joinid').reset_index(drop=True);assert len(a)==4609769 and len(b)==3522180 and len(obs)==6357549
 obs['in_primary_K0']=obs.gsm.isin(pa);obs['in_sensitivity87']=obs.gsm.isin(pb);obs.to_parquet(W/'cache/union_clock_cells.parquet',index=False)
 co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');m=pd.read_csv(O/'gene_mapping_HIHA_Census.csv',keep_default_na=False);targets=set(co.loc[co.gene.ne('(Intercept)')&co.coefficient.ne(0),'gene']);genes={e for s in m.loc[m.model_ensembl.isin(targets),'mapped_feature_ids'] for e in s.split('|') if e};ram=pd.read_csv(O/'RAM_mapped_gene_panels.csv',keep_default_na=False);genes|={e for s in ram.mapped_feature_ids for e in s.split('|') if e};v=pd.read_csv(O/'HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv');v=v[v['_index'].isin(genes)].sort_values('soma_joinid').rename(columns={'_index':'feature_id'});assert len(v)==2957;vid=v.soma_joinid.to_numpy(dtype=np.int64)
 plan={'seed':20261025,'census_release':'2025-11-08','storage_only_union':True,'no_pooled_89_donor_analysis':True,'n_union_cells':len(obs),'n_shared_cells':len(overlap),'n_features':len(v),'primary':{'n_donors':68,'n_visits':320,'n_cells':len(a)},'sensitivity87':{'n_donors':87,'n_visits':236,'n_cells':len(b)},'normalization_denominator':'full Census raw_sum; full-row sums verified before projection','read_config':CFG,'query':'official SparseNDArray.read(coords=(exact wanted soma_joinids, slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1)','version_specific_return_order':'TileDB-SOMA2.3 returns (scipy_matrix,(row_joinids,column_joinids)), verified from installed primary source','population_separation':'Each scoring script filters by its immutable GSM ledger, and each population manifest counts only its selected rows.'};(O/'union_count_cache_plan.json').write_text(json.dumps(plan,indent=2));print('PLAN',plan,flush=True)
 manifests=[];recover= None;recoverqa=[]
 rp=W/'sensitivity87/main_Census_counts/HIHA_K0_000.h5ad'
 if rp.exists():recover=ad.read_h5ad(rp);recover.obs.index=recover.obs.soma_joinid.astype(str);recover.obs.index.name=None;assert np.array_equal(recover.var.soma_joinid.to_numpy(),vid)
 with c.open_soma(census_version='2025-11-08',tiledb_config=CFG) as f:
  X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw'];print('OPEN_RAW',flush=True)
  for i,lo in enumerate(range(0,len(obs),50000)):
   rr=obs.iloc[lo:lo+50000].copy();path=D/f'HIHA_union_{i:03d}.h5ad';started=time.monotonic()
   if not path.exists():
    chunks=[];rowids=[];expected=rr.set_index('soma_joinid');print('FETCH',i,len(rr),flush=True)
    for mat,(rid,_) in X.read(coords=(rr.soma_joinid.to_numpy(dtype=np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
     rid=np.asarray(rid,dtype=np.int64);assert mat.shape[0]==len(rid) and mat.shape[1]>vid.max();full=np.asarray(mat.sum(axis=1)).ravel();assert np.array_equal(full,expected.loc[rid,'raw_sum'].to_numpy()),'Full-row count completeness failed';chunks.append(mat[:,vid].tocsr());rowids.extend(rid.tolist());print('BLOCK',i,len(rowids),flush=True)
    x=sparse.vstack(chunks).tocsr();ix=pd.Index(rowids).get_indexer(rr.soma_joinid);assert (ix>=0).all() and len(rowids)==len(rr) and len(set(rowids))==len(rowids);x=x[ix].tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data));assert np.all(np.asarray(x.sum(axis=1)).ravel()<=rr.raw_sum.to_numpy())
    rr.index=rr.soma_joinid.astype(str);rr.index.name=None;vv=v.copy();vv.index=vv.feature_id.astype(str);vv.index.name=None;obj=ad.AnnData(x,obs=rr,var=vv);obj.write_h5ad(path.with_suffix('.part.h5ad'),compression='gzip');path.with_suffix('.part.h5ad').rename(path)
   obj=ad.read_h5ad(path,backed='r');assert obj.shape==(len(rr),len(v));assert np.array_equal(obj.obs.soma_joinid.to_numpy(),rr.soma_joinid.to_numpy());assert np.array_equal(obj.var.soma_joinid.to_numpy(),vid)
   # Recovering a post-X writer error is only trusted after independent raw
   # counts for the same coordinates have been compared. Recovered data are
   # NOT used as the source of the union cache.
   if recover is not None:
    keep=obj.obs.soma_joinid.isin(recover.obs.soma_joinid).to_numpy()
    if keep.any():
     jj=np.flatnonzero(keep);ridx=pd.Index(recover.obs.soma_joinid).get_indexer(obj.obs.iloc[jj].soma_joinid);d=(obj.X[jj].tocsr()-recover.X[ridx]).tocsr();err=float(abs(d.data).max()) if d.nnz else 0.;recoverqa.append({'union_chunk':i,'n_overlap_cells':len(jj),'max_absolute_count_difference':err,'source':'fresh official blockwise raw rows vs recovered native get_anndata projection'});pd.DataFrame(recoverqa).to_csv(O/'recovered_count_chunk_independent_QA.csv',index=False)
   obj.file.close();row={'path':str(path.relative_to(P)),'physical_n_cells':len(rr),'n_primary_cells':int(rr.in_primary_K0.sum()),'n_sensitivity87_cells':int(rr.in_sensitivity87.sum()),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'seconds':time.monotonic()-started};manifests.append(row);(O/'union_count_cache_manifest.json').write_text(json.dumps({'complete':False,'plan':plan,'chunks':manifests},indent=2));print('DONE',i,row['seconds'],flush=True)
 assert sum(r['physical_n_cells'] for r in manifests)==len(obs)
 (O/'union_count_cache_manifest.json').write_text(json.dumps({'complete':True,'plan':plan,'chunks':manifests},indent=2))
 for folder,field,n in [(O,'n_primary_cells',len(a)),(O/'sensitivity87','n_sensitivity87_cells',len(b))]:
  chunks=[{**r,'n_cells':r[field],'selection':'filter by this population frozen GSM ledger'} for r in manifests if r[field]>0];assert sum(r['n_cells'] for r in chunks)==n
  pop=plan['primary' if field=='n_primary_cells' else 'sensitivity87'];(folder/'main_Census_count_manifest.json').write_text(json.dumps({'complete':True,'exact_union_storage_reuse':True,'plan':{'seed':20261025,'Census_release':'2025-11-08','n_cells':n,'n_features':len(v),**pop},'chunks':chunks},indent=2))
 print('ALL_UNION_COUNTS_COMPLETE',flush=True)
if __name__=='__main__':main()
