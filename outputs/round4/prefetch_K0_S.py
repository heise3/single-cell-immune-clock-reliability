#!/usr/bin/env python3
"""Disjoint chunk download workers; unchanged full-library QC and genes."""
from pathlib import Path
import sys,json,time,gc,fcntl
import numpy as np,pandas as pd,anndata as ad,cellxgene_census as census
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parent));import prepare_round4_inputs as prep
def main():
 worker=int(sys.argv[1]);k,obs,genes=prep.base();obs=obs.sort_values('soma_joinid').reset_index(drop=True);v=pd.read_csv(prep.R/'Census_global_var.csv');mtids=v.loc[v.feature_name.str.startswith('MT-',na=False),'soma_joinid'].to_numpy(np.int64);vid=genes.soma_joinid.to_numpy(np.int64);D=prep.W/'K0_S_evaluation_counts';D.mkdir(exist_ok=True);receipt=[];cfg={**prep.ev.CFG,'sm.compute_concurrency_level':4,'sm.io_concurrency_level':8}
 with census.open_soma(census_version='2025-11-08',tiledb_config=cfg) as f:
  X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw']
  for i,lo in enumerate(range(0,len(obs),40000)):
   if i%2!=worker:continue
   file=D/f'K0_S_{i:03d}.h5ad';qcfile=file.with_suffix('.QC.json');lock=open(D/f'chunk_{i:03d}.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
   if file.exists() and qcfile.exists():lock.close();print('S_REUSE',worker,i,flush=True);continue
   st=time.monotonic();rr=obs.iloc[lo:lo+40000];lk=rr.set_index('soma_joinid');xs=[];os=[]
   for mat,(rid,_) in X.read(coords=(rr.soma_joinid.to_numpy(np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
    rid=np.asarray(rid,np.int64);o=lk.loc[rid].reset_index();lib=np.asarray(mat.sum(axis=1)).ravel();assert np.array_equal(lib,o.raw_sum.to_numpy());nf=np.diff(mat.tocsr().indptr);mt=np.asarray(mat[:,mtids].sum(axis=1)).ravel()/lib;keep=(nf>=200)&(nf<=5000)&(mt<=.25);o=o[keep].copy();o['native_n_features']=nf[keep];o['native_mito_fraction']=mt[keep];os.append(o);xs.append(mat[keep][:,vid].tocsr())
   x=sparse.vstack(xs).tocsr();o=pd.concat(os,ignore_index=True);o.index=o.soma_joinid.astype(str);o.index.name=None;v=genes.copy();v.index=v.feature_id;v.index.name=None;a=ad.AnnData(x,obs=o,var=v);tmp=file.with_suffix(f'.worker{worker}.part.h5ad');a.write_h5ad(tmp,compression='lzf');tmp.rename(file);qc={'input_cells':len(rr),'kept_cells':len(o),'excluded_cells':len(rr)-len(o),'full_row_sum_exact':True};qcfile.write_text(json.dumps(qc,indent=2));lock.close();receipt.append({'path':str(file.relative_to(prep.P)),'sha256':prep.ag.sha(file),'bytes':file.stat().st_size,**qc});(prep.O/f'K0_S_prefetch_worker{worker}.json').write_text(json.dumps(receipt,indent=2));print('S_PREFETCH',worker,i,len(o),round(time.monotonic()-st,1),flush=True);del a,x,o;gc.collect()
 print('S_PREFETCH_COMPLETE',worker,flush=True)
if __name__=='__main__':main()
