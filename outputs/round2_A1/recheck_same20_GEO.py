#!/usr/bin/env python3
"""Isolate feature repair/rounding on immutable round1 cells, labels, QC and RNG."""
from pathlib import Path
import json,sys
import numpy as np,pandas as pd,h5py
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';sys.path.insert(0,str(O));from run_G1_prime import rng_for
from compute_reliability import icc_sem
co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');mp=pd.read_csv(O/'gene_mapping_GEO.csv',keep_default_na=False).set_index('model_ensembl');qc=pd.read_csv(P/'outputs/K4_gene_QC.csv');keepgene=qc.retained.astype(str).str.lower().eq('true').to_numpy();manifest=json.loads((P/'outputs/K4_download_manifest.json').read_text())['files'];prior=pd.read_csv(P/'outputs/K4_predictions.csv').set_index(['gsm','celltype']);cross={'CD4 Memory':'CD4T','CD4 Naive':'CD4T','CD4 Treg':'CD4T','Treg':'CD4T','CD8 effector':'CD8T','CD8 Naive':'CD8T','CD14+ Monocytes':'MONO','CD16+ Monocytes':'MONO','NK cell':'NK','pre-B cell':'B','B cell progenitor':'B','B cell':'B'}
def point(y,donor,gsm,ct):
 rng=rng_for(donor,gsm,ct,'point');return float(np.mean([y[rng.choice(len(y),15,replace=len(y)<=15)].mean() for _ in range(100)]))
rows=[]
for rec in manifest:
 with h5py.File(P/'work/K4_data'/rec['filename'],'r') as h:
  z=h['matrix'];shape=tuple(map(int,z['shape'][:]));x=sparse.csc_matrix((z['data'][:].astype(np.float64),z['indices'][:],z['indptr'][:]),shape=shape);sym=z['features/name'][:].astype(str);ids=z['features/id'][:].astype(str);ids=np.array([v.split('.')[0] for v in ids]);labels=z['observations/seurat_pbmc_type'][:].astype(str)
 assert np.array_equal(ids,qc.gene.to_numpy());lib=np.asarray(x.sum(axis=0)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[np.char.startswith(sym,'MT-')].sum(axis=0)).ravel();keep=(ng>=200)&(ng<=5000)&(lib>0)&(mt/np.maximum(lib,1)<=.25);x=x[keepgene,:][:,keep].T.tocsr();ids=ids[keepgene];lib=np.asarray(x.sum(axis=1)).ravel();broad=pd.Series(labels[keep]).map(cross).fillna('').to_numpy();fi={e:j for j,e in enumerate(ids)};log=x.multiply(10000/lib[:,None]).tocsr();log.data=np.log1p(log.data)
 for ct,g in co[co.gene.ne('(Intercept)')&co.coefficient.ne(0)].groupby('celltype'):
  ii=np.flatnonzero(broad==ct);inter=float(co.loc[co.celltype.eq(ct)&co.gene.eq('(Intercept)'),'coefficient'].iloc[0]);avail=[(j,e) for j,e in enumerate(g.gene) if e in fi];oldy=np.asarray(log[ii][:,[fi[e] for j,e in avail]]@g.coefficient.to_numpy()[[j for j,e in avail]]).ravel()+inter
  oldpoint=point(oldy,rec['donor'],rec['gsm'],ct);oldref=float(prior.loc[(rec['gsm'],ct),'predicted_age_native_unrounded']);assert abs(oldpoint-oldref)<1e-8,(rec['gsm'],ct,oldpoint,oldref)
  rr=[];cc=[];covered=set()
  for j,e in enumerate(g.gene):
   features=[f for f in mp.loc[e,'mapped_feature_ids'].split('|') if f in fi]
   if features:covered.add(e)
   for f in features:rr.append(fi[f]);cc.append(j)
  A=sparse.csc_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ids),len(g)));mapped=(x[ii]@A).tocsr();mapped=mapped.multiply(10000/lib[ii,None]).tocsr();mapped.data=np.log1p(mapped.data);y=np.asarray(mapped@g.coefficient.to_numpy()).ravel()+inter;newpoint=point(y,rec['donor'],rec['gsm'],ct)
  rows.append({'gsm':rec['gsm'],'donor':rec['donor'],'celltype':ct,'n_cells':len(y),'old_unrounded':oldpoint,'old_rounded':float(np.round(oldpoint)),'repaired_unrounded':newpoint,'repaired_rounded':float(np.round(newpoint)),'old_missing':int(prior.loc[(rec['gsm'],ct),'missing_nonzero_features']),'repaired_missing':len(g)-len(covered),'old_native_reproduction_abs_error':abs(oldpoint-oldref)})
 print('SAME20',rec['gsm'],flush=True)
a=pd.DataFrame(rows);a.to_csv(O/'same20_feature_repair_predictions.csv',index=False);pairs=pd.read_csv(P/'outputs/K4_selected_pairs.csv');out=[]
for ct,g in a.groupby('celltype'):
 for metric in ['old_rounded','old_unrounded','repaired_rounded','repaired_unrounded']:
  z=g.set_index('gsm')[metric];v=np.array([[z[r.gsm1],z[r.gsm2]] for r in pairs.itertuples()]);q=icc_sem(v);out.append({'celltype':ct,'analysis':metric,'n_donors':20,'ICC2_1':q[0],'SEM':q[1],'MDC95':q[2],'max_missing_nonzero_genes':int(g.old_missing.max() if metric.startswith('old') else g.repaired_missing.max()),'scope':'same original 20 donors, original GEO/QC/labels/denominator/pseudocell membership; only feature mapping and rounding differ'})
pd.DataFrame(out).to_csv(O/'same20_repair_and_rounding_comparison.csv',index=False);print('SAME20_COMPLETE max native point error',a.old_native_reproduction_abs_error.max(),flush=True)
