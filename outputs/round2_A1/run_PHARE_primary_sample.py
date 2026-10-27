#!/usr/bin/env python3
"""One real K0 sample, native released PHARE preprocessing/annotation/prediction."""
from pathlib import Path
import sys,json,hashlib,random
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';SEED=20261025
random.seed(SEED);np.random.seed(SEED)
import scanpy as sc,torch
torch.manual_seed(SEED)
# Keep native clustering algorithms and explicitly supply the user-fixed RNG.
for owner,name in [(sc.pp,'neighbors'),(sc.pp,'pca'),(sc.tl,'leiden')]:
 original=getattr(owner,name)
 def seeded(*args,_original=original,**kwargs):
  kwargs['random_state']=SEED
  return _original(*args,**kwargs)
 setattr(owner,name,seeded)
sys.path.insert(0,str(P/'work/round2/PHARE_runtime'))
from phare.utils import load_model
r=json.loads((O/'PHARE_input_manifest.json').read_text());p=P/r['input_path'];k0=pd.read_csv(P/'outputs/round2/K0_primary_with_released_age.csv').set_index('gsm');a=ad.read_h5ad(p);x=a.X.tocsr().astype(np.float64);sym=a.var.feature_name.astype(str).to_numpy();sample=[r['specimen']]
# Sum duplicate symbols in raw counts, preserving full library total.
unique=pd.Index(pd.unique(sym));j=unique.get_indexer(sym);A=sparse.csc_matrix((np.ones(len(sym)),(np.arange(len(sym)),j)),shape=(len(sym),len(unique)));x=x@A
obs=pd.DataFrame({'orig.ident':[r['gsm']]*x.shape[0]},index=[str(i) for i in range(x.shape[0])]);a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=unique));assert a.X.sum()>0
modeldir=P/'work/round2/PHARE_runtime/phare/model/fast';model=load_model(str(modeldir));a=model.preprocess(a);pred=model.predict(a,sample_column='orig.ident');assert len(pred)==1 and np.isfinite(pred.select_dtypes(include='number').to_numpy()).all();assert np.isfinite(pred.filter(regex='^predicted_').to_numpy(dtype=float)).all();pred.to_csv(O/'PHARE_real_K0_prediction.csv')
a.obs.to_csv(O/'PHARE_real_K0_cell_annotation.csv',index=False)
report={'seed':SEED,'gsm':r['gsm'],'donor':k0.loc[r['gsm'],'donor'],'specimen':str(sample[0]),'n_cells':a.n_obs,'n_unique_gene_symbols':a.n_vars,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'native_preprocess':'normalize_total target_sum=10000 then log1p','native_annotation':'published CellTypist fairSample classifier; native fast_annotation flag '+str(model.fast_annotation),'randomness':'Python/NumPy/Torch 20261025; any native Scanpy PCA/neighbors/Leiden random_state explicitly 20261025','prediction':pred.reset_index().to_dict('records'),'status':'RUN_COMPLETED_REAL_MAIN_COHORT_SAMPLE','limitation':'One-sample execution validation; no claim of PHARE reliability or accuracy validation'};(O/'PHARE_real_K0_run.json').write_text(json.dumps(report,indent=2,default=str));print(report,flush=True)
