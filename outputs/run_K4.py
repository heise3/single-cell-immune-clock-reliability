#!/usr/bin/env python3
"""Frozen 20-donor K4, sc-ImmuAging pretrained LASSO, seed 20261025.

Native algorithm: RNA@data = log1p(CP10K), 100 pseudocells formed by
averaging 15 cells, predict at lambda.min, round mean donor prediction.
Linearity permits calculating each cell's affine score first; averaging
those scores is algebraically identical to predicting averaged features.
This equivalence is checked against native glmnet by verify_linear.R.

Outer cell bootstrap: 100 resamples of N cells with replacement per
donor/visit/type. Both the conditional mean and a fresh 100x15 native
pseudocell estimate are retained. The conditional-mean variance isolates
cell sampling from the extra Monte Carlo randomness of pseudocells.
"""
from pathlib import Path
import hashlib,json,csv,collections
import numpy as np
import pandas as pd
import h5py
from scipy import sparse

ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent/'work';DATA=WORK/'K4_data'
SEED=20261025;B=100
TYPES=['CD4T','CD8T','MONO','NK','B']
CROSSWALK={
 'CD4 Memory':'CD4T','CD4 Naive':'CD4T','CD4 Treg':'CD4T','Treg':'CD4T',
 'CD8 effector':'CD8T','CD8 Naive':'CD8T',
 'CD14+ Monocytes':'MONO','CD16+ Monocytes':'MONO',
 'NK cell':'NK','pre-B cell':'B','B cell progenitor':'B','B cell':'B',
}
manifest=json.loads((ROOT/'K4_download_manifest.json').read_text())['files']
assert len(manifest)==40 and not any('error' in r for r in manifest)
meta=pd.read_csv(ROOT/'sample_metadata.csv',keep_default_na=False).set_index('gsm')
co=pd.read_csv(ROOT/'scImmuAging_coefficients.csv')

def rng_for(*args):
 entropy=int.from_bytes(hashlib.sha256('|'.join(map(str,args)).encode()).digest()[:8],'little')
 return np.random.default_rng(np.random.SeedSequence([SEED,entropy]))

def load_count(item):
 with h5py.File(DATA/item['filename'],'r') as f:
  g=f['matrix'];shape=tuple(map(int,g['shape'][:]))
  x=sparse.csc_matrix((g['data'][:].astype(np.float64),g['indices'][:].astype(np.int32),g['indptr'][:].astype(np.int64)),shape=shape)
  symbols=g['features/name'][:].astype(str)
  names=np.array([v.split('.')[0] for v in g['features/id'][:].astype(str)])
  labels=g['observations/seurat_pbmc_type'][:].astype(str)
  sample_ids=np.unique(g['observations/pbmc_sample_id'][:].astype(str)).tolist()
  assert len(sample_ids)==1,(item['gsm'],sample_ids)
  assert all(s in item['filename'] for s in sample_ids),(item['gsm'],sample_ids)
 ng=np.diff(x.indptr);total=np.asarray(x.sum(axis=0)).ravel()
 mito=np.asarray(x[np.char.startswith(symbols,'MT-')].sum(axis=0)).ravel()
 keep=(ng>=200)&(ng<=5000)&(total>0)&(mito/np.maximum(total,1)<=.25)
 return x,names,labels,keep,sample_ids[0]

# Gene support is determined across the frozen pilot, not re-fit per visit.
support=None;gene_names=None;qc=[];label_counts=collections.Counter()
for i,item in enumerate(manifest):
 x,genes,labels,keep,sample_id=load_count(item)
 if support is None:gene_names=genes;support=np.zeros(len(genes),dtype=np.int64)
 assert np.array_equal(genes,gene_names)
 support+=np.asarray((x[:,keep]!=0).sum(axis=1)).ravel().astype(np.int64)
 qc.append({**{k:item[k] for k in ['donor','gsm','visit']},'sample_id':sample_id,'cells_raw':x.shape[1],'cells_pass_qc':int(keep.sum())})
 label_counts.update(labels[keep])
 print('QC',i+1,item['gsm'],int(keep.sum()),flush=True)
pd.DataFrame(qc).to_csv(ROOT/'K4_sample_QC.csv',index=False)
gene_keep=support>=3
pd.DataFrame({'gene':gene_names,'cells_detected':support,'retained':gene_keep}).to_csv(ROOT/'K4_gene_QC.csv',index=False)

def native_pseudocell_age(scores,rng):
 n=len(scores);replace=n<=15
 pseudo=np.array([scores[rng.choice(n,15,replace=replace)].mean() for _ in range(100)])
 return pseudo.mean()

predictions=[];bootstrap_rows=[];validation=[];score_cache={}
for i,item in enumerate(manifest):
 x,genes,labels,keep,sample_id=load_count(item)
 x=x[gene_keep,:][:,keep].tocsr();labels=labels[keep];genes=genes[gene_keep]
 totals=np.asarray(x.sum(axis=0)).ravel();x=x.tocsc(copy=False)
 x.data*=np.repeat(10000/totals,np.diff(x.indptr));x.data=np.log1p(x.data)
 x=x.tocsr(); gene_index={g:j for j,g in enumerate(genes)}
 broad=np.array([CROSSWALK.get(label,'UNMAPPED') for label in labels])
 for celltype in TYPES:
  beta=co[co.celltype==celltype].set_index('gene').coefficient
  intercept=float(beta['(Intercept)']);beta=beta.drop('(Intercept)')
  available=[g for g in beta.index if g in gene_index]
  assert len(available)>0,('ENGINEERING_ERROR_NO_MATCHED_MODEL_FEATURES',item['gsm'],celltype)
  missing=[g for g in beta.index if g not in gene_index and beta[g]!=0]
  selected=np.flatnonzero(broad==celltype);assert len(selected)>0,(item['gsm'],celltype)
  sub=x[[gene_index[g] for g in available],:][:,selected]
  scores=np.asarray(sub.T@beta.loc[available].to_numpy()).ravel()+intercept
  assert np.all(np.isfinite(scores))
  rng=rng_for(item['donor'],item['gsm'],celltype,'point')
  unrounded=native_pseudocell_age(scores,rng);rounded=float(np.round(unrounded))
  rng=rng_for(item['donor'],item['gsm'],celltype,'bootstrap')
  cell_means=[];native_boot=[]
  for b in range(B):
   resampled=scores[rng.integers(0,len(scores),size=len(scores))]
   cm=resampled.mean();nb=native_pseudocell_age(resampled,rng)
   cell_means.append(cm);native_boot.append(nb)
   bootstrap_rows.append({'donor':item['donor'],'gsm':item['gsm'],'visit':item['visit'],'celltype':celltype,'replicate':b+1,'cell_conditional_mean':cm,'native_pseudocell_mean':nb,'native_rounded':float(np.round(nb))})
  predictions.append({'donor':item['donor'],'gsm':item['gsm'],'visit':item['visit'],'celltype':celltype,'n_cells':len(scores),'predicted_age_native':rounded,'predicted_age_native_unrounded':unrounded,'predicted_age_exact_expectation':scores.mean(),'cell_boot_variance':np.var(cell_means,ddof=1),'native_boot_variance':np.var(native_boot,ddof=1),'missing_nonzero_features':len(missing),'nonzero_features':int((beta!=0).sum()),'missing_feature_names':'|'.join(missing)})
  score_cache[item['gsm']+'_'+celltype]=scores
  if i==0:
   # Independent native glmnet validation uses actual normalized cells.
   for j in range(min(3,len(selected))):
    vals=dict(zip(available,sub[:,j].toarray().ravel()))
    validation.append({'celltype':celltype,'row':j,'python_prediction':scores[j],**vals})
 print('PRED',i+1,item['gsm'],flush=True)
pd.DataFrame(predictions).to_csv(ROOT/'K4_predictions.csv',index=False)
pd.DataFrame(bootstrap_rows).to_csv(ROOT/'K4_bootstrap100.csv',index=False)
pd.DataFrame(validation).to_csv(ROOT/'K4_linear_validation_inputs.csv',index=False)
np.savez_compressed(WORK/'K4_cell_scores.npz',**score_cache)

def icc2_1(y):
 n,k=y.shape;mean=y.mean();row=y.mean(axis=1);col=y.mean(axis=0)
 msr=k*np.sum((row-mean)**2)/(n-1)
 msc=n*np.sum((col-mean)**2)/(k-1)
 mse=np.sum((y-row[:,None]-col[None,:]+mean)**2)/((n-1)*(k-1))
 value=(msr-mse)/(msr+(k-1)*mse+k*(msc-mse)/n)
 sd=np.std(y.ravel(),ddof=1);sem=sd*np.sqrt(1-value)
 return dict(n_donors=n,ICC_2_1=value,SD_pooled=sd,SEM=sem,MDC95=1.96*np.sqrt(2)*sem,MS_donor=msr,MS_visit=msc,MS_error=mse)
frame=pd.DataFrame(predictions);results=[]
for celltype in TYPES:
 t=frame[frame.celltype==celltype]
 y=t.pivot(index='donor',columns='visit',values='predicted_age_native').reindex(columns=['Immune Variation Day 0','Immune Variation Day 7']).to_numpy()
 assert y.shape==(20,2) and np.all(np.isfinite(y))
 stats=icc2_1(y);error_var=stats['SEM']**2
 results.append({'celltype':celltype,**stats,'cells_min':int(t.n_cells.min()),'cells_median':float(t.n_cells.median()),'cell_sampling_variance_mean':float(t.cell_boot_variance.mean()),'native_boot_variance_mean':float(t.native_boot_variance.mean()),'cell_sampling_to_SEM2_ratio':float(t.cell_boot_variance.mean()/error_var),'native_boot_to_SEM2_ratio':float(t.native_boot_variance.mean()/error_var),'missing_features_max':int(t.missing_nonzero_features.max())})
pd.DataFrame(results).to_csv(ROOT/'K4_reliability.csv',index=False)
(ROOT/'K4_method.json').write_text(json.dumps({'seed':SEED,'n_donors':20,'bootstrap_replicates':B,'crosswalk':CROSSWALK,'source_label_counts_after_QC':dict(label_counts),'QC':'200<=features<=5000, mitochondrial UMI fraction<=25%; genes in >=3 post-QC cells across frozen pilot','normalization':'natural log1p(CP10K), all retained genes in denominator; RNA data slot, no scale.data used','prediction':'public pretrained LASSO at lambda.min; mean of 100 random 15-cell pseudocells per visit/type, rounded as native Age_Donor','SEM':'pooled SD of the 40 donor-visit predictions * sqrt(1-ICC(2,1))','noise_ratio':'mean conditional-mean cell-bootstrap variance divided by SEM^2. This is a diagnostic ratio, not a causal variance decomposition; native_boot ratio additionally includes pseudocell Monte Carlo randomness.','eligibility':'10 young and 10 older donors; both IV0-IV7 visits before first recorded COVID dose; no new vaccine-washout rule','limitations':['provided GEO preliminary Seurat lineage labels, mapped explicitly to five clock lineages; no reannotation or cell-label validation','pilot contains separated young and older age bands; pooled ICC can be high from age-range spread','100 bootstrap replicates estimate conditional cell sampling; no new-draw technical replicate is available']},ensure_ascii=False,indent=2))
print(pd.DataFrame(results).to_string(index=False),flush=True)
