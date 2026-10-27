#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,gc
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';D=W/'training_design';D.mkdir(exist_ok=True);TYPES=['B','CD4T','CD8T','MONO','NK'];SEED=20261025
def rng_for(*keys):
 h=int.from_bytes(hashlib.sha256('|'.join(map(str,keys)).encode()).digest()[:8],'little');return np.random.default_rng(np.random.SeedSequence([SEED,h]))
def rank_gene(x,y,nexpressed,genes):
 xx=x-x.mean(axis=0);yy=y-y.mean();den=np.sqrt(np.square(xx).sum(axis=0)*np.square(yy).sum());r=np.divide(yy@xx,den,out=np.zeros(x.shape[1]),where=den>0);r=np.clip(r,-1+1e-14,1-1e-14);pv=2*stats.t.sf(abs(r)*np.sqrt((len(y)-2)/(1-r*r)),len(y)-2);pv[(den==0)|(nexpressed<3)]=1;q=multipletests(pv,method='fdr_bh')[1];ix=np.flatnonzero(q<.05);ix=ix[np.lexsort((genes[ix],-abs(r[ix])))];return ix,r,q
def prepare_ranks():
 meta=pd.read_csv(W/'OneK1K_groups.csv');genes=pd.read_csv(O/'R2_reference_genes.csv');ranks={};summary=[];rows=[]
 with h5py.File(W/'OneK1K_aggregates.h5') as h:
  for ct in TYPES:
   cm=meta[meta.celltype.eq(ct)].copy();idx=cm.index.to_numpy();means=h['mean_CP10K'][idx];nc=h['expressing_cells'][idx];y=cm.age.to_numpy();training=cm.split.eq('train80').to_numpy();ranks[ct]={}
   for scope in ['full','fold0','fold1','fold2','fold3','fold4']:
    keep=training.copy()
    if scope!='full':keep &= cm.CV_fold.ne(int(scope[-1])).to_numpy()
    ix,rr,q=rank_gene(means[keep],y[keep],nc[keep].sum(axis=0),genes.feature_id.to_numpy());ranked=ix[:2000];ranks[ct][scope]=ranked.tolist();summary.append({'celltype':ct,'scope':scope,'n_training_donors':int(keep.sum()),'n_BH_significant':len(ix),'n_ranked_for_grid':len(ranked)})
    rows.extend({'celltype':ct,'scope':scope,'rank':rank+1,'reference_index':int(j),'feature_id':genes.feature_id.iloc[j],'Pearson_r_training':rr[j],'BH_q_training':q[j]} for rank,j in enumerate(ranked))
   union=sorted(set(j for z in ranks[ct].values() for j in z));genes.iloc[union].reset_index().rename(columns={'index':'reference_index'}).to_csv(D/f'{ct}_design_genes.csv',index=False);cm.reset_index(drop=True).to_csv(D/f'{ct}_donor_metadata.csv',index=False)
 (O/'R2_training_gene_ranks.json').write_text(json.dumps(ranks));pd.DataFrame(summary).to_csv(O/'R2_training_feature_selection_QA.csv',index=False);pd.DataFrame(rows).to_csv(O/'R2_training_ranked_genes.csv',index=False);return meta,ranks
def make_pseudocells():
 if (D/'complete.json').exists():print('REUSE_COMPLETE_DESIGNS',flush=True);return
 manifest=json.loads((O/'R2_OneK1K_raw_manifest.json').read_text());assert manifest['complete'];meta=pd.read_csv(W/'OneK1K_groups.csv');obs=[]
 for row in manifest['chunks']:
  a=ad.read_h5ad(P/row['path'],backed='r');obs.append(a.obs[['soma_joinid','group','celltype','raw_sum']].reset_index(drop=True));a.file.close()
 obs=pd.concat(obs,ignore_index=True);obs.to_parquet(D/'OneK1K_QC_master_obs.parquet',index=False);groups=obs.groupby('group',sort=False).indices;members={};design={};geneindices={};checkpoint=D/'pseudocell_progress.json';start_chunk=0 # Incomplete additive matrices are rebuilt to avoid double counting after a crash.
 for ct in TYPES:
  cm=pd.read_csv(D/f'{ct}_donor_metadata.csv');gg=pd.read_csv(D/f'{ct}_design_genes.csv');geneindices[ct]=gg.reference_index.to_numpy();mr=[];mc=[]
  for i,key in enumerate(cm.group):
   cells=np.asarray(groups[key]);rng=rng_for('OneK1K_R3_pseudocell',key);draw=np.array([rng.choice(cells,15,replace=len(cells)<=15) for _ in range(100)]);mr.extend(np.repeat(np.arange(i*100,(i+1)*100),15));mc.extend(draw.ravel())
  members[ct]=sparse.csr_matrix((np.full(len(mr),1/15.),(mr,mc)),shape=(len(cm)*100,len(obs)));assert np.allclose(np.asarray(members[ct].sum(axis=1)).ravel(),1);sparse.save_npz(D/f'{ct}_membership.npz',members[ct]);path=D/f'{ct}_pseudocells.npy'
  if start_chunk:design[ct]=np.lib.format.open_memmap(path,mode='r+')
  else:design[ct]=np.lib.format.open_memmap(path,mode='w+',dtype=np.float64,shape=(len(cm)*100,len(gg)));design[ct][:]=0
 lo=0
 for i,row in enumerate(manifest['chunks']):
  n=row['kept_cells'];hi=lo+n
  if i<start_chunk:lo=hi;continue
  a=ad.read_h5ad(P/row['path']);assert np.array_equal(a.obs.soma_joinid.to_numpy(),obs.soma_joinid.iloc[lo:hi].to_numpy());x=a.X.astype(np.float64).multiply((10000/a.obs.raw_sum.to_numpy())[:,None]).tocsr();x.data=np.log1p(x.data)
  for ct in TYPES:
   z=(members[ct][:,lo:hi]@x[:,geneindices[ct]]).tocoo();design[ct][z.row,z.col]+=z.data;design[ct].flush()
  lo=hi;checkpoint.write_text(json.dumps({'completed_chunks':i+1,'exact_sample_membership_shared_all_variants':True},indent=2));print('PSEUDOCELL_CHUNK',i+1,'/',len(manifest['chunks']),flush=True);del a,x;gc.collect()
 assert lo==len(obs)
 for ct in TYPES:
  assert np.isfinite(design[ct]).all();np.save(D/f'{ct}_pseudocell_donor_means.npy',design[ct].reshape(-1,100,design[ct].shape[1]).mean(axis=1));print('DESIGN',ct,design[ct].shape,flush=True)
 (D/'complete.json').write_text(json.dumps({'status':'COMPLETE','seed':SEED,'n_valid_OneK1K_cells':len(obs),'E_or_AIDA_outcomes_used':False},indent=2))
def main():
 assert (W/'OneK1K_aggregation_complete.json').exists();prepare_ranks();make_pseudocells();print('TRAINING_DESIGNS_COMPLETE',flush=True)
if __name__=='__main__':main()
