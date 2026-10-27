#!/usr/bin/env python3
"""Freeze50 age-stratified splits and cache OneK-only sufficient statistics."""
from pathlib import Path
import sys,json,gc
import numpy as np,pandas as pd
from sklearn.model_selection import StratifiedKFold
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';R=P/'outputs/round3';V=P/'work/round3';D=V/'training_design';C=W/'OneK_c_stats';C.mkdir(exist_ok=True);sys.path.insert(0,str(R));import train_repair_models as tr
def splits():
 sp=pd.read_csv(R/'R2_K0_donor_split.csv').sort_values('donor').reset_index(drop=True);rows=[]
 bands=(sp.age_at_first_eligible_visit>=45).astype(int)
 for rep in range(50):
  state=int(np.random.SeedSequence([20261025,4401,rep]).generate_state(1)[0]);sk=StratifiedKFold(2,shuffle=True,random_state=state)
  for fold,(_,test) in enumerate(sk.split(sp,bands)):
   for i in test:rows.append({'repeat':rep,'heldout_fold':fold,'donor':sp.donor.iloc[i],'age_band':bands.iloc[i],'split_random_state':state,'training_fold':1-fold})
 out=pd.DataFrame(rows);assert len(out)==3400 and out.groupby(['repeat','heldout_fold']).size().eq(34).all();out.to_csv(O/'c_50x2_donor_fold_ledger.csv',index=False)
def stats():
 ranks=json.loads((R/'R2_training_gene_ranks.json').read_text());models=pd.read_csv(R/'R2_model_training_results.csv');ref=pd.read_csv(R/'R2_reference_genes.csv')
 with threadpool_limits(limits=4):
  for ct in ['B','CD4T','CD8T','MONO','NK']:
   if (C/f'{ct}_complete.json').exists():print('REUSE_ONEK_STATS',ct,flush=True);continue
   m=pd.read_csv(D/f'{ct}_donor_metadata.csv');dg=pd.read_csv(D/f'{ct}_design_genes.csv');lookup={int(k):i for i,k in enumerate(dg.reference_index)};X=np.load(D/f'{ct}_pseudocells.npy',mmap_mode='r');means=np.load(D/f'{ct}_pseudocell_donor_means.npy');age=m.age.to_numpy();training=m.split.eq('train80').to_numpy();K=int(models[models.celltype.eq(ct)&models.variant.eq('a')].baseline_selected_K.iloc[0])
   for fold in ['fold0','fold1','fold2','fold3','fold4','full']:
    keep=training.copy();val=np.zeros(len(m),bool)
    if fold!='full':keep &=m.CV_fold.ne(int(fold[-1])).to_numpy();val=training&m.CV_fold.eq(int(fold[-1])).to_numpy()
    rr=np.flatnonzero(np.repeat(keep,100));ix=np.array(ranks[ct][fold][:K],int);pos=[lookup[i] for i in ix];z=np.asarray(X[np.ix_(rr,pos)]);s=tr.sufficient(z,np.repeat(age[keep],100),age[keep]);np.savez(C/f'{ct}_{fold}.npz',**s,reference_index=ix,validation_mean=means[val][:,pos],validation_age=age[val],test_mean=means[~training][:,pos],test_age=age[~training],test_donor=m.donor[~training].to_numpy(dtype=str),n_training_donors=int(keep.sum()));print('ONEK_STATS',ct,fold,len(ix),flush=True);del z;gc.collect()
   (C/f'{ct}_complete.json').write_text(json.dumps({'status':'COMPLETE','K':K,'scope':'OneK fixed train80 only; no Allen/AIDA outcomes'},indent=2))
if __name__=='__main__':splits();stats()
