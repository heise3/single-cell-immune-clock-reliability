#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np,pandas as pd,h5py
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';A=P/'outputs/round2_A1'
def main():
 assert (W/'short_gene_inputs_complete.json').exists();meta=pd.read_csv(W/'K0_short_groups.csv');genes=pd.read_csv(O/'c_candidate_gene_universe.csv');splits=pd.read_csv(O/'c_50x2_donor_fold_ledger.csv');pairs=pd.read_csv(A/'provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.stratum.eq('7-9d')];lookup={g:i for i,g in enumerate(meta.group)};allmask=[];allval=[];rows=[]
 with h5py.File(W/'K0_short_candidate_pseudobulk.h5') as h:
  x=h['pseudobulk_logCPM'][:]
  for ct in ['B','CD4T','CD8T','MONO','NK']:
   vals=[];masks=[]
   for rep in range(50):
    s=splits[splits.repeat.eq(rep)]
    for held in [0,1]:
     training=set(s.loc[s.heldout_fold.ne(held),'donor']);testing=set(s.loc[s.heldout_fold.eq(held),'donor']);assert len(training)==len(testing)==34 and not training&testing;pp=pairs[pairs.donor.isin(training)];aa=np.array([[x[lookup[a+'|'+ct]],x[lookup[b+'|'+ct]]] for a,b in zip(pp.gsm1,pp.gsm2)]);n=len(aa);assert n>=3;gm=aa.mean(axis=(0,1));rm=aa.mean(axis=1);cm=aa.mean(axis=0);msr=2*np.var(rm,axis=0,ddof=1);msc=n*np.var(cm,axis=0,ddof=1);mse=((aa-rm[:,None,:]-cm[None,:,:]+gm)**2).sum(axis=(0,1))/(n-1);den=msr+mse+2/n*(msc-mse);icc=np.divide(msr-mse,den,out=np.full(len(genes),np.nan),where=den>1e-14);mask=icc>=.7;vals.append(icc);masks.append(mask);rows.append({'celltype':ct,'repeat':rep,'heldout_fold':held,'n_training_donors':34,'n_heldout_donors':34,'n_training_short_pairs':n,'n_candidate_genes':len(genes),'n_ICC_ge_0_7':int(mask.sum()),'heldout_donors_used_for_ICC':False})
   assert 'reference_index' in genes
   np.savez_compressed(W/f'c_{ct}_ICC_masks.npz',ICC=np.asarray(vals),mask=np.asarray(masks),reference_index=genes.reference_index.to_numpy(),feature_id=genes.feature_id.to_numpy(dtype=str));print('C_MASKS',ct,np.asarray(masks).shape,flush=True)
 pd.DataFrame(rows).to_csv(O/'c_crossfit_mask_QA.csv',index=False);(W/'c_masks_complete.json').write_text(json.dumps({'status':'COMPLETE','models':500,'repetitions':50,'folds':2,'gene_ICC_uses_opposite_heldout_fold_only':True},indent=2))
if __name__=='__main__':main()
