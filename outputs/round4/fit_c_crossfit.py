#!/usr/bin/env python3
"""Fit secondary c on OneK only, using each opposite Allen half's ICC mask."""
from pathlib import Path
import sys,json,os,time
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';R=P/'outputs/round3';V=P/'work/round3';C=W/'OneK_c_stats';M=W/'c_models';M.mkdir(exist_ok=True);sys.path.insert(0,str(R));import train_repair_models as tr
def solve(caches,mask):
 total=np.zeros(100);maxgap=0;ns=[];nval=0
 for z in caches[:5]:
  ref=z['reference_index'];ix=np.flatnonzero(mask[ref]);st={k:z[k] for k in ['G','xy','xm','ym','yv']};st=tr.subset_stats(st,ix);b,inter,lam,gap,nit=tr.path(st);pred=z['validation_mean'][:,ix]@b+inter;age=z['validation_age'];total+=np.abs(pred-age[:,None]).sum(axis=0);nval+=len(age);maxgap=max(maxgap,float(np.max(gap)));ns.append(len(ix))
 j=int(total.argmin());z=caches[5];ref=z['reference_index'];ix=np.flatnonzero(mask[ref]);st=tr.subset_stats({k:z[k] for k in ['G','xy','xm','ym','yv']},ix);b,inter,lam,gap,nit=tr.path(st);maxgap=max(maxgap,float(np.max(gap)));return {'reference_index':ref[ix],'coefficients':b[:,j],'intercept':float(inter[j]),'lambda_fraction':float(tr.FRACS[j]),'lambda_value':float(lam[j]),'CV_MAE':float(total[j]/nval),'CV_MAE_path':(total/nval).tolist(),'n_features_CV_folds':ns,'max_dual_gap':maxgap,'n_training_OneK_donors':int(z['n_training_donors'])}
def main():
 assert (W/'c_masks_complete.json').exists();ref=pd.read_csv(R/'R2_reference_genes.csv');iccold=pd.read_csv(R/'R2_S_gene_ICC.csv');types=os.environ.get('R4_C_TYPES','B,CD4T,CD8T,MONO,NK').split(',');ledger=pd.read_csv(O/'c_crossfit_mask_QA.csv')
 with threadpool_limits(limits=int(os.environ.get('R4_BLAS_THREADS','2'))):
  for ct in types:
   if (M/f'{ct}_complete.json').exists():print('REUSE_C_MODELS',ct,flush=True);continue
   caches=[dict(np.load(C/f'{ct}_{f}.npz')) for f in ['fold0','fold1','fold2','fold3','fold4','full']];masks=np.load(W/f'c_{ct}_ICC_masks.npz');pos=masks['reference_index'];values=masks['mask'];assert values.shape==(100,len(pos));start=time.monotonic()
   # Recover the known round3 c model using its old S-only mask. Compare only
   # fixed-model predictions, never holdout age outcomes for tuning.
   oldmask=np.zeros(len(ref),bool);s=iccold[iccold.celltype.eq(ct)].set_index('feature_id').ICC2_1.reindex(ref.feature_id).to_numpy();oldmask[s>=.7]=True;control=solve(caches,oldmask);z=caches[-1];ii=np.flatnonzero(oldmask[z['reference_index']]);pr=z['test_mean'][:,ii]@control['coefficients']+control['intercept'];old=np.load(R/'models'/f'{ct}_c.npz');assert np.array_equal(ref.feature_id.iloc[control['reference_index']].to_numpy(dtype=str),old['feature_ids']);oldpred=z['test_mean'][:,ii]@old['coefficients']+float(old['intercept']);err=float(np.max(abs(pr-oldpred)));assert err<1e-3,(ct,err);(O/f'c_{ct}_round3_recovery_QA.json').write_text(json.dumps({'status':'PASS','max_fixed_prediction_difference_years':err,'same_OneK_mask_training_protocol':True,'does_not_use_holdout_age_for_selection':True},indent=2));print('C_CONTROL',ct,err,flush=True)
   for rep in range(50):
    for held in [0,1]:
     tag=f'{ct}_r{rep:02d}_h{held}';path=M/f'{tag}.npz';receipt=M/f'{tag}.json'
     if path.exists() and receipt.exists():continue
     mask=np.zeros(len(ref),bool);mask[pos]=values[rep*2+held];result=solve(caches,mask);ix=result.pop('reference_index');beta=result.pop('coefficients');np.savez_compressed(path,feature_ids=ref.feature_id.iloc[ix].to_numpy(dtype=str),reference_index=ix,coefficients=beta,intercept=result['intercept'],lambda_fraction=result['lambda_fraction'],expression_unit='mean_native_pseudocell_log1p_CP10K');m=ledger[ledger.celltype.eq(ct)&ledger.repeat.eq(rep)&ledger.heldout_fold.eq(held)].iloc[0];result.update(celltype=ct,repeat=rep,heldout_fold=held,n_features=len(ix),n_nonzero=int((beta!=0).sum()),n_Allen_training_short_pairs=int(m.n_training_short_pairs),seed=20261025,Allen_heldout_used_for_training=False);receipt.write_text(json.dumps(result,indent=2))
    print('C_REPEAT',ct,rep,round(time.monotonic()-start,1),flush=True)
   (M/f'{ct}_complete.json').write_text(json.dumps({'status':'COMPLETE','models':100,'seconds':time.monotonic()-start},indent=2))
 if all((M/f'{ct}_complete.json').exists() for ct in ['B','CD4T','CD8T','MONO','NK']):
  rows=[json.loads(f.read_text()) for f in sorted(M.glob('*_r*_h*.json'))];cv=[]
  for z in rows:
   curve=z.pop('CV_MAE_path');cv.extend({'celltype':z['celltype'],'repeat':z['repeat'],'heldout_fold':z['heldout_fold'],'lambda_fraction':float(l),'CV_MAE':float(v)} for l,v in zip(tr.FRACS,curve));z['n_features_CV_folds']='|'.join(map(str,z['n_features_CV_folds']))
  assert len(rows)==500;pd.DataFrame(rows).to_csv(O/'c_500_model_training_QA.csv',index=False);pd.DataFrame(cv).to_csv(O/'c_OneK_CV_paths.csv',index=False);(M/'all_complete.json').write_text(json.dumps({'status':'COMPLETE','models':500,'secondary_only':True},indent=2));print('C_ALL_MODELS_COMPLETE',flush=True)
if __name__=='__main__':main()
