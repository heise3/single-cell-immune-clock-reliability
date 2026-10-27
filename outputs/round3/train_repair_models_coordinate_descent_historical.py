#!/usr/bin/env python3
"""Leakage-safe donor CV; exact weighted sufficient-statistic LASSO paths."""
from pathlib import Path
import json,sys,subprocess,gc,time,warnings,os,fcntl
import numpy as np,pandas as pd,h5py
from sklearn.linear_model import enet_path
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';D=W/'training_design';M=W/'models';M.mkdir(exist_ok=True);SEED=20261025;FRACS=np.geomspace(1,1e-4,100);TYPES=['B','CD4T','CD8T','MONO','NK'];sys.path.insert(0,str(O));from prepare_training_design import rng_for
def sufficient(X,y,donor_age):
 age_counts=pd.Series(donor_age).value_counts();wd=np.array([1/age_counts[a] for a in donor_age]);w=np.repeat(wd/100,100) if len(X)==len(donor_age)*100 else wd.copy();assert len(w)==len(X);w=w/w.sum();xm=w@X;ym=float(w@y);dx=X-xm;dy=y-ym;wx=dx*np.sqrt(w[:,None]);G=wx.T@wx;xy=dx.T@(w*dy);yv=float(w@(dy*dy));return {'G':G,'xy':xy,'xm':xm,'ym':ym,'yv':yv}
def path(st):
 p=len(st['xy']);ym=st['ym']
 if not p:return np.zeros((0,100)),np.full(100,ym),np.zeros(100),np.zeros(100),np.zeros(100,dtype=int)
 maxlam=float(np.max(np.abs(st['xy'])));alphas=maxlam*FRACS
 if maxlam<1e-14:return np.zeros((p,100)),np.full(100,ym),alphas,np.zeros(100),np.zeros(100,dtype=int)
 if p==1:
  beta=np.sign(st['xy'][0])*np.maximum(abs(st['xy'][0])-alphas,0)/st['G'][0,0];beta=beta[None,:];return beta,ym-st['xm']@beta,alphas,np.zeros(100),np.ones(100,dtype=int)
 # The objective uses normalized weighted covariance/cross-product. Passing
 # n=1 and sqrt(weighted y variance) is an API carrier for these sufficient
 # statistics, not a one-subject training sample. Native glmnet QA is required.
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always',ConvergenceWarning)
  a,beta,gaps,nit=enet_path(np.zeros((1,p),order='F'),np.array([np.sqrt(st['yv'])]),l1_ratio=1.,alphas=alphas,precompute=np.ascontiguousarray(st['G']),Xy=np.ascontiguousarray(st['xy']),check_input=False,return_n_iter=True,tol=1e-8,max_iter=50000,selection='cyclic')
 if any(issubclass(x.category,ConvergenceWarning) for x in caught):
  print('NUMERICAL_RETRY',p,float(max(gaps)),flush=True)
  a,beta,gaps,nit=enet_path(np.zeros((1,p),order='F'),np.array([np.sqrt(st['yv'])]),l1_ratio=1.,alphas=alphas,precompute=np.ascontiguousarray(st['G']),Xy=np.ascontiguousarray(st['xy']),check_input=False,return_n_iter=True,tol=1e-8,max_iter=250000,selection='cyclic')
 assert np.isfinite(beta).all() and np.max(gaps)<max(1e-5,st['yv']*2e-8),('Unconverged path',p,float(max(gaps)),st['yv']);return beta,ym-st['xm']@beta,a,gaps,np.asarray(nit)
def subset_stats(st,ix):return {'G':st['G'][np.ix_(ix,ix)],'xy':st['xy'][ix],'xm':st['xm'][ix],'ym':st['ym'],'yv':st['yv']}
def native_validate(X,meta):
 if (O/'R2_native_Gram_validation.json').exists():assert json.loads((O/'R2_native_Gram_validation.json').read_text())['status']=='PASS';return
 di=np.flatnonzero(meta.split.eq('train80'))[:40];n=len(di);rr=(di[:,None]*100+np.arange(100)).ravel();z=np.asarray(X[np.ix_(rr,np.arange(min(100,X.shape[1])))]);ages=meta.age.iloc[di].to_numpy();y=np.repeat(ages,100);st=sufficient(z,y,ages);b,inter,aa,gaps,nit=path(st);inds=[10,30,60];pred=z@b[:,inds]+inter[inds];pd.DataFrame(z).to_csv(W/'native_Gram_X.tsv',sep='\t',index=False,header=False);cnt=pd.Series(ages).value_counts();ww=np.repeat(np.array([1/cnt[a] for a in ages]),100);pd.DataFrame({'y':y,'w':ww}).to_csv(W/'native_Gram_yw.tsv',sep='\t',index=False);pd.DataFrame({'lambda':aa[inds]}).to_csv(W/'native_Gram_lambda.tsv',sep='\t',index=False);pd.DataFrame(pred).to_csv(W/'native_Gram_python.tsv',sep='\t',index=False,header=False)
 subprocess.run(['Rscript',str(O/'validate_Gram_native.R')],cwd=P,check=True);native=pd.read_csv(W/'native_Gram_R.tsv',sep='\t',header=None).to_numpy();err=float(np.max(abs(native-pred)));assert err<1e-3,('Native mismatch',err);(O/'R2_native_Gram_validation.json').write_text(json.dumps({'status':'PASS','n_real_OneK_training_donors':n,'n_real_pseudocells':len(z),'n_genes':z.shape[1],'lambda_indices':inds,'max_abs_native_prediction_difference_years':err,'native_glmnet':'alpha1,standardizeFALSE,interceptTRUE,weights inverse age frequency, explicit matching lambda','sufficient_statistics_are_real_data':True},indent=2));print('NATIVE_GRAM_QA',err,flush=True)
def main():
 assert (D/'complete.json').exists() and (O/'R2_S_gene_ICC.csv').exists();rank=json.loads((O/'R2_training_gene_ranks.json').read_text());genes=pd.read_csv(O/'R2_reference_genes.csv');icctab=pd.read_csv(O/'R2_S_gene_ICC.csv');allmeta=pd.read_csv(W/'OneK1K_groups.csv');cvrows=[];modelrows=[];predrows=[];coefrows=[];solverrows=[]
 requested=os.environ.get('R3_TRAIN_TYPES',','.join(TYPES)).split(',');assert set(requested).issubset(TYPES)
 with h5py.File(W/'OneK1K_aggregates.h5') as h,threadpool_limits(limits=int(os.environ.get('R3_BLAS_THREADS','8'))):
  for ct in requested:
   done=M/f'{ct}_complete.json'
   if done.exists():print('REUSE_MODELS',ct,flush=True);continue
   started=time.monotonic();meta=pd.read_csv(D/f'{ct}_donor_metadata.csv');dg=pd.read_csv(D/f'{ct}_design_genes.csv');lookup={int(k):i for i,k in enumerate(dg.reference_index)};X=np.load(D/f'{ct}_pseudocells.npy',mmap_mode='r');means=np.load(D/f'{ct}_pseudocell_donor_means.npy');native_validate(X,meta) if ct=='B' and not (O/'R2_native_Gram_validation.json').exists() else None
   ai=allmeta.index[allmeta.celltype.eq(ct)].to_numpy();assert np.array_equal(allmeta.iloc[ai].group,meta.group);pb=h['pseudobulk_logCPM'][ai];ml=h['mean_logCP10K'][ai];age=meta.age.to_numpy();train=meta.split.eq('train80').to_numpy();test=~train;geneicc=icctab[icctab.celltype.eq(ct)].set_index('feature_id').ICC2_1.reindex(genes.feature_id).to_numpy();cache={};loss={k:np.zeros(100) for k in range(100,2001,100)};valn=0
   for fold in range(5):
    tr=train&meta.CV_fold.ne(fold).to_numpy();va=train&meta.CV_fold.eq(fold).to_numpy();ref=np.array(rank[ct][f'fold{fold}'],int);pos=np.array([lookup[k] for k in ref]);rr=np.flatnonzero(np.repeat(tr,100));xx=np.asarray(X[np.ix_(rr,pos)]);st=sufficient(xx,np.repeat(age[tr],100),age[tr]);del xx;cache[fold]={'st':st,'ref':ref,'tr':tr,'va':va};vm=means[va][:,pos];valn+=int(va.sum())
    unique={}
    for k in loss:
     kk=min(k,len(ref))
     if kk not in unique:
      bst=subset_stats(st,np.arange(kk));b,inter,aa,gaps,nit=path(bst);pred=vm[:,:kk]@b+inter;unique[kk]=np.abs(pred-age[va,None]).sum(axis=0);solverrows.append({'celltype':ct,'variant':'a','fold':fold,'K':k,'n_features':kk,'max_dual_gap':max(gaps),'max_iterations':max(nit)})
     loss[k]+=unique[kk]
    print('BASELINE_CV_FOLD',ct,fold,flush=True)
   candidates=[(float((vv/valn).min()),k,int((vv/valn).argmin())) for k,vv in loss.items()];bestloss,bestK,bestidx=min(candidates);cvrows.extend({'celltype':ct,'variant':'a','K':k,'lambda_fraction':FRACS[j],'CV_MAE':float(z/valn),'selected':k==bestK and j==bestidx} for k,vv in loss.items() for j,z in enumerate(vv));print('BASELINE_SELECTED',ct,bestK,bestidx,bestloss,flush=True)
   choices={'a':(bestidx,bestloss)}
   for variant in ['b','c','d','e']:
    total=np.zeros(100)
    for fold,z in cache.items():
     ref=z['ref'][:bestK];ix=np.arange(len(ref))
     if variant in ['b','c']:ix=ix[geneicc[ref]>=({'b':.5,'c':.7}[variant])];ref=ref[ix]
     tr=z['tr'];va=z['va']
     if variant=='e':bst=sufficient(pb[tr][:,ref],age[tr],age[tr]);vm=pb[va][:,ref];b,inter,aa,gaps,nit=path(bst)
     elif variant=='d':
      base=subset_stats(z['st'],ix);fit=ml[tr][:,ref];sd=fit.std(axis=0,ddof=1);valid=sd>1e-12;ix=ix[valid];ref=ref[valid];base=subset_stats(z['st'],ix);fit=fit[:,valid];sd=sd[valid];pc=PCA(n_components=min(50,len(fit)-1,len(ref)),svd_solver='randomized',random_state=int(rng_for('PCA',ct,fold).integers(0,2**31-1))).fit((fit-fit.mean(axis=0))/sd);A=pc.components_.T/sd[:,None];bst={'G':A.T@base['G']@A,'xy':A.T@base['xy'],'xm':np.zeros(A.shape[1]),'ym':base['ym'],'yv':base['yv']};bp,_,aa,gaps,nit=path(bst);b=A@bp;inter=base['ym']-base['xm']@b;vm=means[va][:,[lookup[k] for k in ref]]
     else:bst=subset_stats(z['st'],ix);b,inter,aa,gaps,nit=path(bst);vm=means[va][:,[lookup[k] for k in ref]]
     total+=abs(vm@b+inter-age[va,None]).sum(axis=0);solverrows.append({'celltype':ct,'variant':variant,'fold':fold,'K':bestK,'n_features':len(ref),'max_dual_gap':max(gaps),'max_iterations':max(nit)})
    j=int(total.argmin());choices[variant]=(j,float(total[j]/valn));cvrows.extend({'celltype':ct,'variant':variant,'K':bestK,'lambda_fraction':FRACS[k],'CV_MAE':float(z/valn),'selected':k==j} for k,z in enumerate(total));print('VARIANT_CV',ct,variant,choices[variant],flush=True)
   ref0=np.array(rank[ct]['full'][:bestK],int);pos0=np.array([lookup[k] for k in ref0]);rr=np.flatnonzero(np.repeat(train,100));xx=np.asarray(X[np.ix_(rr,pos0)]);fullst=sufficient(xx,np.repeat(age[train],100),age[train]);del xx
   for variant in ['a','b','c','d','e']:
    ref=ref0.copy();ix=np.arange(len(ref))
    if variant in ['b','c']:ix=ix[geneicc[ref]>=({'b':.5,'c':.7}[variant])];ref=ref[ix]
    if variant=='e':st=sufficient(pb[train][:,ref],age[train],age[train]);b,inter,aa,gaps,nit=path(st);vm=pb[test][:,ref];unit='log2(1+pooled_CPM)'
    elif variant=='d':
     fit=ml[train][:,ref];sd=fit.std(axis=0,ddof=1);valid=sd>1e-12;ix=ix[valid];ref=ref[valid];base=subset_stats(fullst,ix);fit=fit[:,valid];sd=sd[valid];pc=PCA(n_components=min(50,len(fit)-1,len(ref)),svd_solver='randomized',random_state=int(rng_for('PCA',ct,'full').integers(0,2**31-1))).fit((fit-fit.mean(axis=0))/sd);A=pc.components_.T/sd[:,None];st={'G':A.T@base['G']@A,'xy':A.T@base['xy'],'xm':np.zeros(A.shape[1]),'ym':base['ym'],'yv':base['yv']};bp,_,aa,gaps,nit=path(st);b=A@bp;inter=base['ym']-base['xm']@b;vm=means[test][:,[lookup[k] for k in ref]];unit='mean_native_pseudocell_log1p_CP10K'
    else:st=subset_stats(fullst,ix);b,inter,aa,gaps,nit=path(st);vm=means[test][:,[lookup[k] for k in ref]];unit='mean_native_pseudocell_log1p_CP10K'
    j,cvloss=choices[variant];beta=b[:,j];intercept=float(inter[j]);pred=vm@beta+intercept;featureids=genes.feature_id.iloc[ref].to_numpy(dtype=str);np.savez(M/f'{ct}_{variant}.npz',feature_ids=featureids,coefficients=beta,intercept=intercept,expression_unit=unit,lambda_value=aa[j],lambda_fraction=FRACS[j]);modelrows.append({'celltype':ct,'variant':variant,'baseline_selected_K':bestK,'n_features':len(ref),'n_nonzero':int((beta!=0).sum()),'CV_MAE':cvloss,'lambda_fraction':FRACS[j],'lambda_value':aa[j],'n_training_donors':int(train.sum()),'n_holdout_donors':int(test.sum()),'expression_unit':unit,'n_training_pseudocells':int(train.sum())*100 if variant!='e' else 0,'intercept':intercept});coefrows.extend({'celltype':ct,'variant':variant,'feature_id':f,'coefficient':float(co),'unit':unit} for f,co in zip(featureids,beta));predrows.extend({'dataset':'OneK1K_test20','celltype':ct,'variant':variant,'donor':don,'age':float(ag),'prediction':float(pr)} for don,ag,pr in zip(meta.loc[test,'donor'],age[test],pred));print('MODEL_SAVED',ct,variant,len(ref),'nonzero',int((beta!=0).sum()),flush=True)
   pd.DataFrame([z for z in modelrows if z['celltype']==ct]).to_csv(M/f'{ct}_model_results.csv',index=False);pd.DataFrame([z for z in cvrows if z['celltype']==ct]).to_csv(M/f'{ct}_CV_results.csv',index=False);pd.DataFrame([z for z in predrows if z['celltype']==ct]).to_csv(M/f'{ct}_holdout_predictions.csv',index=False);pd.DataFrame([z for z in coefrows if z['celltype']==ct]).to_csv(M/f'{ct}_coefficients.csv',index=False);pd.DataFrame([z for z in solverrows if z['celltype']==ct]).to_csv(M/f'{ct}_solver_QA.csv',index=False);done.write_text(json.dumps({'status':'COMPLETE','seconds':time.monotonic()-started,'seed':SEED},indent=2));print('CELLTYPE_COMPLETE',ct,round(time.monotonic()-started,1),flush=True);del X,pb,ml,cache;gc.collect()
 if not all((M/f'{ct}_complete.json').exists() for ct in TYPES):print('REQUESTED_TYPES_COMPLETE_OTHER_TYPES_PENDING',flush=True);return
 lock=open(M/'aggregate.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
 for suffix,name in [('model_results','R2_model_training_results.csv'),('CV_results','R2_CV_results.csv'),('holdout_predictions','R2_OneK1K_holdout_predictions.csv'),('coefficients','R2_new_model_coefficients.csv'),('solver_QA','R2_solver_QA.csv')]:pd.concat([pd.read_csv(M/f'{ct}_{suffix}.csv') for ct in TYPES],ignore_index=True).to_csv(O/name,index=False)
 (M/'all_models_complete.json').write_text(json.dumps({'status':'COMPLETE','seed':SEED,'models':25,'E_or_AIDA_used_for_training':False},indent=2));lock.close();print('ALL_MODELS_COMPLETE',flush=True)
if __name__=='__main__':main()
