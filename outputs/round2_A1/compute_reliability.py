#!/usr/bin/env python3
"""Donor-cluster CIs, leakage-safe LODO, and conditional variance components."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import linalg,optimize
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';SEED=20261025

def rng_for(*args):
 h=int.from_bytes(hashlib.sha256('|'.join(map(str,args)).encode()).digest()[:8],'little');return np.random.default_rng(np.random.SeedSequence([SEED,h]))

def icc_sem(a):
 a=np.asarray(a,dtype=float);n,k=a.shape
 if n<3 or k!=2 or not np.isfinite(a).all():return np.array([np.nan]*3)
 mean=a.mean();rm=a.mean(axis=1);cm=a.mean(axis=0);msr=k*np.var(rm,ddof=1);msc=n*np.var(cm,ddof=1);mse=np.sum((a-rm[:,None]-cm[None,:]+mean)**2)/((n-1)*(k-1))
 den=msr+(k-1)*mse+k/n*(msc-mse)
 icc=(msr-mse)/den if den>0 else np.nan
 sem=np.std(a.ravel(),ddof=1)*np.sqrt(1-icc) if np.isfinite(icc) and icc<=1 else np.nan
 return np.array([icc,sem,1.96*np.sqrt(2)*sem])

def lodo(g,y,multiplicity=None):
 # All visits of all copies of the held ORIGINAL donor are excluded.
 ids=sorted(g.donor.unique());di={d:i for i,d in enumerate(ids)};idx=g.donor.map(di).to_numpy();x=np.c_[np.ones(len(g)),g.age.to_numpy(dtype=float)];n=np.bincount(idx,minlength=len(ids));w=1/n[idx]
 xx=np.array([np.einsum('i,ij,ik->jk',w[idx==i],x[idx==i],x[idx==i]) for i in range(len(ids))]);xy=np.array([x[idx==i].T@(w[idx==i]*y[idx==i]) for i in range(len(ids))])
 m=np.ones(len(ids)) if multiplicity is None else np.asarray(multiplicity,dtype=float)
 total_x=np.einsum('i,ijk->jk',m,xx);total_y=np.einsum('i,ij->j',m,xy);coef=np.empty((len(ids),2))
 A=total_x[None,:,:]-m[:,None,None]*xx;bb=total_y[None,:]-m[:,None]*xy
 try:coef=np.linalg.solve(A,bb[...,None])[...,0]
 except np.linalg.LinAlgError:
  for i in range(len(ids)):coef[i]=np.linalg.lstsq(A[i],bb[i],rcond=None)[0]
 return y-np.sum(x*coef[idx],axis=1)

def batched_lodo(g,y,boot):
 ids=sorted(g.donor.unique());idx=pd.Categorical(g.donor,categories=ids).codes;x=np.c_[np.ones(len(g)),g.age.to_numpy(dtype=float)];n=np.bincount(idx);w=1/n[idx]
 xx=np.array([np.einsum('i,ij,ik->jk',w[idx==i],x[idx==i],x[idx==i]) for i in range(len(ids))]);xy=np.array([x[idx==i].T@(w[idx==i]*y[idx==i]) for i in range(len(ids))])
 tx=np.einsum('bd,dij->bij',boot,xx);ty=boot@xy;A=tx[:,None,:,:]-boot[:,:,None,None]*xx[None,:,:,:];b=ty[:,None,:]-boot[:,:,None]*xy[None,:,:]
 try:coef=np.linalg.solve(A,b[...,None])[...,0]
 except np.linalg.LinAlgError:coef=np.einsum('bdij,bdj->bdi',np.linalg.pinv(A),b)
 return y[None,:]-np.einsum('ni,bni->bn',x,coef[:,idx,:])

def pair_array(g,values,pairs):
 z=dict(zip(g.gsm,values));p=pairs[pairs.gsm1.isin(z)&pairs.gsm2.isin(z)]
 return p,np.array([[z[a],z[b]] for a,b in zip(p.gsm1,p.gsm2)],dtype=float).reshape(-1,2)

def reliability(g,pairs,clock,celltype,B=2000,ram=False):
 g=g.sort_values(['donor','gsm']).reset_index(drop=True);ids=sorted(g.donor.unique());raw=g.prediction.to_numpy(dtype=float);acc=lodo(g,raw);out=[]
 boot=rng_for('donor_bootstrap',clock,celltype).multinomial(len(ids),np.full(len(ids),1/len(ids)),size=B)
 boot_acc=batched_lodo(g,raw,boot)
 for metric,values in [('prediction',raw),('age_adjusted_RAM_residual' if ram else 'age_acceleration',acc)]:
  for stratum,pp in pairs.groupby('stratum',sort=True):
   q,a=pair_array(g,values,pp);point=icc_sem(a);samples=[]
   for ib,m in enumerate(boot):
    bv=values if metric=='prediction' else boot_acc[ib]
    qq,ba=pair_array(g,bv,pp);repeat=np.array([m[ids.index(d)] for d in qq.donor],dtype=int);ba=np.repeat(ba,repeat,axis=0);samples.append(icc_sem(ba))
   samples=np.array(samples);lo,hi=np.nanpercentile(samples,[2.5,97.5],axis=0)
   row={'clock':clock,'celltype':celltype,'metric':metric,'stratum':stratum,'n_donors':len(a),'ICC2_1':point[0],'ICC_CI_low':lo[0],'ICC_CI_high':hi[0],'bootstrap_replicates':B,'bootstrap_valid':np.isfinite(samples[:,0]).sum()}
   if not ram:row.update(SEM=point[1],SEM_CI_low=lo[1],SEM_CI_high=hi[1],MDC95=point[2],MDC95_CI_low=lo[2],MDC95_CI_high=hi[2])
   out.append(row)
 return out,acc

def variance_components(g,y):
 n=len(g);dn=g.donor.to_numpy();bn=g.batch.to_numpy();kd=(dn[:,None]==dn[None,:]).astype(float);kb=(bn[:,None]==bn[None,:]).astype(float);eye=np.eye(n);K=np.array([kd,kb,eye]);rank=np.linalg.matrix_rank(K.reshape(3,-1).T)
 known=g.cell_sampling_var.to_numpy(dtype=float)+g.algorithm_var.to_numpy(dtype=float);X=np.ones((n,1));scale=max(np.var(y,ddof=1),1e-6);yn=(y-y.mean())/np.sqrt(scale);known=known/scale
 def objective(v,ret=False):
  V=np.einsum('i,ijk->jk',v,K)+np.diag(known+1e-10)
  try:
   chol=linalg.cho_factor(V,lower=True,check_finite=False);iy=linalg.cho_solve(chol,yn,check_finite=False);ix=linalg.cho_solve(chol,X,check_finite=False);s=float((X.T@ix)[0,0]);res=yn-X[:,0]*float((X.T@iy)[0]/s);quad=float(res@linalg.cho_solve(chol,res,check_finite=False));val=2*np.log(np.diag(chol[0])).sum()+np.log(s)+quad
   return val
  except (linalg.LinAlgError,ValueError,FloatingPointError):return 1e100
 fit=optimize.minimize(objective,[.4,.1,.4],method='L-BFGS-B',bounds=[(0,None)]*3,options={'maxiter':400,'ftol':1e-10});vv=fit.x*scale;vc=float(g.cell_sampling_var.mean());va=float(g.algorithm_var.mean());total=vv.sum()+vc+va
 out={'n_visits':n,'n_donors':len(set(dn)),'n_batches':len(set(bn)),'covariance_design_rank':int(rank),'fit_converged':bool(fit.success),'optimizer_message':str(fit.message),'donor_variance':vv[0],'batch_variance':vv[1],'visit_nested_plus_unmodelled_residual_variance':vv[2],'cell_sampling_variance':vc,'pseudocell_algorithm_variance':va,'total_variance':total}
 for key in ['donor','batch','visit_nested_plus_unmodelled_residual','cell_sampling','pseudocell_algorithm']:out[key+'_fraction']=out[key+'_variance']/total if total>0 else np.nan
 out['interpretation']='Conditional variance model; visit and unmodelled technical residual cannot be separated with one measurement per visit. Boundary zero is not proof of no component.'
 return out

def main():
 data=pd.read_csv(O/'main_predictions.csv');pairs=pd.read_csv(P/'outputs/round2/K0_primary_pair_ledger.csv');rows=[];var=[];res=[]
 for (clock,ct),g in data.groupby(['clock','celltype'],sort=True):
  print('STATS',clock,ct,len(g),flush=True);rr,acc=reliability(g,pairs,clock,ct,ram=clock.startswith('RAM'));rows+=rr;g=g.sort_values(['donor','gsm']).reset_index(drop=True);g['age_adjusted']=acc;res.append(g)
  for metric,y in [('prediction',g.prediction.to_numpy()),('age_adjusted',acc)]:
   z=variance_components(g,y);z.update(clock=clock,celltype=ct,metric=metric);var.append(z)
  pd.DataFrame(rows).to_csv(O/'reliability_results.csv',index=False);pd.DataFrame(var).to_csv(O/'variance_decomposition.csv',index=False)
 pd.concat(res).to_csv(O/'main_predictions_with_LODO.csv',index=False)
 print('RELIABILITY_COMPLETE',flush=True)
if __name__=='__main__':main()
