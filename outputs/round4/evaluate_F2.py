#!/usr/bin/env python3
"""Primary e: full68 fixed-model paired difference bootstrap, untouched by c/T."""
from pathlib import Path
import sys,json,datetime
import numpy as np,pandas as pd,h5py
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';R=P/'outputs/round3';V=P/'work/round3';A=P/'outputs/round2_A1';sys.path.insert(0,str(A));from compute_reliability import lodo,batched_lodo,icc_sem,pair_array
def assemble():
 assert (W/'K0_S_evaluation_input_complete.json').exists();genes=pd.read_csv(R/'R2_evaluation_gene_universe.csv');m1=pd.read_csv(W/'K0_S_groups.csv');m2=pd.read_csv(V/'E_groups.csv');m=pd.concat([m1,m2],ignore_index=True).sort_values('group').reset_index(drop=True);assert m.donor.nunique()==68 and m.gsm.nunique()==320 and len(m)==1600;lookup={g:i for i,g in enumerate(m.group)};m.to_csv(W/'K0_full_groups.csv',index=False)
 with h5py.File(W/'K0_full_inputs.h5','w') as h:
  for name in ['pseudocell_mean_logCP10K','pseudobulk_logCPM']:
   out=np.zeros((len(m),len(genes)))
   for prefix,mm in [(W/'K0_S_aggregates.h5',m1),(V/'E_aggregates.h5',m2)]:
    with h5py.File(prefix) as src:out[[lookup[g] for g in mm.group]]=src[name][:]
   h.create_dataset(name,data=out,compression='lzf')
 genes.to_csv(O/'K0_prediction_gene_universe.csv',index=False);return m,genes
def predict(m,genes):
 lk={f:i for i,f in enumerate(genes.feature_id)};rows=[]
 with h5py.File(W/'K0_full_inputs.h5') as h:
  ps=h['pseudocell_mean_logCP10K'][:];pb=h['pseudobulk_logCPM'][:]
  for ct in ['B','CD4T','CD8T','MONO','NK']:
   take=m.celltype.eq(ct).to_numpy();g=m[take]
   for var in ['a','e']:
    model=np.load(R/'models'/f'{ct}_{var}.npz');ix=[lk[f] for f in model['feature_ids'].astype(str)];x=(pb if var=='e' else ps)[take][:,ix];p=x@model['coefficients']+float(model['intercept']);rows.extend({'dataset':'Allen_K0_full','celltype':ct,'variant':var,'gsm':gsm,'donor':don,'age':float(age),'prediction':float(y)} for gsm,don,age,y in zip(g.gsm,g.donor,g.age,p))
 out=pd.DataFrame(rows);old=pd.read_csv(R/'R2_all_predictions.csv');old=old[old.dataset.eq('Allen_E')&old.variant.isin(['a','e'])];q=out.merge(old[['celltype','variant','gsm','prediction']],on=['celltype','variant','gsm'],suffixes=['','_round3'],validate='one_to_one');err=float(np.max(abs(q.prediction-q.prediction_round3)));assert len(q)==1650 and err<1e-10;(O/'F2_round3_E_prediction_identity_QA.json').write_text(json.dumps({'status':'PASS','matched_predictions':len(q),'max_absolute_difference':err,'fixed_a_e_models_unchanged':True},indent=2));out.to_csv(O/'F2_a_e_full68_predictions.csv',index=False);return out
def main():
 m,genes=assemble();pred=predict(m,genes);pairs=pd.read_csv(A/'provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.stratum.isin(['75-105d','330-400d'])];ids=sorted(pred.donor.unique());boots=np.random.default_rng(np.random.SeedSequence([20261025,4202])).multinomial(68,np.full(68,1/68),size=2000);rows=[];bs={};allres=[]
 for (ct,var),g in pred.groupby(['celltype','variant']):
  g=g.sort_values(['donor','gsm']).reset_index(drop=True);y=g.prediction.to_numpy();acc=lodo(g,y);sample=batched_lodo(g,y,boots);g['age_acceleration']=acc;allres.append(g)
  for st,pp in pairs.groupby('stratum'):
   q,a=pair_array(g,acc,pp);pt=icc_sem(a)[0];vals=[]
   for i,multi in enumerate(boots):
    qq,b=pair_array(g,sample[i],pp);rep=np.array([multi[ids.index(d)] for d in qq.donor]);b=np.repeat(b,rep,axis=0);vals.append(icc_sem(b)[0])
   v=np.array(vals);lo,hi=np.nanpercentile(v,[2.5,97.5]);assert np.isfinite(v).all();rows.append({'celltype':ct,'variant':var,'stratum':st,'n_donors_full':68,'n_pairs':len(q),'ICC2_1':pt,'ICC_CI_low':lo,'ICC_CI_high':hi,'bootstrap_replicates':2000,'bootstrap_valid':int(np.isfinite(v).sum()),'metric':'LODO_age_acceleration','CI_scope':'Fixed trained model; full68 donor-cluster bootstrap, refitLODO; a/e share every draw'});bs[ct,var,st]=v
  print('F2_FULL_ICC',ct,var,flush=True)
 r=pd.DataFrame(rows);r.to_csv(O/'Figure4_F2_full68_a_e_ICC.csv',index=False);delta=[]
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  for st in ['75-105d','330-400d']:
   a=r[r.celltype.eq(ct)&r.variant.eq('a')&r.stratum.eq(st)].iloc[0];e=r[r.celltype.eq(ct)&r.variant.eq('e')&r.stratum.eq(st)].iloc[0];d=bs[ct,'e',st]-bs[ct,'a',st];lo,hi=np.percentile(d,[2.5,97.5]);delta.append({'celltype':ct,'stratum':st,'n_pairs':int(a.n_pairs),'baseline_a_ICC':a.ICC2_1,'variant_e_ICC':e.ICC2_1,'delta_ICC_e_minus_a':e.ICC2_1-a.ICC2_1,'delta_CI_low':lo,'delta_CI_high':hi,'paired_bootstrap_replicates':2000,'paired_bootstrap_valid':int(np.isfinite(d).sum()),'F2_target_type':ct in ['B','MONO','NK'],'CI_excludes_zero':bool(lo>0 or hi<0),'CI_confirms_improvement':bool(lo>0),'CI_confirms_worsening':bool(hi<0),'multiplicity':'pointwise95%CI; six F2 target comparisons, no familywise adjustment prescribed'})
 d=pd.DataFrame(delta);d.to_csv(O/'Figure4_F2_paired_delta_ICC.csv',index=False);passed=d[d.F2_target_type&d.CI_excludes_zero];positive=d[d.F2_target_type&d.CI_confirms_improvement];label='确认' if len(passed) else '探索性';decision={'status':'PASS' if len(passed) else 'FAIL','F2_result':label,'primary_variant':'e','secondary_c_does_not_decide_F2':True,'n_target_CIs_excluding_zero':len(passed),'n_target_positive_CIs':len(positive),'qualifying_comparisons':passed.to_dict('records'),'gate_unchanged':True,'scope':'Internal full-cohort extension, includes previously evaluated E; pointwise CI, not an independent replication or global multiple-testing guarantee','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(O/'F2_decision.json').write_text(json.dumps(decision,indent=2,ensure_ascii=False));pd.concat(allres,ignore_index=True).to_csv(O/'F2_a_e_predictions_with_LODO.csv',index=False);np.savez_compressed(W/'F2_paired_ICC_bootstrap.npz',**{'__'.join(k):v for k,v in bs.items()});print(json.dumps(decision,indent=2,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
