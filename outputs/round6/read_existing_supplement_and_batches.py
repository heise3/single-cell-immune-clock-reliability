#!/usr/bin/env python3
"""Read the existing Allen r CSV; summarize released source batch identifiers."""
from pathlib import Path
import json,hashlib,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6';CTS=['B','CD4T','CD8T','MONO','NK']
def classify(a,b):
 a=set(a);b=set(b)
 if not a or not b:return 'unknown'
 if len(a)==len(b)==1 and a==b:return 'same_single_batch'
 if a.isdisjoint(b):return 'cross_disjoint_batches'
 return 'mixed_or_shared_multibatch'
def main():
 f=P/'outputs/round5/P2_chronological_age_signal.csv';z=pd.read_csv(f);z=z[z.model.isin(['published','e'])];w=z.pivot(index='celltype',columns='model',values='Pearson_r').loc[CTS].reset_index().rename(columns={'published':'published_Pearson_r','e':'e_Pearson_r'});w['n_donors']=68;w['source']='round5/P2_chronological_age_signal.csv;read_only_no_recalculation';w=w[['celltype','n_donors','published_Pearson_r','e_Pearson_r','source']];w.to_csv(O/'Allen_published_vs_e_Pearson_r.csv',index=False);(O/'Allen_supplement_read_receipt.json').write_text(json.dumps({'source':str(f.relative_to(P)),'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'new_Pearson_calculation':False,'n_rows':5},indent=2))
 if not (O/'K_c_status.json').exists() or json.loads((O/'K_c_status.json').read_text())['status']!='PASS':return
 obs=pd.read_parquet(W/'primary_selected_cells.parquet');pairs=pd.read_csv(O/'external_primary_pair_ledger.csv');rows=[]
 for ct in CTS:
  g=obs[obs.celltype.eq(ct)];bs=g.groupby('gsm').Batch.apply(lambda s:set(map(str,s)));fs=g.groupby('gsm').File_name.apply(lambda s:set(map(str,s)))
  for p in pairs.itertuples():rows.append({'celltype':ct,'donor':p.donor,'visit1':p.gsm1,'visit2':p.gsm2,'batch_set1':'|'.join(sorted(bs[p.gsm1])),'batch_set2':'|'.join(sorted(bs[p.gsm2])),'batch_relation':classify(bs[p.gsm1],bs[p.gsm2]),'technical_library_overlap':bool(fs[p.gsm1]&fs[p.gsm2]),'technical_library_set1':'|'.join(sorted(fs[p.gsm1])),'technical_library_set2':'|'.join(sorted(fs[p.gsm2]))})
 rows=pd.DataFrame(rows);rows.to_csv(O/'same_cross_batch_pair_ledger.csv',index=False);rows.groupby(['celltype','batch_relation']).size().unstack(fill_value=0).reset_index().to_csv(O/'same_cross_batch_pair_counts.csv',index=False)
if __name__=='__main__':main()
