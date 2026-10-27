#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4'
def main():
 present=set(pd.read_csv(P/'outputs/round3/AIDA_source_var.csv')['_index']);rows=[]
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  for var in ['a','c','e']:
   files=list((P/'work/round4/c_models').glob(ct+'_r*_h*.npz')) if var=='c' else [P/'outputs/round3/models'/f'{ct}_{var}.npz'];nonzero=set();missing=[]
   for f in files:
    z=np.load(f);ids=z['feature_ids'].astype(str);b=z['coefficients'];nonzero.update(ids[b!=0]);missing.extend({'model':f.name,'feature_id':str(i),'coefficient':float(co)} for i,co in zip(ids,b) if co!=0 and i not in present)
   rows.append({'celltype':ct,'variant':var,'n_models':len(files),'n_unique_nonzero_genes':len(nonzero),'n_annotation_absent_unique_genes':len(nonzero-present),'nonzero_gene_union_coverage':1-len(nonzero-present)/len(nonzero),'missing_genes':'|'.join(sorted(nonzero-present)),'missing_model_feature_occurrences':len(missing)})
   if missing:pd.DataFrame(missing).to_csv(O/f'AIDA_missing_{ct}_{var}_features.csv',index=False)
 pd.DataFrame(rows).to_csv(O/'AIDA_a_c_e_feature_coverage.csv',index=False);print('AIDA_COVERAGE_COMPLETE')
if __name__=='__main__':main()
