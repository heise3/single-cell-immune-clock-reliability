#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1/cache'
v=pd.read_csv(O/'HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv',keep_default_na=False);g=pd.read_csv(P/'outputs/round2/GEO_gene_reference.csv',keep_default_na=False);q=pd.read_csv(O/'RAM_fixed_gene_panels.csv');h=pd.read_csv(W/'hgnc_current_20261026.tsv',sep='\t',keep_default_na=False)
byname=defaultdict(set);approved=defaultdict(set);ref=defaultdict(set);geo=defaultdict(set)
for r in h.itertuples():
 approved[r.symbol].add(r.hgnc_id)
 for n in [r.symbol]+r.prev_symbol.split('|')+r.alias_symbol.split('|'):
  if n:byname[n].add(r.hgnc_id)
for e,n in zip(v['_index'],v.feature_name):ref[n].add(e)
for e,n in zip(g.geo_feature_id,g.geo_symbol):geo[n].add(e)
h=h.set_index('hgnc_id');vs=set(v['_index']);out=[]
for r in q.to_dict('records'):
 name=r['Gene'];ids=set(ref[name]);route='EXACT_SOURCE_SYMBOL'
 if not ids:
  ids=geo[name]&vs;route='EXACT_GEO_GENE_ID_BRIDGE'
 if not ids:
  candidates=approved[name] or byname[name]
  if len(candidates)==1:
   z=h.loc[next(iter(candidates))];ids=ref[z.symbol]|(set(z.ensembl_gene_id.split('|'))&vs);route='UNIQUE_HGNC_CURRENT_OR_HISTORICAL_SYMBOL'
 if not ids:route='UNMAPPED'
 r['mapped_feature_ids']='|'.join(sorted(ids));r['mapping_route']=route;out.append(r)
a=pd.DataFrame(out);a.to_csv(O/'RAM_mapped_gene_panels.csv',index=False)
print(a.groupby(['RAM_celltype','direction']).mapped_feature_ids.agg(['size',lambda x:x.ne('').sum()]).to_string());print(a[a.mapped_feature_ids.eq('')][['RAM_celltype','Gene']].to_string(index=False))
