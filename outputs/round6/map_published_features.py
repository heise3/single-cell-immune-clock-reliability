#!/usr/bin/env python3
"""A1 stable-ID identity chain to source symbols; never choose ambiguous aliases."""
from pathlib import Path
from collections import defaultdict
import json,pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';A=P/'outputs/round2_A1';CTS=['B','CD4T','CD8T','MONO','NK']
def main():
 assert json.loads((O/'all_count_inputs_ready.json').read_text())['status'].startswith('PASS');co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');ident=pd.read_csv(A/'model_gene_identity_recheck.csv').fillna('').set_index('model_ensembl');hg=pd.read_csv(P/'work/round2_A1/cache/hgnc_current_20261026.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('hgnc_id');approved=defaultdict(set)
 for h,r in hg.iterrows():approved[r.symbol.upper()].add(h)
 one=pd.read_csv(A/'OneK1K_var.csv',keep_default_na=False);frozen=dict(zip(one['_index'],one.feature_name));paper=pd.read_csv(A/'published_model_gene_annotation.csv',keep_default_na=False);paper_names=paper.groupby('model_ensembl').published_symbol.apply(lambda s:set(s)-{''}).to_dict();coverage=[];missing=[];allrows=[]
 for ct in CTS:
  var=pd.read_csv(O/f'source_{ct}_gene_annotation.csv');symidx=defaultdict(set)
  for r in var.itertuples():symidx[str(r.feature_name).upper()].add(int(r.feature_index))
  lookup=var.set_index('feature_index').feature_name.to_dict();g=co[co.celltype.eq(ct)&co.gene.ne('(Intercept)')&co.coefficient.ne(0)].copy().reset_index(drop=True);g['model_index']=np.arange(1,len(g)+1);g.to_csv(O/f'published_{ct}_weights.csv',index=False);inter=float(co.loc[co.celltype.eq(ct)&co.gene.eq('(Intercept)'),'coefficient'].iloc[0]);(O/f'published_{ct}_intercept.txt').write_text(repr(inter));rows=[];trip=[]
  for r in g.itertuples():
   e=r.gene;h=ident.loc[e,'hgnc_id'];cand=set();status='UNMAPPED_NATIVE_ZERO'
   if e.upper()in symidx:cand=symidx[e.upper()];status='EXACT_STABLE_ENSEMBL_ID'
   elif h in hg.index:
    rr=hg.loc[h];canon=symidx[rr.symbol.upper()]
    if canon:cand=canon;status='APPROVED_HGNC_SYMBOL'
    else:
     al=set().union(*(symidx[n.upper()] for n in (rr.prev_symbol+'|'+rr.alias_symbol).split('|') if n));al={i for i in al if not(approved[str(lookup[i]).upper()]-{h})}
     if len(al)==1:cand=al;status='UNIQUE_HGNC_HISTORICAL_ALIAS'
     elif len(al)>1:status='AMBIGUOUS_HGNC_ALIAS_NATIVE_ZERO'
   if not cand and status=='UNMAPPED_NATIVE_ZERO':
    names=set(paper_names.get(e,set()))|({frozen[e]} if e in frozen and frozen[e] else set());options=set().union(*(symidx[n.upper()] for n in names)) if names else set()
    if h:options={i for i in options if not(approved[str(lookup[i]).upper()]-{h})}
    if len(options)==1:cand=options;status='UNIQUE_FROZEN_STABLE_ID_TO_SYMBOL_ANNOTATION'
    elif len(options)>1:status='AMBIGUOUS_FROZEN_ANNOTATION_NATIVE_ZERO'
   rec={'celltype':ct,'gene':e,'coefficient':r.coefficient,'model_index':r.model_index,'hgnc_id':h,'display_symbol':ident.loc[e,'display_symbol'],'mapped_source_feature_indices':'|'.join(map(str,sorted(cand))),'mapped_source_symbols':'|'.join(str(lookup[i]) for i in sorted(cand)),'mapping_status':status};rows.append(rec)
   for i in sorted(cand):trip.append({'model_index':r.model_index,'feature_index':i,'source_symbol':lookup[i]})
   if not cand:missing.append(rec)
  allrows+=rows;pd.DataFrame(trip).to_csv(O/f'published_{ct}_source_mapping_triplets.csv',index=False);n=sum(bool(r['mapped_source_feature_indices']) for r in rows);coverage.append({'celltype':ct,'nonzero_model_genes':len(rows),'mapped_genes':n,'residual_missing_genes':len(rows)-n,'residual_missing_fraction':1-n/len(rows),'native_missing_policy':'zero only for unresolved/unmeasured sourcefeatures;no cache-column substitution','source_gene_annotations':'RNA RDS has gene symbols;stableEnsembl identity resolved via frozenHGNC and modelannotation before source-symbol match'})
 pd.DataFrame(allrows).to_csv(O/'published_gene_mapping_all_types.csv',index=False);pd.DataFrame(missing).to_csv(O/'residual_missing_gene_weights.csv',index=False);pd.DataFrame(coverage).to_csv(O/'model_gene_coverage.csv',index=False);print(pd.DataFrame(coverage).to_string(index=False),flush=True)
if __name__=='__main__':main()
