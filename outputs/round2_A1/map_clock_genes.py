#!/usr/bin/env python3
"""Build-aware, provenance-preserving feature matching; no outcome-based choices.

Exact stable IDs take precedence. A single HGNC identity can link multiple
features with the SAME approved symbol, which must be summed before log(CP10K).
Shared aliases never resolve a multi-gene ambiguity by choosing the first row.
"""
from pathlib import Path
from collections import defaultdict
import re,json,hashlib
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';OLD=P/'outputs/round2';W=P/'work/round2_A1/cache'
def clean(x):return re.sub(r'\.[0-9]+(_PAR_Y)?$',lambda z:z.group(1) or '',str(x))
def parts(x):return [z for z in str(x).split('|') if z]
def main():
 current=pd.read_csv(W/'hgnc_current_20261026.tsv',sep='\t',dtype=str,keep_default_na=False)
 by_hgnc={r.hgnc_id:r for r in current.itertuples()};by_name=defaultdict(set);by_approved=defaultdict(set);by_ens=defaultdict(set);by_entrez=defaultdict(set)
 for r in current.itertuples():
  by_approved[r.symbol.upper()].add(r.hgnc_id)
  for s in [r.symbol]+parts(r.prev_symbol)+parts(r.alias_symbol):by_name[s.upper()].add(r.hgnc_id)
  for e in parts(r.ensembl_gene_id):by_ens[clean(e)].add(r.hgnc_id)
  for e in parts(r.entrez_id):by_entrez[e].add(r.hgnc_id)
 histories=[W/'hgnc_complete_set_2020-07-01.tsv']+list((P/'work/round2/cache').glob('hgnc_complete_set_*.txt'))
 for fp in histories:
  if not fp.exists():continue
  h=pd.read_csv(fp,sep='\t',dtype=str,keep_default_na=False)
  for r in h.itertuples():
   if r.hgnc_id in by_hgnc and r.ensembl_gene_id:by_ens[clean(r.ensembl_gene_id)].add(r.hgnc_id)
 ncbi=pd.read_csv(OLD/'NCBI_orgdb_bidirectional_gene_identity.csv',dtype=str,keep_default_na=False)
 ncbi_by_e={k:g for k,g in ncbi.groupby('ENSEMBL')}
 bm=pd.read_csv(P/'work/round2/cache/biomart_remaining_response.tsv',sep='\t',dtype=str,keep_default_na=False)
 bm_by_e={k:g for k,g in bm.groupby('Gene stable ID')}
 paper=pd.read_csv(O/'published_model_gene_annotation.csv',keep_default_na=False)
 paper_by_e={e:set(g.published_symbol)-{''} for e,g in paper.groupby('model_ensembl')}
 geo=pd.read_csv(OLD/'GEO_gene_reference.csv',keep_default_na=False)
 geo_name={r.geo_feature_id:r.geo_symbol for r in geo.itertuples()}
 prior=pd.read_csv(OLD/'scImmuAging_gene_mapping.csv',keep_default_na=False).set_index('model_ensembl')
 coef=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');targets=sorted(set(coef.gene)-{'(Intercept)'})
 one=pd.read_csv(O/'OneK1K_var.csv',keep_default_na=False)
 one_name=dict(zip(one['_index'],one.feature_name))
 identity={};id_rows=[]
 for e in targets:
  strong=set(by_ens[e]);names=set(paper_by_e.get(e,set()))
  if e in geo_name:names.add(geo_name[e])
  if e in one_name:names.add(one_name[e])
  b=bm_by_e.get(e)
  if b is not None:
   approved=set(b['HGNC symbol'])-{''}
   for s in approved:strong.update(by_approved[s.upper()])
   names.update(set(b['Gene name'])-{''});names.update(approved)
  n=ncbi_by_e.get(e)
  if not strong and n is not None:
   ents=set(n.ENTREZID)-{''}
   if len(ents)==1:strong.update(by_entrez[next(iter(ents))])
  if not strong:
   for s in names:strong.update(by_approved[s.upper()])
  if not strong and n is not None:
   for s in set(n.SYMBOL)-{''}:strong.update(by_approved[s.upper()])
  if not strong:
   for s in names:strong.update(by_name[s.upper()])
  # Only a unique approved HGNC identity is accepted for identity-based fallback.
  h=next(iter(strong)) if len(strong)==1 else ''
  if h:
   r=by_hgnc[h];symbols={r.symbol}|set(parts(r.prev_symbol))|set(parts(r.alias_symbol));canonical=r.symbol
  else:symbols=set();canonical='|'.join(sorted(names)) or 'UNRESOLVED_ID'
  identity[e]={'hgnc_id':h,'canonical':canonical,'names':symbols,'candidate_hgnc_ids':'|'.join(sorted(strong)),'source_names':'|'.join(sorted(names))}
  id_rows.append({'model_ensembl':e,'hgnc_id':h,'display_symbol':canonical,'candidate_hgnc_ids':'|'.join(sorted(strong)),'source_names':'|'.join(sorted(names))})
 pd.DataFrame(id_rows).to_csv(O/'model_gene_identity_recheck.csv',index=False)
 def match(ref,label):
  ens=set(ref);symidx=defaultdict(set)
  for e,s in ref.items():symidx[s.upper()].add(e)
  out=[]
  for e in targets:
   info=identity[e];cand=set();status='UNMAPPED';accepted=False
   if e in ens:cand={e};status='EXACT_STABLE_ID';accepted=True
   else:
    h=info['hgnc_id']
    if h:
     r=by_hgnc[h];canon=set(symidx.get(r.symbol.upper(),set()))
     canon.update(k for k in parts(r.ensembl_gene_id) if k in ens)
     if len(canon)==1:cand=canon;status='UNIQUE_APPROVED_HGNC';accepted=True
     elif canon and all(ref[k].upper()==r.symbol.upper() for k in canon):cand=canon;status='SAME_APPROVED_GENE_SUM_RAW_COUNTS';accepted=True
     else:
      alias=set()
      for s in info['names']:alias.update(symidx.get(s.upper(),set()))
      # Reject aliases which are another gene's approved symbol.
      alias={k for k in alias if not(by_approved[ref[k].upper()]-{h})}
      if len(alias)==1:cand=alias;status='UNIQUE_HGNC_HISTORICAL_ALIAS';accepted=True
      elif canon or alias:cand=canon|alias;status='AMBIGUOUS'
    if not accepted and e in prior.index and prior.loc[e,'mapped_geo_ids']:
     prev=set(parts(prior.loc[e,'mapped_geo_ids']));pc=prev&ens
     if pc and label=='OneK1K':cand=pc;status='FROZEN_GEO_ID_BRIDGE';accepted=True
     elif pc==prev and pc and label in ['GEO','HIHA_Census']:cand=pc;status='FROZEN_GEO_MAPPING_ID_BRIDGE';accepted=True
   out.append({'model_ensembl':e,'display_symbol':info['canonical'],'hgnc_id':info['hgnc_id'],
    'mapped_feature_ids':'|'.join(sorted(cand)) if accepted else '', 'status':status,'all_candidates':'|'.join(sorted(cand))})
  frame=pd.DataFrame(out);frame.to_csv(O/f'gene_mapping_{label}.csv',index=False);return frame
 gm=match(geo_name,'GEO');om=match(one_name,'OneK1K');
 cvar=pd.read_csv(O/'HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv',keep_default_na=False);cm=match(dict(zip(cvar['_index'],cvar.feature_name)),'HIHA_Census');summary=[]
 for name,m in [('GEO',gm),('OneK1K',om),('HIHA_Census',cm)]:
  covered=set(m.loc[m.mapped_feature_ids.ne(''),'model_ensembl'])
  for ct,z in coef[coef.gene.ne('(Intercept)')].groupby('celltype'):
   nz=set(z.loc[z.coefficient.ne(0),'gene']);allg=set(z.gene)
   summary.append({'dataset':name,'celltype':ct,'nonzero_gene_count':len(nz),'nonzero_mapped':len(nz&covered),'nonzero_coverage':len(nz&covered)/len(nz),'all_model_input_count':len(allg),'all_model_input_mapped':len(allg&covered),'all_model_input_coverage':len(allg&covered)/len(allg),'G1_prime_coverage_pass':len(nz&covered)/len(nz)>=0.99})
 pd.DataFrame(summary).to_csv(O/'G1_prime_gene_coverage.csv',index=False)
 residual=coef[coef.coefficient.ne(0)&coef.gene.ne('(Intercept)')].merge(cm,left_on='gene',right_on='model_ensembl');residual=residual[residual.mapped_feature_ids.eq('')].drop(columns=['mapped_feature_ids']).merge(om[['model_ensembl','mapped_feature_ids','status']].rename(columns={'status':'validation_mapping_status'}),on='model_ensembl',suffixes=('','_validation'))
 residual['expression_mean']=float('nan');residual['expression_between_donor_SD']=float('nan');residual['expression_status']='PENDING_VALIDATION_EXPRESSION'
 residual.to_csv(O/'G1_prime_residual_genes_pending_expression.csv',index=False)
 old=set(pd.read_csv(OLD/'G1_residual_missing_gene_weights.csv').gene);new=set(residual.gene)
 (O/'gene_repair_recheck_summary.json').write_text(json.dumps({'additional_genes_recovered':sorted(old-new),'unresolved_or_absent_genes':sorted(new),'source_annotation_snapshots':[str(p.relative_to(P)) for p in histories],'native_coefficients_used':True,'published_coefficients_not_substituted':True},indent=2))
 print(pd.DataFrame(summary).to_string(index=False),flush=True)
if __name__=='__main__':main()
