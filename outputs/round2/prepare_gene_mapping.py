#!/usr/bin/env python3
"""Prepare exact Ensembl + authoritative HGNC alias mapping without guessing."""
from pathlib import Path
import json,re,subprocess,hashlib,collections
import pandas as pd
import h5py
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parents[1];WORK=PROJECT/'work/round2';CACHE=WORK/'cache'
plan=json.loads((PROJECT/'outputs/K4_download_plan.json').read_text())
with h5py.File(PROJECT/'work/K4_data'/plan['files'][0]['filename']) as f:
 ids=f['matrix/features/id'][:].astype(str);symbols=f['matrix/features/name'][:].astype(str)
ref=pd.DataFrame({'geo_feature_id':ids,'geo_symbol':symbols});ref.to_csv(ROOT/'GEO_gene_reference.csv',index=False)
co=pd.read_csv(PROJECT/'outputs/scImmuAging_coefficients.csv');targets=sorted(set(co.loc[(co.coefficient!=0)&co.gene.ne('(Intercept)'),'gene']))
missing=sorted(set(targets)-set(ids));(CACHE/'ensembl_lookup_request.json').write_text(json.dumps({'ids':missing}))
HGNC=CACHE/'hgnc_complete_set.tsv'
if not HGNC.exists():subprocess.run(['curl','-fsSL','--max-time','60','-o',str(HGNC),'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'],check=True)
lookup=CACHE/'ensembl_lookup_missing.json'
if not lookup.exists() and not (CACHE/'ensembl_REST_attempted.txt').exists():
 p=subprocess.run(['curl','-fsSL','--max-time','60','-H','Content-Type: application/json','-H','Accept: application/json','-X','POST','--data-binary','@'+str(CACHE/'ensembl_lookup_request.json'),'https://rest.ensembl.org/lookup/id','-o',str(lookup)],capture_output=True,text=True)
 (CACHE/'ensembl_REST_attempted.txt').write_text(str(p.returncode)+' '+p.stderr)
 if p.returncode:print('Ensembl lookup access failure:',p.stderr,flush=True)
ens=json.loads(lookup.read_text()) if lookup.exists() else {}
h=pd.read_csv(HGNC,sep='\t',dtype=str,keep_default_na=False)
bm_path=CACHE/'biomart_remaining_response.tsv'
bm=pd.read_csv(bm_path,sep='\t',dtype=str,keep_default_na=False) if bm_path.exists() else pd.DataFrame()
bm_by_e={k:g for k,g in bm.groupby('Gene stable ID')} if len(bm) else {}
supp_path=ROOT/'NCBI_orgdb_mapping_supplement.csv'
supp=pd.read_csv(supp_path,dtype=str,keep_default_na=False) if supp_path.exists() else pd.DataFrame(columns=['ENSEMBL','SYMBOL'])
supp_symbols={k:set(g.SYMBOL) for k,g in supp.groupby('ENSEMBL')}
bidirectional_path=ROOT/'NCBI_orgdb_bidirectional_gene_identity.csv'
bi=pd.read_csv(bidirectional_path,dtype=str,keep_default_na=False)
entrez_by_e={k:set(g.ENTREZID)-{''} for k,g in bi.groupby('ENSEMBL')}
geo_by_entrez=collections.defaultdict(set)
for _,r in bi[bi.ENSEMBL.isin(set(ids))].iterrows():
 if r.ENTREZID:geo_by_entrez[r.ENTREZID].add(r.ENSEMBL)
h_by_e=collections.defaultdict(list);h_by_sym=collections.defaultdict(list)
for idx,row in h.iterrows():
 for e in row.ensembl_gene_id.split('|'):
  if e:h_by_e[e].append(idx)
 for name in [row.symbol]+row.alias_symbol.split('|')+row.prev_symbol.split('|'):
  if name:h_by_sym[name].append(idx)
historical_paths=sorted(CACHE.glob('hgnc_complete_set_*.txt'))
historical_bridge=[]
by_hgnc={r.hgnc_id:i for i,r in h.iterrows()}
for histpath in historical_paths:
 hist=pd.read_csv(histpath,sep='\t',dtype=str,keep_default_na=False)
 for _,r in hist.iterrows():
  if r.ensembl_gene_id and r.hgnc_id in by_hgnc:
   h_by_e[r.ensembl_gene_id].append(by_hgnc[r.hgnc_id])
   if r.ensembl_gene_id in targets:historical_bridge.append({'model_ensembl':r.ensembl_gene_id,'historical_snapshot':histpath.name,'hgnc_id':r.hgnc_id,'historical_symbol':r.symbol,'current_symbol':h.loc[by_hgnc[r.hgnc_id],'symbol']})
pd.DataFrame(historical_bridge).to_csv(ROOT/'historical_HGNC_gene_identity_bridge.csv',index=False)
sym_to_ids=collections.defaultdict(set)
for e,s in zip(ids,symbols):sym_to_ids[s].add(e)
rows=[]
for e in targets:
 if e in set(ids):rows.append({'model_ensembl':e,'mapped_geo_ids':e,'status':'EXACT_ENSEMBL','ensembl_symbol':'','hgnc_symbols':'','candidates':e});continue
 display=(ens.get(e) or {}).get('display_name','')
 bm_rows=bm_by_e.get(e,pd.DataFrame())
 approved_biomart=set(bm_rows['HGNC symbol'])-{''} if len(bm_rows) else set()
 if len(approved_biomart)==1:display=next(iter(approved_biomart))
 hits=set(h_by_e[e])
 if display:hits.update(h_by_sym.get(display,[]))
 # The locally version-pinned NCBI org.Hs.eg.db cross-reference can connect
 # an historical/alternative Ensembl ID to an HGNC-approved gene identity.
 for sym in supp_symbols.get(e,set()):hits.update(h_by_sym.get(sym,[]))
 candidates=set();syms=set();hgncsyms=[]
 # Prioritize a unique stable identity, then approved symbols, then aliases.
 # An alias shared by another gene must not override a canonical identity.
 canonical_candidates=set()
 for idx in hits:
  row=h.loc[idx]
  if row.ensembl_gene_id in set(ids):canonical_candidates.add(row.ensembl_gene_id)
  canonical_candidates.update(sym_to_ids.get(row.symbol,set()))
 entrez_candidates=set()
 for ent in entrez_by_e.get(e,set()):entrez_candidates.update(geo_by_entrez.get(ent,set()))
 if display:syms.add(display)
 for idx in hits:
  row=h.loc[idx];hgncsyms.append(row.symbol)
  syms.update([row.symbol]+row.alias_symbol.split('|')+row.prev_symbol.split('|'))
  if row.ensembl_gene_id in set(ids):candidates.add(row.ensembl_gene_id)
 for s in syms:candidates.update(sym_to_ids.get(s,set()))
 approved_hits=h.index[h.symbol.isin(approved_biomart)].tolist()
 approved_candidates=set()
 for idx in approved_hits:
  row=h.loc[idx]
  if row.ensembl_gene_id in set(ids):approved_candidates.add(row.ensembl_gene_id)
  approved_candidates.update(sym_to_ids.get(row.symbol,set()))
 if len(approved_biomart)==1 and len(approved_candidates)==1:
  candidates=approved_candidates;method='ENSEMBL_BIOMART_APPROVED_HGNC'
 elif len(approved_biomart)==1 and len(approved_candidates)>1 and all(sym in approved_biomart for eid,sym in zip(ids,symbols) if eid in approved_candidates):
  candidates=approved_candidates;method='SAME_APPROVED_GENE_MULTIPLE_GEO_FEATURES_SUM_COUNTS'
 elif len(canonical_candidates)==1:
  candidates=canonical_candidates;method='UNIQUE_CANONICAL_HGNC'
 elif len(entrez_by_e.get(e,set()))==1 and len(entrez_candidates)==1:
  candidates=entrez_candidates;method='UNIQUE_NCBI_ENTREZ_BRIDGE'
 else:method='UNIQUE_HGNC_OR_ENSEMBL_SYMBOL'
 # Only unique GEO-feature replacements are used. Keep every ambiguity.
 accepted=(len(candidates)==1 or method=='SAME_APPROVED_GENE_MULTIPLE_GEO_FEATURES_SUM_COUNTS')
 status=method if accepted else ('AMBIGUOUS' if candidates else 'UNMAPPED')
 rows.append({'model_ensembl':e,'mapped_geo_ids':'|'.join(sorted(candidates)) if accepted else '', 'status':status,'ensembl_symbol':display,'hgnc_symbols':'|'.join(sorted(set(hgncsyms))),'candidates':'|'.join(sorted(candidates))})
mapping=pd.DataFrame(rows);mapping.to_csv(ROOT/'scImmuAging_gene_mapping.csv',index=False)
used=mapping[mapping.mapped_geo_ids.ne('')]
collisions=used[used.duplicated('mapped_geo_ids',keep=False)]
collisions.to_csv(ROOT/'gene_mapping_many_model_features_to_one_GEO.csv',index=False)
covered=set(mapping.loc[mapping.mapped_geo_ids.ne(''),'model_ensembl']);summary=[]
for ct,g in co[(co.coefficient!=0)&co.gene.ne('(Intercept)')].groupby('celltype'):
 req=set(g.gene);remain=req-covered;summary.append({'clock':'scImmuAging','celltype':ct,'nonzero_genes':len(req),'missing_round1':len(req-set(ids)),'missing_after_Ensembl_HGNC':len(remain),'residual_missing_fraction':len(remain)/len(req)})
pd.DataFrame(summary).to_csv(ROOT/'G1_gene_coverage_GEO.csv',index=False)
receipts=[]
for p in [HGNC,lookup,CACHE/'ensembl_lookup_request.json',bidirectional_path,CACHE/'biomart_CD8_residual.tsv',bm_path]+historical_paths:
 if p.exists():receipts.append({'path':str(p.relative_to(PROJECT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(ROOT/'gene_annotation_receipts.json').write_text(json.dumps(receipts,indent=2))
print(pd.DataFrame(summary).to_string(index=False),flush=True)
