#!/usr/bin/env python3
"""One frozen20-donor AIDA technical-replicate run; fixed pretrained models only."""
from pathlib import Path
from collections import defaultdict
import sys,json,gc,time
import numpy as np,pandas as pd,anndata as ad,cellxgene_census as census
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';W=P/'work/round5';A=P/'outputs/round2_A1';R=P/'outputs/round4';sys.path.insert(0,str(A));from run_G1_prime import native_indices,rng_for
from compute_reliability import reliability
CTS=['B','CD4T','CD8T','MONO','NK'];SEED=20261025
def mapping():
 co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');ids=pd.read_csv(A/'model_gene_identity_recheck.csv').fillna('').set_index('model_ensembl');hg=pd.read_csv(P/'work/round2_A1/cache/hgnc_current_20261026.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('hgnc_id');var=pd.read_csv(P/'outputs/round3/AIDA_source_var.csv');ref=dict(zip(var._index,var.feature_name));sym=defaultdict(set);approved=defaultdict(set)
 for e,n in ref.items():sym[str(n).upper()].add(e)
 for h,g in hg.iterrows():approved[g.symbol.upper()].add(h)
 rows=[]
 for e in sorted(set(co.gene)-{'(Intercept)'}):
  cand=set();status='UNMAPPED_NATIVE_ZERO'
  if e in ref:cand={e};status='EXACT_STABLE_ID'
  else:
   h=ids.loc[e,'hgnc_id']
   if h in hg.index:
    g=hg.loc[h];canon=sym[g.symbol.upper()]|{x for x in g.ensembl_gene_id.split('|') if x in ref}
    if len(canon)==1 or (canon and all(ref[x].upper()==g.symbol.upper() for x in canon)):cand=canon;status='UNIQUE_APPROVED_GENE_RAW_SUM'
    else:
     names=[g.symbol]+g.prev_symbol.split('|')+g.alias_symbol.split('|');al=set().union(*(sym[n.upper()] for n in names if n));al={x for x in al if not(approved[ref[x].upper()]-{h})}
     if len(al)==1:cand=al;status='UNIQUE_HISTORICAL_ALIAS'
  rows.append({'model_ensembl':e,'mapped_feature_ids':'|'.join(sorted(cand)),'status':status,'hgnc_id':ids.loc[e,'hgnc_id'],'display_symbol':ids.loc[e,'display_symbol']})
 m=pd.DataFrame(rows);m.to_csv(O/'P4_AIDA_published_gene_mapping.csv',index=False);md=m.set_index('model_ensembl').mapped_feature_ids.to_dict();cov=[]
 for ct,g in co[co.gene.ne('(Intercept)')&co.coefficient.ne(0)].groupby('celltype'):cov.append({'model':'published','celltype':ct,'nonzero_genes':len(g),'mapped_genes':sum(bool(md[e]) for e in g.gene),'coverage':sum(bool(md[e]) for e in g.gene)/len(g),'unmapped_model_genes':'|'.join(e for e in g.gene if not md[e]),'missing_policy':'native zero; measured features absent only from4044localcache fetched separately; this is source-specific coverage, not full model coverage assertion'})
 mods={ct:dict(np.load(R/'models'/f'{ct}_e.npz')) for ct in CTS}
 for ct,model in mods.items():
  nz=model['feature_ids'].astype(str)[model['coefficients']!=0];assert set(nz)<=set(ref);cov.append({'model':'e','celltype':ct,'nonzero_genes':len(nz),'mapped_genes':len(nz),'coverage':1.,'unmapped_model_genes':'','missing_policy':'all nonzero features available'})
 pd.DataFrame(cov).to_csv(O/'P4_model_feature_coverage.csv',index=False);return co,md,mods
def main():
 contract=json.loads((O/'P4_external_selection_contract.json').read_text());assert contract['n_selected_donors']==20;obs=pd.read_parquet(W/'P4_selected_cells.parquet').sort_values('soma_joinid');co,md,mods=mapping();man=json.loads((P/'outputs/round3/R2_AIDA_evaluation_raw_manifest.json').read_text());pieces=[]
 for row in man['chunks']:
  aa=ad.read_h5ad(P/row['path'],backed='r');idx=np.flatnonzero(aa.obs.soma_joinid.isin(obs.soma_joinid));
  if len(idx):pieces.append(aa[idx].to_memory())
  aa.file.close()
 genes=pieces[0].var.copy();assert all(list(q.var.feature_id)==list(genes.feature_id) for q in pieces);a=ad.concat(pieces,join='inner') if len(pieces)>1 else pieces[0];a.var=genes;del pieces;gc.collect();perm=pd.Series(np.arange(len(a)),index=a.obs.soma_joinid).loc[obs.soma_joinid].to_numpy();x=a.X[perm].tocsr();sourceobs=a.obs.iloc[perm].reset_index(drop=True);lib=sourceobs.raw_sum.to_numpy(float);assert np.array_equal(obs.soma_joinid,sourceobs.soma_joinid);geneids=list(a.var.feature_id.astype(str));fi={f:i for i,f in enumerate(geneids)};needed=sorted(set(f for g in co[co.gene.ne('(Intercept)')&co.coefficient.ne(0)].gene for f in md[g].split('|') if f)-set(geneids));print('P4_INPUT',len(obs),'cells',len(needed),'additional measured genes',flush=True)
 if needed:
  file=W/'P4_additional_measured_counts.h5ad'
  if file.exists():extra=ad.read_h5ad(file);assert list(extra.obs.soma_joinid)==list(obs.soma_joinid) and list(extra.var.feature_id)==needed
  else:
   gv=pd.read_csv(P/'outputs/round3/Census_global_var.csv').set_index('feature_id');vid=gv.loc[needed,'soma_joinid'].to_numpy(np.int64);cfg={'sm.compute_concurrency_level':4,'sm.io_concurrency_level':8};xs=[];ridall=[]
   with census.open_soma(census_version='2025-11-08',tiledb_config=cfg) as f:
    print('P4_CENSUS_OPENED',flush=True)
    X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw']
    # Read complete sparse rows then project: thousands of separate column
    # ranges form a costly Cartesian range query in TileDB. No new cell/QC
    # selection; only the1145 additional measured genes are retained on disk.
    for xx,(rid,_) in X.read(coords=(obs.soma_joinid.to_numpy(np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
     xs.append(xx[:,vid].tocsr());ridall.extend(rid);print('P4_EXTRA_BLOCK',len(ridall),flush=True)
   ex=sparse.vstack(xs).tocsr();ix=pd.Series(np.arange(len(ridall)),index=ridall).loc[obs.soma_joinid].to_numpy();extra=ad.AnnData(ex[ix],obs=pd.DataFrame({'soma_joinid':obs.soma_joinid.to_numpy()},index=obs.soma_joinid.astype(str)),var=pd.DataFrame({'feature_id':needed},index=needed));extra.obs.index.name=None;extra.var.index.name=None;extra.write_h5ad(file,compression='lzf')
  x=sparse.hstack([x,extra.X],format='csr');geneids+=needed;fi={f:i for i,f in enumerate(geneids)}
 rows=[]
 for ct in CTS:
  ii=np.flatnonzero(obs.celltype.eq(ct));o=obs.iloc[ii].copy();c=co[co.celltype.eq(ct)&co.gene.ne('(Intercept)')&co.coefficient.ne(0)];mg=c.gene.tolist();rr=[];cc=[]
  for j,g in enumerate(mg):
   for f in md[g].split('|'):
    if f:rr.append(fi[f]);cc.append(j)
  transform=sparse.csc_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(fi),len(mg)));cp=(x[ii]@transform).astype(float).multiply((10000/lib[ii])[:,None]).tocsr();cp.data=np.log1p(cp.data);score=np.asarray(cp@c.coefficient.to_numpy()).ravel()+float(co.loc[co.celltype.eq(ct)&co.gene.eq('(Intercept)'),'coefficient'].iloc[0]);model=mods[ct];beta=model['coefficients'];nz=beta!=0;xx=x[ii][:,[fi[f] for f in model['feature_ids'].astype(str)[nz]]];groups=o.groupby('gsm',sort=True).indices
  for gsm,jj in groups.items():
   cellscore=score[jj];ix=native_indices(len(jj),rng_for('R5_AIDA_external',gsm,ct));pub=float(cellscore[ix].mean());ee=float(np.log2(1+np.asarray(xx[jj].sum(axis=0)).ravel()*1e6/lib[ii[jj]].sum())@beta[nz]+float(model['intercept']));first=o.iloc[jj[0]];base={'gsm':gsm,'donor':first.donor,'age':first.age,'celltype':ct,'replicate':int(first.replicate),'library_id':first.library_id,'n_cells':len(jj),'repeat_type':'same_draw_technical_GEX_library','interval_days':0}
   rows.extend([{**base,'model':'published','prediction':pub},{**base,'model':'e','prediction':ee}])
 pred=pd.DataFrame(rows);assert len(pred)==400;pred.to_csv(O/'P4_external_technical_predictions.csv',index=False);pair=pred[['donor','gsm','replicate']].drop_duplicates().pivot(index='donor',columns='replicate',values='gsm').reset_index().rename(columns={1:'gsm1',2:'gsm2'});pair['stratum']='same-draw_technical0d';pair.to_csv(O/'P4_external_technical_pair_ledger.csv',index=False);out=[];res=[]
 for (ct,model),g in pred.groupby(['celltype','model']):
  rr,acc=reliability(g,pair,model,ct,B=2000);out.extend(rr);gg=g.sort_values(['donor','gsm']).copy();gg['age_acceleration']=acc;res.append(gg);print('P4_ICC',ct,model,flush=True)
 pd.DataFrame(out).to_csv(O/'P4_external_technical_ICC.csv',index=False);pd.concat(res).to_csv(O/'P4_external_technical_predictions_with_LODO.csv',index=False);(O/'P4_external_run_status.json').write_text(json.dumps({'status':'PASS_TECHNICAL_REPLICATION_EXECUTED','n_donors':20,'n_library_units':40,'n_cells':len(obs),'n_models':2,'n_celltypes':5,'ICC_rows':20,'seed':SEED,'biological_longitudinal_reproduction':'NOT_RUN; technical experiment cannot establish90day or1year reliability','AIDA_not_new_training_or_model_selection':True,'source_feature_coverage_file':'P4_model_feature_coverage.csv','full_library_denominator_from_prior_verified_raw_counts':True},indent=2))
if __name__=='__main__':main()
