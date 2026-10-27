#!/usr/bin/env python3
"""Reuse verified E; obtain remaining K0 inputs and extra short-pair genes."""
from pathlib import Path
import sys,json,os,time,gc
import numpy as np,pandas as pd,h5py,anndata as ad,cellxgene_census as census
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';W.mkdir(exist_ok=True);(W/'cache').mkdir(exist_ok=True);A=P/'outputs/round2_A1';R=P/'outputs/round3';V=P/'work/round3';sys.path.insert(0,str(R));import aggregate_training_inputs as ag;import prepare_evaluation_inputs as ev
ag.O=O;ag.W=W;ev.O=O;ev.W=W
def base():
 k=pd.read_csv(A/'provenance_inputs/K0_primary_with_released_age.csv');split=pd.read_csv(R/'R2_K0_donor_split.csv');s=set(split.loc[split.split.eq('S'),'donor']);obs=pd.read_parquet(P/'work/round2_A1/cache/HIHA_K0_clock_cells.parquet');obs=obs[obs.donor.isin(s)].rename(columns={'clock_type':'celltype'});obs=obs.merge(k[['gsm','sample.subjectAgeAtDraw']].rename(columns={'sample.subjectAgeAtDraw':'age'}),on='gsm',validate='many_to_one');obs['split']='full_K0_cache_S';obs['CV_fold']=-1;obs['group']=obs.gsm+'|'+obs.celltype;genes=pd.read_csv(R/'R2_evaluation_gene_universe.csv');return k,obs,genes
def S():
 k,obs,genes=base();v=pd.read_csv(R/'Census_global_var.csv');ev.run('K0_S',obs,genes,v)
def extra():
 if (W/'short_gene_inputs_complete.json').exists():print('REUSE_SHORT_GENES');return
 k,obs,genes=base();ref=pd.read_csv(R/'R2_reference_genes.csv');ranks=json.loads((R/'R2_training_gene_ranks.json').read_text());allref=sorted(set(i for ct in ranks.values() for z in ct.values() for i in z));allgenes=ref.iloc[allref].reset_index().rename(columns={'index':'reference_index'});allgenes.to_csv(O/'c_candidate_gene_universe.csv',index=False);new=allgenes[~allgenes.feature_id.isin(genes.feature_id)].copy();new.to_csv(O/'c_additional_E_short_genes.csv',index=False)
 pairs=pd.read_csv(A/'provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.stratum.eq('7-9d')];es=set(pd.read_csv(V/'E_groups.csv').gsm);ep=pairs[pairs.gsm1.isin(es)&pairs.gsm2.isin(es)];gsm=set(ep.gsm1)|set(ep.gsm2);meta=pd.read_csv(V/'E_groups.csv');master=pd.read_parquet(V/'cache/E_evaluation_master_obs.parquet');e=master.merge(meta[['group','donor','gsm','celltype','age']],on='group',validate='many_to_one');e=e[e.gsm.isin(gsm)].copy();e['split']='E_short_QC_reused';e['CV_fold']=-1;e=e.sort_values('soma_joinid').reset_index(drop=True);agg=ag.Aggregate(e,new.feature_id.to_numpy());vid=new.soma_joinid.to_numpy(np.int64);D=W/'E_short_additional_counts';D.mkdir(exist_ok=True);man=[]
 with census.open_soma(census_version='2025-11-08',tiledb_config=ev.CFG) as f:
  X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw']
  for i,lo in enumerate(range(0,len(e),40000)):
   rr=e.iloc[lo:lo+40000].copy();path=D/f'extra_{i:03d}.h5ad';st=time.monotonic()
   if path.exists():a=ad.read_h5ad(path);x=a.X;o=a.obs.reset_index(drop=True)
   else:
    xx=[];oo=[];lk=rr.set_index('soma_joinid')
    # Read contiguous full-feature rows: thousands of disjoint gene coordinate
    # ranges are much slower for this TileDB source. Project the identical1287
    # genes locally, additionally checking full library sums again.
    for mat,(rid,_) in X.read(coords=(rr.soma_joinid.to_numpy(np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
     rid=np.asarray(rid,np.int64);o=lk.loc[rid].reset_index();assert np.array_equal(np.asarray(mat.sum(axis=1)).ravel(),o.raw_sum.to_numpy());xx.append(mat[:,vid].tocsr());oo.append(o)
    x=sparse.vstack(xx).tocsr();o=pd.concat(oo,ignore_index=True);o.index=o.soma_joinid.astype(str);o.index.name=None;vv=new.copy();vv.index=vv.feature_id;vv.index.name=None;a=ad.AnnData(x,obs=o,var=vv);a.write_h5ad(path.with_suffix('.part.h5ad'),compression='lzf');path.with_suffix('.part.h5ad').rename(path)
   agg.add(x,o);man.append({'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':ag.sha(path),'n_cells':len(o),'QC_provenance':'exact fixed-version E post-native-QC source row IDs; full raw sums verified in round3, not inferred from this partial feature matrix'});(O/'E_short_additional_download_manifest.json').write_text(json.dumps({'complete':False,'chunks':man},indent=2));print('E_EXTRA',i,len(o),round(time.monotonic()-st,1),flush=True);del a,x,o;gc.collect()
 agg.save('E_short_extra');(O/'E_short_additional_download_manifest.json').write_text(json.dumps({'complete':True,'chunks':man},indent=2))
 # Assemble all33 paired donors, using exact full library denominators in every source.
 sm=pd.read_csv(V/'S_groups.csv');em=pd.read_csv(V/'E_groups.csv');exm=pd.read_csv(W/'E_short_extra_groups.csv');slk={f:i for i,f in enumerate(ref.feature_id)};elk={f:i for i,f in enumerate(genes.feature_id)};xlk={f:i for i,f in enumerate(new.feature_id)};si=[slk[f] for f in allgenes.feature_id];order=pd.concat([sm,em[em.gsm.isin(gsm)]],ignore_index=True).sort_values('group').reset_index(drop=True);assert len(order)==330 and order.donor.nunique()==33;out=np.zeros((len(order),len(allgenes)));lookup={g:i for i,g in enumerate(order.group)}
 with h5py.File(V/'S_aggregates.h5') as h:out[[lookup[g] for g in sm.group]]=h['pseudobulk_logCPM'][:][:,si]
 with h5py.File(V/'E_aggregates.h5') as h,h5py.File(W/'E_short_extra_aggregates.h5') as hx:
  erow={g:i for i,g in enumerate(em.group)};xrow={g:i for i,g in enumerate(exm.group)}
  for j,g in enumerate(order.group):
   if g not in erow:continue
   i=erow[g];ii=xrow[g];assert np.isclose(em.full_library_sum.iloc[i],exm.full_library_sum.iloc[ii]);u=h['pseudobulk_logCPM'][i];v=hx['pseudobulk_logCPM'][ii];out[j]=[u[elk[f]] if f in elk else v[xlk[f]] for f in allgenes.feature_id]
 with h5py.File(W/'K0_short_candidate_pseudobulk.h5','w') as h:h.create_dataset('pseudobulk_logCPM',data=out,compression='lzf')
 order.to_csv(W/'K0_short_groups.csv',index=False);(W/'short_gene_inputs_complete.json').write_text(json.dumps({'status':'COMPLETE','n_short_donors':33,'n_visits':66,'n_genes':len(allgenes),'additional_E_genes':len(new),'S_prior_aggregation':'14 paired donors,28 visits','E_aggregation':'19 paired donors,38 visits','full_libraries_preserved':True},indent=2));print('SHORT_GENES_COMPLETE',out.shape,flush=True)
if __name__=='__main__':
 assert (O/'round4_execution_contract.json').exists();{'S':S,'extra':extra}[sys.argv[1]]()
