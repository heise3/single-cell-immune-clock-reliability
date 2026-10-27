#!/usr/bin/env python3
"""Full-gene OneK1K training inputs and S-only gene reliability; no E outcomes."""
from pathlib import Path
import json,hashlib,time,sys,gc
import numpy as np,pandas as pd,h5py,anndata as ad,cellxgene_census as c
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';A=P/'outputs/round2_A1';TYPES=['B','CD4T','CD8T','MONO','NK']
CFG={'soma.init_buffer_bytes':134217728,'py.init_buffer_bytes':134217728,'sm.mem.total_budget':17179869184,'sm.compute_concurrency_level':16,'sm.io_concurrency_level':32}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
class Aggregate:
 def __init__(self,obs,genes):
  self.groups=obs[['group','donor','gsm','celltype','age','split','CV_fold']].drop_duplicates('group').sort_values('group').reset_index(drop=True);self.lookup=pd.Index(self.groups.group);self.n=np.zeros(len(self.groups),np.int64);self.lib=np.zeros(len(self.groups));shape=(len(self.groups),len(genes));self.raw=np.zeros(shape);self.cp=np.zeros(shape);self.log=np.zeros(shape);self.expr=np.zeros(shape,dtype=np.uint32);self.genes=genes
 def add(self,x,obs):
  ix=self.lookup.get_indexer(obs.group);assert (ix>=0).all();unique,inv=np.unique(ix,return_inverse=True);G=sparse.csr_matrix((np.ones(len(ix)),(inv,np.arange(len(ix)))),shape=(len(unique),len(ix)));self.n+=np.bincount(ix,minlength=len(self.n));lib=obs.raw_sum.to_numpy(dtype=float);self.lib+=np.bincount(ix,weights=lib,minlength=len(self.n))
  def update(target,mat):
   q=(G@mat).tocoo();target[unique[q.row],q.col]+=q.data.astype(target.dtype,copy=False)
  update(self.raw,x);norm=x.astype(np.float64).multiply((10000/lib)[:,None]).tocsr();update(self.cp,norm);norm.data=np.log1p(norm.data);update(self.log,norm);binary=x.copy();binary.data=np.ones_like(binary.data);update(self.expr,binary)
 def save(self,name):
  valid=self.n>0;g=self.groups.copy();g['n_cells']=self.n;g['full_library_sum']=self.lib;g.to_csv(O/f'R2_{name}_group_QC.csv',index=False);g=g[valid].reset_index(drop=True);g.to_csv(W/f'{name}_groups.csv',index=False)
  with h5py.File(W/f'{name}_aggregates.h5','w') as h:
   for key,a in [('mean_CP10K',self.cp[valid]/self.n[valid,None]),('mean_logCP10K',self.log[valid]/self.n[valid,None]),('pseudobulk_logCPM',np.log2(1+self.raw[valid]*1e6/self.lib[valid,None])),('expressing_cells',self.expr[valid])]:h.create_dataset(key,data=a,compression='lzf',chunks=(min(64,len(a)),min(1024,a.shape[1])))
  print('AGGREGATES_SAVED',name,len(g),int(self.n.sum()),flush=True)
def prepare():
 assert (O/'R2_execution_contract.json').exists();v=pd.read_csv(O/'Census_global_var.csv');shared=set(pd.read_csv(A/'OneK1K_analysis_var.csv')['_index'])
 for f in A.glob('HIHA_Census_var_*.csv'):shared&=set(pd.read_csv(f)['_index'])
 genes=v[v.feature_id.isin(shared)].sort_values('soma_joinid').reset_index(drop=True);genes.to_csv(O/'R2_reference_genes.csv',index=False);print('SHARED_GENES',len(genes),flush=True)
 cross=json.loads((A/'G1_prime_execution_plan.json').read_text())['celltype_crosswalk'];onek=pd.read_parquet(P/'work/round2_A1/cache/OneK1K_census_obs.parquet');onek['celltype']=onek.cell_type.astype(str).map(cross);onek=onek[onek.celltype.notna()&onek.is_primary_data&onek.raw_sum.gt(0)].copy();onek=onek.rename(columns={'donor_id':'donor'});sp=pd.read_csv(O/'R2_OneK1K_donor_split.csv').rename(columns={'donor_id':'donor'});onek=onek.merge(sp[['donor','age','split','CV_fold']],on='donor',validate='many_to_one');onek['gsm']=onek.donor;onek['group']=onek.gsm+'|'+onek.celltype
 sdon=set(pd.read_csv(O/'R2_K0_donor_split.csv').loc[lambda d:d.split.eq('S'),'donor']);pairs=pd.read_csv(A/'provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.donor.isin(sdon)&pairs.stratum.eq('7-9d')];pairs.to_csv(O/'R2_S_short_pair_ledger.csv',index=False);gsms=set(pairs.gsm1)|set(pairs.gsm2);sobs=pd.read_parquet(P/'work/round2_A1/cache/HIHA_K0_clock_cells.parquet');sobs=sobs[sobs.gsm.isin(gsms)].copy().rename(columns={'clock_type':'celltype'});k0=pd.read_csv(A/'provenance_inputs/K0_primary_with_released_age.csv');sobs=sobs.merge(k0[['gsm','sample.subjectAgeAtDraw']].rename(columns={'sample.subjectAgeAtDraw':'age'}),on='gsm',validate='many_to_one');sobs['split']='S';sobs['CV_fold']=-1;sobs['group']=sobs.gsm+'|'+sobs.celltype
 for name,obs in [('OneK1K',onek),('S',sobs)]:obs.sort_values('soma_joinid').to_parquet(W/'cache'/f'{name}_stage1_cells.parquet',index=False)
 (O/'R2_stage1_plan.json').write_text(json.dumps({'OneK1K_cells_pre_native_QC':len(onek),'OneK1K_donors_pre_QC':onek.donor.nunique(),'S_7_9d_paired_donors':len(pairs),'S_cells':len(sobs),'common_source_genes':len(genes),'OneK1K_split_frozen':True,'E_outcomes_used':False},indent=2));return genes,v,onek,sobs
def run_population(name,obs,genes,globalvar,save_raw):
 out=W/f'{name}_aggregates.h5'
 if out.exists() and (W/f'{name}_aggregation_complete.json').exists():print('REUSE_COMPLETE',name,flush=True);return
 obs=obs.sort_values('soma_joinid').reset_index(drop=True);agg=Aggregate(obs,genes.feature_id.to_numpy());vid=genes.soma_joinid.to_numpy(np.int64);mito=globalvar.loc[globalvar.feature_name.str.startswith('MT-',na=False),'soma_joinid'].to_numpy(np.int64);D=W/'OneK1K_full_gene_counts';D.mkdir(exist_ok=True);man=[];excluded=0
 with c.open_soma(census_version='2025-11-08',tiledb_config=CFG) as f:
  X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw']
  for i,lo in enumerate(range(0,len(obs),40000)):
   rr=obs.iloc[lo:lo+40000].copy();file=D/f'OneK1K_full_{i:03d}.h5ad';qcfile=D/f'OneK1K_full_{i:03d}.QC.json';start=time.monotonic()
   if save_raw and file.exists() and qcfile.exists():a=ad.read_h5ad(file);xx=a.X;oo=a.obs.reset_index(drop=True);qc=json.loads(qcfile.read_text());excluded+=qc['excluded_cells'];assert list(a.var.feature_id)==list(genes.feature_id)
   else:
    chunks=[];metadata=[];expected=rr.set_index('soma_joinid');qcs=[]
    for mat,(rid,_) in X.read(coords=(rr.soma_joinid.to_numpy(np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
     rid=np.asarray(rid,np.int64);o=expected.loc[rid].reset_index();lib=np.asarray(mat.sum(axis=1)).ravel();assert np.array_equal(lib,o.raw_sum.to_numpy());nf=np.diff(mat.tocsr().indptr);mt=np.asarray(mat[:,mito].sum(axis=1)).ravel()/lib;keep=(nf>=200)&(nf<=5000)&(mt<=.25);o=o[keep].copy();o['native_n_features']=nf[keep];o['native_mito_fraction']=mt[keep];metadata.append(o);chunks.append(mat[keep][:,vid].tocsr());qcs.append({'input':len(rid),'kept':int(keep.sum())});print('RAW_BLOCK',name,i,len(rid),int(keep.sum()),flush=True)
    xx=sparse.vstack(chunks).tocsr();oo=pd.concat(metadata,ignore_index=True);qc={'input_cells':len(rr),'kept_cells':len(oo),'excluded_cells':len(rr)-len(oo),'full_row_sum_exact':True};excluded+=qc['excluded_cells']
    if save_raw:
     oo.index=oo.soma_joinid.astype(str);oo.index.name=None;vv=genes.copy();vv.index=vv.feature_id;vv.index.name=None;a=ad.AnnData(xx,obs=oo,var=vv);a.write_h5ad(file.with_suffix('.part.h5ad'),compression='lzf');file.with_suffix('.part.h5ad').rename(file);qcfile.write_text(json.dumps(qc,indent=2))
   agg.add(xx,oo)
   if save_raw:man.append({'path':str(file.relative_to(P)),'bytes':file.stat().st_size,'sha256':sha(file),**qc});(O/'R2_OneK1K_raw_manifest.json').write_text(json.dumps({'complete':False,'chunks':man},indent=2))
   print('CHUNK_DONE',name,i,len(oo),'seconds',round(time.monotonic()-start,2),flush=True);del xx,oo;gc.collect()
 agg.save(name);summary={'population':name,'input_cells':len(obs),'kept_cells':int(agg.n.sum()),'excluded_native_QC':excluded,'E_outcomes_used':False};(W/f'{name}_aggregation_complete.json').write_text(json.dumps(summary,indent=2));(O/f'R2_{name}_aggregation_summary.json').write_text(json.dumps(summary,indent=2))
 if save_raw:(O/'R2_OneK1K_raw_manifest.json').write_text(json.dumps({'complete':True,'chunks':man},indent=2))
def gene_icc():
 genes=pd.read_csv(O/'R2_reference_genes.csv');meta=pd.read_csv(W/'S_groups.csv');pairs=pd.read_csv(O/'R2_S_short_pair_ledger.csv');rows=[]
 with h5py.File(W/'S_aggregates.h5') as h:
  x=h['pseudobulk_logCPM'][:]
  for ct in TYPES:
   ix={g:i for i,g in enumerate(meta.group)};pp=pairs[pairs.gsm1.map(lambda g:g+'|'+ct in ix)&pairs.gsm2.map(lambda g:g+'|'+ct in ix)];a=np.array([[x[ix[g1+'|'+ct]],x[ix[g2+'|'+ct]]] for g1,g2 in zip(pp.gsm1,pp.gsm2)]);n=len(a);assert n>=3;gm=a.mean(axis=(0,1));rm=a.mean(axis=1);cm=a.mean(axis=0);msr=2*np.var(rm,axis=0,ddof=1);msc=n*np.var(cm,axis=0,ddof=1);mse=((a-rm[:,None,:]-cm[None,:,:]+gm)**2).sum(axis=(0,1))/(n-1);den=msr+mse+2/n*(msc-mse);icc=np.divide(msr-mse,den,out=np.full(len(genes),np.nan),where=den>1e-14);z=genes[['feature_id','feature_name']].copy();z['celltype']=ct;z['n_S_paired_donors']=n;z['ICC2_1']=icc;z['pass_0_5']=icc>=.5;z['pass_0_7']=icc>=.7;rows.append(z);print('S_GENE_ICC',ct,n,int((icc>=.5).sum()),int((icc>=.7).sum()),flush=True)
 pd.concat(rows,ignore_index=True).to_csv(O/'R2_S_gene_ICC.csv',index=False)
def main():
 genes,v,onek,s=prepare();run_population('OneK1K',onek,genes,v,True);run_population('S',s,genes,v,False);gene_icc();print('STAGE1_COMPLETE',flush=True)
if __name__=='__main__':main()
