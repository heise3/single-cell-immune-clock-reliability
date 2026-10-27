import json,pathlib,concurrent.futures,csv,collections,time
import h5py,fsspec,numpy as np
root=pathlib.Path('work');j=json.loads((root/'cellxgene_collection.txt').read_text())
cols=['subject.subjectGuid','sample.visitName','specimen.specimenGuid','batch_id','subject.cmv','sex']
def dec(x):return x.decode() if isinstance(x,bytes) else x.item() if isinstance(x,np.generic) else x
def run(d):
 out=root/('k0_obs_'+d['dataset_id']+'.json')
 if out.exists():return json.loads(out.read_text())
 u=d['assets'][0]['url'];print('start',d['title'],flush=True)
 with fsspec.open(u,'rb',block_size=2**20,cache_type='blockcache') as stream:
  with h5py.File(stream,'r') as f:
   ob=f['obs'];sc=ob['sample.sampleKitGuid']['codes'][:];cats=ob['sample.sampleKitGuid']['categories'][:];print('sample_codes',d['title'],flush=True)
   unique,idx,counts=np.unique(sc,return_index=True,return_counts=True)
   rr=[{'sample.sampleKitGuid':dec(cats[c]),'n_cells':int(n),'dataset_id':d['dataset_id'],'dataset_version_id':d['dataset_version_id']} for c,n in zip(unique,counts)]
   for col in cols:
    o=ob[col]
    if isinstance(o,h5py.Group):
     codes=o['codes'][:];categories=o['categories'][:];values=categories[codes[idx]]
    else:values=o[:][idx]
    for r,val in zip(rr,values):r[col]=dec(val)
   result={'samples':rr,'X_dtype':str(f['X/data'].dtype),'X_shape':[int(x) for x in f['X'].attrs['shape']],'layers':list(f['layers'].keys()),'url':u}
 out.write_text(json.dumps(result,indent=2));print('done',d['title'],len(rr),sum(r['n_cells'] for r in rr),flush=True);return result
if __name__ == '__main__':
 import multiprocessing
 with concurrent.futures.ProcessPoolExecutor(4,mp_context=multiprocessing.get_context('spawn')) as ex:results=list(ex.map(run,j['datasets']))
 rr=[r for d in results for r in d['samples']]
 with (root/'k0_sample_cells.csv').open('w') as f:
  w=csv.DictWriter(f,list(rr[0]));w.writeheader();w.writerows(rr)
 n=np.array([r['n_cells'] for r in rr]);print('SUMMARY',json.dumps({'samples':len(rr),'cells':int(n.sum()),'quantiles':dict(zip(['min','q1','median','q3','max'],np.quantile(n,[0,.25,.5,.75,1]).tolist()))}),flush=True)
