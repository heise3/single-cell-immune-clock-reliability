import json,pathlib,time,collections
import h5py,fsspec,numpy as np
j=json.loads(pathlib.Path('work/cellxgene_collection.txt').read_text())
d=j['datasets'][0];u=d['assets'][0]['url'];print('opening',u,flush=True)
with fsspec.open(u,'rb',block_size=2**20,cache_type='blockcache') as stream:
 with h5py.File(stream,'r') as f:
  print('keys',list(f.keys()),flush=True);print('obs keys',list(f['obs'].keys()),flush=True)
  for k in f['obs']:
   o=f['obs'][k]
   if isinstance(o,h5py.Group):
    print(k,'categorical',len(o['categories']),o['categories'][:20],flush=True)
   else:print(k,o.shape,o.dtype,flush=True)
