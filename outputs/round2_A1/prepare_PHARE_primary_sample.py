#!/usr/bin/env python3
"""Real K0 specimen: exact same released cells, original full-gene GEO counts."""
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np,anndata as ad,h5py
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';k=pd.read_csv(P/'outputs/round2/K0_primary_with_released_age.csv').set_index('gsm');m=json.loads((O/'main_GEO_download_manifest.json').read_text());r=sorted([z for z in m['files'] if z['status']=='DOWNLOADED_VERIFIED' and z['gsm'] in k.index],key=lambda z:z['gsm'])[0];q=k.loc[r['gsm']]
meta=pd.read_parquet(W/'cache/HIHA_K0_cells.parquet');meta=meta[meta.gsm.eq(r['gsm'])].sort_values('soma_joinid');path=P/r['path']
with h5py.File(path,'r') as h:
 z=h['matrix'];shape=tuple(map(int,z['shape'][:]));x=sparse.csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=shape).T.tocsr();bar=z['barcodes'][:].astype(str);ix=pd.Index(bar).get_indexer(meta.barcodes.astype(str));assert (ix>=0).all() and len(np.unique(ix))==len(meta);x=x[ix].astype(np.float64);sym=z['features/name'][:].astype(str);ids=z['features/id'][:].astype(str);sample=np.unique(z['observations/pbmc_sample_id'][:].astype(str));assert list(sample)==[q['specimen.specimenGuid']]
 assert np.allclose(np.asarray(x.sum(axis=1)).ravel(),meta.n_umis.to_numpy()),'Source-GEO library totals must match the exact released cells'
a=ad.AnnData(x,obs=meta.set_index('soma_joinid').copy(),var=pd.DataFrame({'feature_name':sym,'feature_id':ids},index=ids));p=W/'PHARE_primary_full_gene_sample.h5ad';a.write_h5ad(p,compression='gzip')
report={'gsm':r['gsm'],'donor':q.donor,'specimen':q['specimen.specimenGuid'],'n_cells':a.n_obs,'n_features':a.n_vars,'input_path':str(p.relative_to(P)),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_GEO_path':str(path.relative_to(P)),'source_GEO_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'selection':'lexicographically first available verified cached K0 GSM, fixed before PHARE prediction; all exact final CZI released cells, source QC retained','identity_QA':'100% CZI barcodes matched one-to-one to original GEO; per-cell full original library totals match released n_umis','gene_universe':'Full original GEO feature universe for native PHARE/CellTypist, not the 2957-gene sc-ImmuAging+RAM projection; disclose small gene-universe difference from Census main input'};(O/'PHARE_input_manifest.json').write_text(json.dumps(report,indent=2));print('PREPARED',report,flush=True)
