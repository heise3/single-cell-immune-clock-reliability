#!/usr/bin/env python3
"""Engineering equivalence QA, not a reliability subsample."""
from pathlib import Path
import time,json
import numpy as np,pandas as pd,anndata as ad,cellxgene_census as c
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';t=time.monotonic();a=ad.read_h5ad(W/'OneK1K_counts/OneK1K_000.h5ad');a=a[:1000].copy();wanted=a.obs.soma_joinid.to_numpy(dtype=np.int64);vv=a.var.soma_joinid.to_numpy(dtype=np.int64);cfg={'soma.init_buffer_bytes':134217728,'py.init_buffer_bytes':134217728,'sm.mem.total_budget':17179869184,'sm.compute_concurrency_level':16,'sm.io_concurrency_level':32};rows=[]
with c.open_soma(census_version='2025-11-08',tiledb_config=cfg) as f:
 print('OPEN',time.monotonic()-t,flush=True);X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw'];ys=[];ids=[]
 for mat,(coord,_) in X.read(coords=(wanted,slice(None))).blockwise(axis=0,size=1000,reindex_disable_on_axis=1,eager=False).scipy():
  coord=np.asarray(coord);assert mat.shape[0]==len(coord) and mat.shape[1]>vv.max();ref=a.obs.set_index('soma_joinid').loc[coord];fullsum=np.asarray(mat.sum(axis=1)).ravel();assert np.allclose(fullsum,ref.raw_sum.to_numpy());ys.append(mat[:,vv].tocsr());ids.extend(coord.tolist());print('BLOCK',len(coord),mat.shape,mat.nnz,time.monotonic()-t,flush=True)
 y=sparse.vstack(ys).tocsr();ix=pd.Index(ids).get_indexer(wanted);assert (ix>=0).all();d=(y[ix]-a.X).tocsr();err=float(abs(d.data).max()) if d.nnz else 0.;assert err==0
 q={'status':'PASS_EXACT_RAW_COUNT_EQUIVALENCE','n_cells':len(wanted),'n_projected_features':len(vv),'maximum_absolute_count_difference':err,'full_row_UMI_totals_match_Census_raw_sum':True,'seconds':time.monotonic()-t,'config':cfg,'scope':'known cached OneK1K raw counts vs official blockwise full-row read; no scientific reliability analysis on this subset'};(O/'Census_blockwise_equivalence_QA.json').write_text(json.dumps(q,indent=2));print(q,flush=True)
