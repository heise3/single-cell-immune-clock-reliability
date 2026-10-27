#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,importlib.metadata as md,sys
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';a=pd.read_parquet(W/'cache/HIHA_K0_cells.parquet');k=pd.read_csv(P/'outputs/round2/K0_primary_with_released_age.csv');assert a.gsm.nunique()==320 and a.donor.nunique()==68 and a.soma_joinid.is_unique
q=a.groupby(['donor','gsm'],observed=True).agg(n_cells=('soma_joinid','size'),source_UMIs=('n_umis','sum'),Census_UMIs=('raw_sum','sum')).reset_index();q['Census_over_source_UMI_ratio']=q.Census_UMIs/q.source_UMIs;q.to_csv(O/'K0_library_universe_QA.csv',index=False)
versions={'seed':20261025,'python':sys.version,'packages':{}}
for package in ['numpy','pandas','scipy','h5py','anndata','scanpy','statsmodels','cellxgene-census','tiledbsoma','pyarrow','fsspec','s3fs']:
 try:versions['packages'][package]=md.version(package)
 except md.PackageNotFoundError:versions['packages'][package]=None
versions['Census_release']='2025-11-08';versions['scImmuAging_commit']='112f8fab10c72aec6e7a7d83762a40117a1e6e24';versions['RAM_IHA_Figure_commit']='ad523b6a72ebf9394c1dcc80289e5cf7c6edb85a';versions['PHARE_commit']='0ab277b38ddac534d6d91e4580134cf6c78822ca';versions['siAge_commit']='3e8e17567afa296213b414083eaca59b30455b81';(O/'versions_and_seed.json').write_text(json.dumps(versions,indent=2));print('IDENTITY_VALIDATED',len(a),q.Census_over_source_UMI_ratio.describe().to_dict())
