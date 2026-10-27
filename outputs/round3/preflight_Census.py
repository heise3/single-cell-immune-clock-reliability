#!/usr/bin/env python3
from pathlib import Path
import json,pandas as pd,cellxgene_census as c
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';(W/'cache').mkdir(exist_ok=True)
def main():
 with c.open_soma(census_version='2025-11-08') as f:
  ds=f['census_info']['datasets'].read().concat().to_pandas();a=ds[ds.collection_id.eq('ced320a1-29f3-47c1-a735-513c7084d508')];a.to_csv(O/'AIDA_Census_dataset_records.csv',index=False);print(a[['dataset_id','dataset_version_id','dataset_title','dataset_total_cell_count']].to_string(index=False),flush=True)
  a=a[a.dataset_title.str.contains('Data Freeze v2',case=False,regex=False)];assert len(a)==1,'Require exact Freezev2 record';a.to_csv(O/'AIDA_Census_selected_v2.csv',index=False)
  did=a.iloc[0].dataset_id;obs=f['census_data']['homo_sapiens'].obs.read(value_filter=f"dataset_id == '{did}'",column_names=['soma_joinid','donor_id','cell_type','raw_sum','is_primary_data','observation_joinid','dataset_id','disease','sex']).concat().to_pandas();obs.to_parquet(W/'cache/AIDA_census_obs.parquet',index=False);print('AIDAOBS',len(obs),obs.donor_id.nunique(),obs.cell_type.value_counts().to_dict(),flush=True)
  var=f['census_data']['homo_sapiens'].ms['RNA'].var.read().concat().to_pandas();var.to_csv(O/'Census_global_var.csv',index=False)
 print('PREFLIGHT_DONE',flush=True)
if __name__=='__main__':main()
