#!/usr/bin/env python3
"""Read only AIDA source donor/sample/library metadata, never the expression matrix."""
from pathlib import Path
import sys,json
import h5py,fsspec,aiohttp,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';sys.path.insert(0,str(P/'outputs/round3'))
from metadata_range_reader import MetadataRangeReader
from read_AIDA_source_metadata import col
def main():
 pre=json.loads((P/'outputs/round3/AIDA_source_metadata_preflight.json').read_text());local=P/'work/round3/AIDA_metadata_segments'/f"{pre['dataset_id']}.METADATA_ONLY.sparse.h5ad";receipt=pre['metadata_segment_receipt'];fs=fsspec.filesystem('https',client_kwargs={'trust_env':True,'timeout':aiohttp.ClientTimeout(total=180)})
 # Save additional receipt in round5; prior manifests are not rewritten.
 reader=MetadataRangeReader(local,pre['url'],fs,receipt,O/'P4_AIDA_additional_metadata_receipt.json')
 with reader,h5py.File(reader,'r') as h:
  keys=['donor_id','sample_id','library_id','suspension_uuid','sample_collection_method','disease','assay'];d=pd.DataFrame({k:col(h['obs'][k]) for k in keys})
 d=d.groupby(keys,dropna=False,observed=True).size().reset_index(name='n_cells');d.to_csv(O/'P4_AIDA_donor_sample_library.csv',index=False);z=d.groupby('donor_id').agg(n_samples=('sample_id','nunique'),n_libraries=('library_id','nunique'),n_suspensions=('suspension_uuid','nunique'));z.to_csv(O/'P4_AIDA_repeat_summary.csv');print(z.describe().to_string());print('MULTI_SAMPLE',(z.n_samples>1).sum(),'MULTI_LIBRARY',(z.n_libraries>1).sum(),'MULTI_SUSPENSION',(z.n_suspensions>1).sum(),flush=True)
if __name__=='__main__':main()
