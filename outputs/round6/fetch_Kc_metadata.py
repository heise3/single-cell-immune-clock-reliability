#!/usr/bin/env python3
"""Read only required same-source metadata from tar prefixes, before K-c."""
from check_Ka_access import *
import tarfile,io,shutil,os
class AuditStream(io.RawIOBase):
 def __init__(self,stream):self.stream=stream;self.n=0;self.h=hashlib.sha256()
 def readable(self):return True
 def read(self,size=-1):
  b=self.stream.read(size);self.n+=len(b);self.h.update(b);return b
def one(item):
 ct,record=item;url=record['azul_url'];url=('https://humancellatlas.s3.amazonaws.com/'+record['azul_mirror_uri'].split('s3://humancellatlas/',1)[1]) if os.environ.get('R6_USE_DECLARED_MIRROR') else url;z={'celltype':ct,'url':url,'purpose':'metadataonly beforeKc;stoparchive stream afterCSV,not completearchive','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with requests.get(url,stream=True,timeout=(45,75)) as r:
   z['http_status']=r.status_code;r.raise_for_status();r.raw.decode_content=False;s=AuditStream(r.raw);tf=tarfile.open(fileobj=s,mode='r|gz');members=[]
   for m in tf:
    members.append({'name':m.name,'bytes':m.size})
    if m.isfile() and m.name.endswith('_metadata.csv'):
     out=W/'extracted'/f'{ct}_metadata.csv';out.parent.mkdir(exist_ok=True)
     with out.open('wb') as f:shutil.copyfileobj(tf.extractfile(m),f)
     z.update(metadata_path=str(out.relative_to(P)),metadata_bytes=out.stat().st_size,metadata_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),archive_bytes_streamed=s.n,archive_prefix_sha256=s.h.hexdigest(),members_before_metadata=members);break
   else:raise RuntimeError('No metadataCSV beforearchiveend')
 except Exception as e:z['error']=repr(e)
 return z
def main():
 assert json.loads((O/'K_b_status.json').read_text())['status']=='PASS';refs=json.loads((O/'HCA_declared_file_records.json').read_text());jobs=[(ct,refs[name]) for ct,name in [('CD4T','cd4_t_cells.tar.gz'),('CD8T','conventional_cd8_t_cells.tar.gz'),('MONO','myeloid_cells.tar.gz'),('NK','nk_cells.tar.gz')]];jobs=[j for j in jobs if not os.environ.get('R6_CT') or j[0]==os.environ['R6_CT']];rows=json.loads((O/'K_c_metadata_download_receipts.json').read_text()) if (O/'K_c_metadata_download_receipts.json').exists() else []
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(one,x) for x in jobs]):
   z=f.result();rows.append(z);(O/'K_c_metadata_download_receipts.json').write_text(json.dumps(rows,indent=2));print('METADATA',z['celltype'],z.get('metadata_bytes'),z.get('archive_bytes_streamed'),z.get('error',''),flush=True)
if __name__=='__main__':main()
