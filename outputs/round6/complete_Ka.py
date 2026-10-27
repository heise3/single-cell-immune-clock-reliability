#!/usr/bin/env python3
"""Complete the anonymous access gate before any biological-pair analysis."""
from check_Ka_access import *
import pandas as pd
def probe(item):
 name,ref=item;url=ref['azul_url'];z={'file':name,'request_url':url,'anonymous':True,'credential_signature_or_DUA_supplied':False}
 try:
  with requests.get(url,stream=True,headers={'Range':'bytes=0-1023'},timeout=(20,45)) as r:
   b=next(r.iter_content(1024));z.update(http_status=r.status_code,bytes_inspected=len(b),gzip_magic=b[:2].hex(),final_host=urlsplit(r.url).netloc,content_range=r.headers.get('Content-Range'),download_access_verified=r.status_code in [200,206] and b[:2]==b'\x1f\x8b')
 except Exception as e:z.update(error=repr(e),download_access_verified=False)
 return z
def main():
 ref=json.loads((O/'HCA_declared_file_records.json').read_text());jobs=[(n,ref[n]) for n in ['cd4_t_cells.tar.gz','conventional_cd8_t_cells.tar.gz','myeloid_cells.tar.gz','nk_cells.tar.gz']];probes=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(probe,x) for x in jobs]):probes.append(f.result())
 (O/'K_a_remaining_type_access_probes.json').write_text(json.dumps(probes,indent=2));assert all(z['download_access_verified'] for z in probes),probes
 receipts=json.loads((O/'K_a_anonymous_access_receipts.json').read_text());r=next(x for x in receipts if x['source']=='b_cells.tar.gz');assert r['http_status']==200 and r['bytes']==ref['b_cells.tar.gz']['size'] and r['sha256']==ref['b_cells.tar.gz']['sha256'];fmt=pd.read_csv(O/'K_a_B_matrix_format_QA.csv').set_index('property').value;assert float(fmt.reconstructed_fraction_integer_tol1e6 if 'reconstructed_fraction_integer_tol1e6' in fmt.index else fmt['reconstructed_fraction_integer_tol1e-6'])==1.;assert float(fmt.reconstructed_max_rounding_error)<1e-8 and float(fmt.max_absolute_library_sum_error)<1e-8 and float(fmt.max_abs_feature_count_vs_metadata_error)==0
 status={'gate':'K-a','status':'PASS','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'access':'Barchive222791692bytes anonymouslyGET200 andSHA256 matches HCA declared object;all four other required type archives anonymousGET200/206 withgzipmagic;no login/institutional signature/applicationDUA','license':'HCA portal generic DataUseAgreement page specifiesCC-BY4.0;not an institution-signed/restrictedaccess application','matrix':'ReleasedBplain dgCMatrix is log1pCP10K;all86491453nonzero values losslessly recover integer counts usingexpm1(x)*nCount_RNA/10000;maxroundingerror1.819e-12;per-cell recoveredlibrarytotal matchesmetadata atsame tolerance;genes expressed exactlymatchesnFeature_RNA. This is numerically verified count reconstruction,not treating logvalues ascounts.','scope':'K-a source access plus A1-compatible count availability verified;other type files must pass same count-inversion QA before any clock scoring.','metadata':'Anonymous CSV containsDonor_id,Tube_id,Age,File_name,Batch;nCount_RNA,nFeature_RNA,percent.mt;donor-biological-visit counts/time semantics belong to K-b,not assumed here.'};(O/'K_a_status.json').write_text(json.dumps(status,indent=2));f=P/'outputs/gate_ledger.json';z=json.loads(f.read_text());z['round6']['K_a_result']=status;f.write_text(json.dumps(z,ensure_ascii=False,indent=2));print(json.dumps(status,indent=2),flush=True)
if __name__=='__main__':main()
