#!/usr/bin/env python3
"""Read the public DEG Explorer through its normal Dash UI callbacks.
Opaque serialized Store values are relayed unchanged, NEVER deserialized,
executed, or modified to name a different server file or query plan.
"""
from pathlib import Path
import json,requests,pandas as pd
P=Path(__file__).resolve().parents[2];W=P/'work/round2_A1/cache';O=P/'outputs/round2_A1'
URL='https://allenimmunology.org/public/proxy/dash/3512cd7c-1633-4b5b-a431-82c3f7bf8cf3/_dash-update-component'
def main():
 layout=json.loads((W/'RAM_dash-layout.json').read_text());deps=json.loads((W/'RAM_dash-dependencies.json').read_text());nodes={}
 def walk(z):
  if isinstance(z,dict):
   p=z.get('props',{})
   if p.get('id'):nodes[str(p['id'])]=p
   for v in z.values():walk(v)
  elif isinstance(z,list):
   for v in z:walk(v)
 walk(layout);session=requests.Session()
 def call(output,values,state=None):
  dep=next(d for d in deps if d['output']==output);multi=output.startswith('..')
  ops=[x.rsplit('.',1) for x in output.strip('.').split('...')] if multi else [output.rsplit('.',1)]
  out=[{'id':i,'property':p} for i,p in ops]
  inputs=[{**i,'value':values.get(i['id']+'.'+i['property'])} for i in dep['inputs']]
  states=[{**i,'value':(state or {}).get(i['id']+'.'+i['property'])} for i in dep['state']]
  body={'output':output,'outputs':out if multi else out[0],'inputs':inputs,'state':states,'changedPropIds':[i['id']+'.'+i['property'] for i in dep['inputs']]}
  r=session.post(URL,json=body,timeout=45);print('CALL',output[:35],r.status_code,len(r.content),flush=True);r.raise_for_status();return r.json()['response']
 celltypes=pd.read_excel(P/'work/sources/41586_2025_9686_MOESM5_ESM.xlsx',sheet_name='Fig1e_2');celltypes=celltypes[(celltypes.Expression=='Up-regulated')&(celltypes['count']>20)].celltype.str.strip().to_list()
 grid=next(d['output'] for d in deps if 'volcano_deg_grid.rowData' in d['output'])
 all_rows=[]
 for ct in celltypes:
  dest=W/('RAM_DEG_'+ct.replace('/','_').replace(' ','_')+'.json')
  if dest.exists():rows=json.loads(dest.read_text())
  else:
   b=call('barcode.data',{'covariate-labels.data':nodes['covariate-labels']['data'],'covariate-values.data':['Subject Age: Younger | Older',ct]}, {'covariates-df-store.data':nodes['covariates-df-store']['data']})['barcode']['data']
   raw=call('volcano-deg-df-1.data',{'barcode.data':b})['volcano-deg-df-1']['data']
   # The public UI's .2 effect cutoff colors/labels the full DEG table.
   # The paper's .1 RAM filter is applied locally to raw effect sizes.
   df=call('volcano-deg-df.data',{'volcano-deg-df-1.data':raw,'p-cutoff-store.data':.05,'es-cutoff-store.data':.2})['volcano-deg-df']['data']
   ans=call(grid,{'volcano-deg-df.data':df,'gene-selection-store.data':nodes['gene-selection-store']['data'],'toggle-dark-mode.value':True,'tertiary-tabs.active_tab':'tab-1','volcano_deg_grid.columnState':None,'volcano_deg_grid.filterModel':None})
   rows=ans['volcano_deg_grid']['rowData'];dest.write_text(json.dumps(rows));print('ROWS',ct,len(rows),'keys',list(rows[0]) if rows else [],flush=True)
  f=pd.DataFrame(rows);f['RAM_celltype']=ct;all_rows.append(f)
 pd.concat(all_rows,ignore_index=True).to_csv(O/'RAM_public_DEG_results.csv',index=False)
if __name__=='__main__':main()
