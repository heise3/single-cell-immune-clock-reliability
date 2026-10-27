#!/usr/bin/env python3
"""Summarize observed anonymous access; never infer DUA from HTTP403."""
from pathlib import Path
import json, datetime, hashlib
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1'
def main():
 api=json.loads((O/'siAge_alternate_API_queries.json').read_text())
 syn=json.loads((O/'siAge_fresh_Synapse_access.json').read_text())
 model=P/'work/siAgeModel.rds'
 d={'seed':20261025,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'status':'SKIPPED_NATIVE_ANNOTATION_REFERENCE_NOT_ANONYMOUSLY_ACQUIRED',
 'native_weights_available':True,'native_weight_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),
 'native_repository':'https://github.com/shanzha9/siAgeModel','pinned_commit':'3e8e17567afa296213b414083eaca59b30455b81',
 'required_reference':'syn61609846 scRNA-seqProcessedLabelledObject.rds, required by native Readme.MD and predictCellType preprocessing',
 'routes':[{'channel':'GitHub releases','url':'https://api.github.com/repos/shanzha9/siAgeModel/releases','observed_response':[],'n_releases':0},
 {'channel':'GitHub model directory','url':'https://github.com/shanzha9/siAgeModel/tree/3e8e17567afa296213b414083eaca59b30455b81/model','files':['all_selectedModel.rds'],'annotation_reference_present':False}]+api['routes']+syn['routes'],
 'third_party_caveat':'IDEAL-Age checkpoints labeled siAge are IDEAL-Age models tuned on the siAge dataset, not the original siAge randomForest model or native cell-label reference; no substitution.',
 'access_interpretation':'Public metadata accessible; anonymous file request403 with only READ permission. This establishes anonymous download failure, not institutional-signature or DUA requirements. No login, terms acceptance or access request performed.',
 'search_scope':'Zero results refer only to the recorded Zenodo/figshare API queries, not universal nonexistence.',
 'authorization':'A1 explicitly allows recording and skipping siAge if alternate public routes remain unavailable.',
 'downstream':'Native preprocessing and main-cohort siAge reliability NOT_RUN; no biological negative inferred.'}
 assert syn['routes'][1]['http_status']==403
 (O/'siAge_alternate_routes.json').write_text(json.dumps(d,indent=2,ensure_ascii=False))
 print(d['status'])
if __name__=='__main__':main()
