#!/usr/bin/env python3
"""Metadata-only postprocessing of a completed PHARE prediction."""
from pathlib import Path
import json,pandas as pd,anndata as ad
from celltypist.models import Model
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1'
def main():
 f=O/'PHARE_real_K0_run.json';d=json.loads(f.read_text());manifest=json.loads((O/'PHARE_input_manifest.json').read_text());k0=pd.read_csv(O/'provenance_inputs/K0_primary_with_released_age.csv').set_index('gsm');a=ad.read_h5ad(P/manifest['input_path'],backed='r');symbols=set(a.var.feature_name.astype(str));a.file.close();model=Model.load(str(P/'work/round2/PHARE_runtime/phare/model/celltypist_model_from_panage_fairSample(10w).pkl'));required=set(map(str,model.features));matched=len(symbols&required);d.update(released_chronological_age=int(k0.loc[d['gsm'],'sample.subjectAgeAtDraw']),annotation_genes_matched=matched,annotation_gene_count=len(required),annotation_gene_coverage=matched/len(required));f.write_text(json.dumps(d,indent=2));print(d['released_chronological_age'],matched,len(required),d['prediction'][0]['predicted_Age'])
if __name__=='__main__':main()
