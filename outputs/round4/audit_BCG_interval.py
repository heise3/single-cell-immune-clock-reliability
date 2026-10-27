#!/usr/bin/env python3
"""Verify archived sc-ImmuAging Methods and bind the curated primary review."""
from pathlib import Path
import hashlib,json,html,re
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4'
def main():
 source=P/'work/sources/scImmuAging.html';s=source.read_text();text=html.unescape(re.sub('<[^>]+>','',s));text=' '.join(text.split());assert 'before (T0) and three months after (T3m) BCG vaccination' in text;assert 'EGAS00001006990' in text;assert 'Differentially expressed genes identification' in text
 evidence=pd.read_csv(O/'B_sampling_interval_primary_sources.csv');assert len(evidence)==2 and evidence.nominal_interval_months.eq(3).all();r=pd.read_csv(P/'outputs/round3/R1_critical_SEM_vs_CD8_estimates.csv');z=r[r.metric.eq('age_acceleration')&r.stratum.eq('75-105d')].iloc[0];out=pd.read_csv(O/'B_SEM_interval_match.csv');assert abs(out.SEM_CD8_age_acceleration.iloc[0]-z.SEM)<1e-12;assert z.SEM_CI_low>z.critical_SEM_Holm;(O/'B_interval_verification_QA.json').write_text(json.dumps({'status':'PASS','archived_scImmuAging_Methods_phrase_verified':True,'scImmuAging_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'original300BCG_Methods_review':'Curated primary-source review recorded in B_sampling_interval_primary_sources.csv; script does not pretend to automate full literature reading','matched_interval':'75-105d','exact_individual_dates':'not publicly supplied; nominal3months only','R1_values_copied_without_resimulation':True},indent=2));print('BCG_INTERVAL_QA_PASS')
if __name__=='__main__':main()
