"""Metadata-only frozen-window audit; no expression/model outcomes read."""
import collections
import csv
import hashlib
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / 'outputs/HIHA_gate_audit/K0_v2'
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / 'work/geo_sample_metadata.csv'
rows = list(csv.DictReader(source.open()))
donors = collections.defaultdict(list)
for row in rows:
    donors[row['subject id']].append(row)
eligibility, donor_records, summaries = [], [], []
for window in [30, 90, 180]:
    by_donor = []
    for donor, samples in sorted(donors.items()):
        fields = [k for k in samples[0] if k.startswith('covid vax_')]
        covid = {int(s[k]) for s in samples for k in fields if s[k] != 'NA'}
        complete = all(any(s[k] != 'NA' for s in samples) for k in fields)
        flu = {int(s['days since_first_visit']) for s in samples
               if s['visit'].startswith('Flu') and s['visit'].endswith('Day 0')}
        labels = collections.Counter(s['visit'] for s in samples)
        ambiguous = any(v > 1 for v in labels.values())
        days = set()
        for s in samples:
            day = int(s['days since_first_visit'])
            followup = bool(s['visit'].startswith('Flu') and re.search(r'Day (7|90)$', s['visit']))
            flu_window = any(0 < day - v <= window for v in flu)
            covid_window = any(0 <= day - v <= window for v in covid)
            eligible = complete and not ambiguous and not followup and not flu_window and not covid_window
            if eligible:
                days.add(day)
            eligibility.append({'window_days': window, 'donor': donor, 'gsm': s['gsm'],
                                'visit': s['visit'], 'relative_day': day,
                                'complete_covid_dates': complete,
                                'ambiguous_visit_labels': ambiguous,
                                'vaccine_followup': followup, 'in_flu_window': flu_window,
                                'in_covid_window': covid_window, 'eligible_visit': eligible})
        by_donor.append({'window_days': window, 'donor': donor, 'n_distinct_eligible_days': len(days),
                         'donor_ge4': len(days) >= 4, 'complete_covid_dates': complete,
                         'ambiguous_visit_labels': ambiguous})
    donor_records.extend(by_donor)
    eligible_donors = {r['donor'] for r in by_donor if r['donor_ge4']}
    visit_count = sum(r['n_distinct_eligible_days'] for r in by_donor if r['donor_ge4'])
    summaries.append({'window_days': window, 'role': 'PRIMARY' if window == 90 else 'SENSITIVITY',
                      'n_donors_total': len(donors), 'n_eligible_donors_ge4': len(eligible_donors),
                      'n_distinct_visits_in_ge4_donors': visit_count,
                      'n_eligible_samples_before_ge4': sum(r['eligible_visit'] for r in eligibility if r['window_days'] == window),
                      'eligible_visit_distribution_in_ge4': dict(sorted(collections.Counter(
                          r['n_distinct_eligible_days'] for r in by_donor if r['donor_ge4']).items())),
                      'original_K0_ge60': len(eligible_donors) >= 60})
for name, records in [('sample_eligibility.csv', eligibility), ('donor_eligibility.csv', donor_records)]:
    with (OUT / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
result = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
          'metadata_only': True, 'windows': summaries,
          'n_donors_missing_covid_dose_dates': sum(not r['complete_covid_dates'] for r in donor_records if r['window_days'] == 90),
          'ambiguous_donors': [r['donor'] for r in donor_records if r['window_days'] == 90 and r['ambiguous_visit_labels']],
          'interpretation': 'The >=60 threshold applies to the primary 90-day K0. Sensitivity failure does not replace or reverse the frozen primary gate.'}
(OUT / 'window_summary.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
