import csv
import collections
import json
import pathlib
import re

root = pathlib.Path(__file__).resolve().parent
out = root.parent / 'outputs' / 'HIHA_gate_audit'
out.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader((root / 'geo_sample_metadata.csv').open()))
donors = collections.defaultdict(list)
for row in rows:
    donors[row['subject id']].append(row)

eligible = []
donor_rows = []
for donor, samples in sorted(donors.items()):
    label_counts = collections.Counter(s['visit'] for s in samples)
    ambiguous_visit_labels = any(n > 1 for n in label_counts.values())
    covid_fields = [k for k in samples[0] if k.startswith('covid vax_')]
    covid = {int(s[k]) for s in samples for k in covid_fields if s[k] != 'NA'}
    complete_covid = all(any(s[k] != 'NA' for s in samples) for k in covid_fields)
    flu = {int(s['days since_first_visit']) for s in samples
           if s['visit'].startswith('Flu') and s['visit'].endswith('Day 0')}
    counts = collections.Counter()
    for sample in samples:
        day = int(sample['days since_first_visit'])
        flu_followup = bool(sample['visit'].startswith('Flu') and re.search(r'Day (7|90)$', sample['visit']))
        in_flu_window = any(0 < day - v <= 90 for v in flu)
        in_covid_window = any(0 <= day - v <= 90 for v in covid)
        nominal = not flu_followup
        known_events = nominal and not in_flu_window and not in_covid_window
        conservative = known_events and complete_covid and not ambiguous_visit_labels
        counts['nominal'] += int(nominal)
        counts['known_events'] += int(known_events)
        counts['conservative'] += int(conservative)
        eligible.append({**sample, 'nominal_steady': nominal,
                         'known_event_steady': known_events,
                         'complete_covid_dates': complete_covid,
                         'ambiguous_visit_labels': ambiguous_visit_labels,
                         'conservative_steady': conservative})
    donor_rows.append({'donor': donor, 'cmv': samples[0]['cmv status'],
                       'n_visits': len({s['days since_first_visit'] for s in samples}),
                       'n_unique_visit_labels': len(label_counts),
                       'ambiguous_visit_labels': ambiguous_visit_labels,
                       'complete_covid_dates': complete_covid, **dict(counts)})

def write_csv(name, records):
    with (out / name).open('w') as handle:
        writer = csv.DictWriter(handle, list(records[0]))
        writer.writeheader()
        writer.writerows(records)

write_csv('K0_sample_eligibility.csv', eligible)
write_csv('K0_donor_visits.csv', donor_rows)
summary = {
    'n_donors': len(donors), 'n_samples': len(rows),
    'n_batches': len({s['batch id'] for s in rows}),
    'batches': sorted({s['batch id'] for s in rows}),
    'cmv_donors': dict(collections.Counter(s[0]['cmv status'] for s in donors.values())),
    'visit_distribution': dict(sorted(collections.Counter(r['n_visits'] for r in donor_rows).items())),
    'n_ge4': {key: sum(r[key] >= 4 for r in donor_rows)
              for key in ['nominal', 'known_events', 'conservative']},
    'conservative_visit_distribution': dict(sorted(collections.Counter(r['conservative'] for r in donor_rows).items())),
    'conservative_n_samples_all_donors': sum(r['conservative'] for r in donor_rows),
    'conservative_n_samples_ge4_donors': sum(r['conservative'] for r in donor_rows if r['conservative'] >= 4),
    'window_definition': 'Exclude Flu Day 7/90 labels, all recorded influenza post-dose days (0,90], and recorded COVID dose days [0,90]. Influenza Day 0 blood draws precede vaccination. Conservative eligibility also requires non-missing first and second COVID dose dates for the donor.',
    'window_status': 'Analyst operational assumption, user did not specify a numeric post-vaccine window.',
    'metadata_ambiguity': 'BR2021 has two samples labeled Immune Variation Day 90 at study-relative days 0 and 364 (GSM8388368 and GSM8388374). Conservative K0 lower bound excludes this donor pending resolution.',
    'geo_bridging_control_flag': dict(collections.Counter(s['sample bridging_control'] for s in rows)),
    'bridge_interpretation': 'All 868 analyzed donor samples are marked FALSE. Separate batch-control aliquots from a common leukapheresis draw are documented and downloadable by Allen; these are not included among the 868 donor samples.',
}
(out / 'K0_metadata_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
