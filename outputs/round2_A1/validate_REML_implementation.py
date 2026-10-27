#!/usr/bin/env python3
"""Method QA only: synthetic no-known-noise special case vs native lme4 REML."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';sys.path.insert(0,str(O));from compute_reliability import variance_components
rng=np.random.default_rng(20261025);n=120;di=np.repeat(np.arange(24),5);bi=np.arange(n)%12;don=rng.normal(0,3,24);batch=rng.normal(0,1.5,12);y=don[di]+batch[bi]+rng.normal(0,2,n);g=pd.DataFrame({'donor':di.astype(str),'batch':bi.astype(str),'cell_sampling_var':0.,'algorithm_var':0.,'y':y});g.to_csv(O/'synthetic_REML_validation_input.csv',index=False);out=variance_components(g,y);(O/'synthetic_REML_python.json').write_text(json.dumps(out,indent=2));print(out)
