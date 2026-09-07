from pathlib import Path
import sys,os
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'baseline/BBI_v2.0.0/code'))
import run_bbi_extension as old
old.DATA=Path(os.environ.get('BBI_RAW_DIR',str(ROOT/'raw')))
for acc,func in [('GSE95401',old.analyze_gse95401),('GSE199460',old.analyze_gse199460),('GSE279183',old.analyze_gse279183),('GSE210776',old.analyze_gse210776)]:
    dest=ROOT/'derived'/f'{acc}_payload.npz'
    if dest.exists():print(acc,'cached',flush=True);continue
    print(acc,'starting',flush=True)
    scores,coverage,qc,payload=func()
    np.savez_compressed(dest,**payload)
    qc.to_csv(ROOT/'derived'/f'{acc}_reconstructed_QC.csv',index=False)
    scores.to_csv(ROOT/'derived'/f'{acc}_reconstructed_scores.csv',index=False)
    print(acc,payload['expression'].shape,payload['groups'],flush=True)
