from pathlib import Path
import gzip,csv,sys,os
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'baseline/BBI_v2.0.0/code'))
import run_bbi_extension as old
p=Path(os.environ.get('BBI_RAW_DIR',str(R/'raw')))/'GSE199460'
m=pd.read_csv(p/'GSE199460_cell_annotation.meta_data.cd31_selection.csv.gz',index_col=0)
m.index=m.index.str.replace('-1$','.1',regex=True)
values=[];genes=[]
with gzip.open(p/'GSE199460_normalized_expr.cd31_selection.sctransform.csv.gz','rt') as f:
    header=next(csv.reader([f.readline()]))
    cm=m.reindex(header)
    samples=sorted(cm['orig.ident'].dropna().unique())
    indices=[np.flatnonzero(cm['orig.ident'].eq(s)) for s in samples]
    for j,line in enumerate(f):
        gene,_,nums=line.partition(',');a=np.fromstring(nums,sep=',')
        assert len(a)==len(header)
        genes.append(gene.strip('"').upper());values.append([np.nanmean(a[ix]) for ix in indices])
        if (j+1)%5000==0:print('rows',j+1,flush=True)
groups=['Acute' if s.startswith('EAE') else 'Control' for s in samples]
payload=old.gene_payload('GSE199460',genes,np.array(values),samples,groups,'Acute','Control')
np.savez_compressed(R/'derived/GSE199460_payload.npz',**payload)
print('DONE',len(genes),samples,flush=True)
