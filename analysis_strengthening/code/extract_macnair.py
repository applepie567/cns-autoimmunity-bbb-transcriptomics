"""Stream published Matrix Market counts into donor/lesion/endothelial-subtype sums.
Only source-labelled endothelial nuclei in white matter are selected. Mixed
clusters and source exclude_pseudobulk nuclei are excluded independently of
the outcome genes. No cell or nucleus is treated as a biological replicate.
"""
from pathlib import Path
import gzip, json, time, os
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
P=Path(os.environ.get('BBI_RAW_DIR',str(ROOT/'raw')))/'Macnair'
meta=pd.read_csv(P/'ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz')
genes=pd.read_csv(P/'ms_lesions_snRNAseq_row_data_2023-09-12.txt.gz')
sel=meta.type_fine.str.startswith('Endocyte') & ~meta.exclude_pseudobulk & meta.matter.eq('WM')
sub=meta[sel].copy()
keycols=['individual_id_anon','lesion_type','type_fine']
keys=list(map(tuple,sub[keycols].values))
key_order=list(dict.fromkeys(keys)); key_idx={k:i for i,k in enumerate(key_order)}
mapcol=np.full(len(meta)+1,-1,dtype=np.int32)
mapcol[sub.index.to_numpy()+1]=[key_idx[k] for k in keys]
groups=pd.DataFrame(key_order,columns=keycols)
groups=groups.merge(sub.groupby(keycols).agg(nuclei=('cell_id','size'),sample_source=('sample_source','first'),age=('age_at_death','first'),sex=('sex','first'),pmi=('pmi_minutes','first')).reset_index(),on=keycols,validate='1:1')
counts=np.zeros((len(genes),len(groups)),dtype=np.int64)
cell_sums=np.zeros(len(meta)+1,dtype=np.int64)
expected=P/'ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz'
assert expected.stat().st_size==4563990935, 'Download not complete'
start=time.time();nlines=0;selected_entries=0
with gzip.open(expected,'rb') as f:
    header=f.readline();print(header.decode().strip(),flush=True)
    assert header.startswith(b'%%MatrixMarket matrix coordinate integer general')
    while True:
        line=f.readline()
        if not line.startswith(b'%'):break
    dims=list(map(int,line.split()));assert dims[:2]==[len(genes),len(meta)],dims
    print('DIMS',dims,'groups',len(groups),'selected nuclei',len(sub),flush=True)
    rest=b'';last=0
    while True:
        chunk=f.read(32*1024*1024)
        if not chunk:
            assert not rest.strip();break
        chunk=rest+chunk;cut=chunk.rfind(b'\n');rest=chunk[cut+1:]
        arr=np.fromstring(chunk[:cut].decode('ascii'),dtype=np.int64,sep=' ').reshape(-1,3)
        nlines+=len(arr)
        keep=mapcol[arr[:,1]]>=0
        a=arr[keep];selected_entries+=len(a)
        np.add.at(counts,(a[:,0]-1,mapcol[a[:,1]]),a[:,2])
        np.add.at(cell_sums,a[:,1],a[:,2])
        if nlines-last>=50000000:
            last=nlines;print('entries',nlines,'/',dims[2],'seconds',int(time.time()-start),flush=True)
    assert nlines==dims[2],(nlines,dims)
sub['reconstructed_UMI']=cell_sums[sub.index.to_numpy()+1]
sub['expected_UMI_from_metadata']=10**sub.log10_counts
ratio=sub.reconstructed_UMI/sub.expected_UMI_from_metadata
print('UMI_RATIO',ratio.describe().to_dict(),flush=True)
# Original QC counts can include transcripts subsequently removed by gene QC.
# Check that the file ordering gives close agreement with per-nucleus metadata.
assert ratio.between(.85,1.02).mean()>.99,ratio.describe()
groups.to_csv(ROOT/'derived/Macnair_pseudobulk_groups.csv',index=False)
sub.to_csv(ROOT/'derived/Macnair_selected_nuclei_QC.csv',index=False)
np.savez_compressed(ROOT/'derived/Macnair_WM_endothelial_counts.npz',counts=counts,genes=genes.symbol.fillna(genes.ensembl).str.upper().to_numpy(dtype=str))
(ROOT/'derived/Macnair_extraction_audit.json').write_text(json.dumps(dict(matrix_shape=dims[:2],entries=nlines,selected_entries=selected_entries,nuclei=len(sub),groups=len(groups),umi_ratio_min=float(ratio.min()),umi_ratio_median=float(ratio.median()),seconds=time.time()-start),indent=2))
print('DONE',flush=True)
