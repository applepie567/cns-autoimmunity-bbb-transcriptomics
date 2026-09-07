"""Download the source inputs listed in output/INPUT_FILES.csv, checking MD5.

Usage: python code/download_inputs.py --raw-dir /absolute/path/to/source_data
Optional: --dataset Macnair or --dataset GSE199460
The derived matrices in this package suffice for downstream reproduction.
"""
from pathlib import Path
import argparse, csv, hashlib, urllib.request

def digest(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'md5').hexdigest()

def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,required=True)
    parser.add_argument('--dataset')
    args=parser.parse_args()
    with (root/'output/INPUT_FILES.csv').open() as f: rows=list(csv.DictReader(f))
    for row in rows:
        if args.dataset and row['dataset']!=args.dataset: continue
        dest=args.raw_dir/row['relative_download_path']; dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists() and dest.stat().st_size==int(row['bytes']) and digest(dest)==row['md5']:
            print('Verified existing',row['filename'],flush=True);continue
        partial=dest.with_name(dest.name+'.download')
        print('Downloading',row['url'],flush=True)
        with urllib.request.urlopen(row['url'],timeout=120) as src,partial.open('wb') as out:
            while block:=src.read(8*1024*1024): out.write(block)
        if partial.stat().st_size!=int(row['bytes']) or digest(partial)!=row['md5']:
            raise ValueError(f'Integrity check failed: {partial}')
        partial.replace(dest)
        print('Verified',dest,flush=True)

if __name__=='__main__': main()
