"""Verify the frozen v3 manifest and all table/figure source paths."""
from pathlib import Path
import csv,hashlib,json
root=Path(__file__).resolve().parents[1]
fail=[];count=0
for line in (root/'MANIFEST.sha256').read_text().splitlines():
    digest,rel=line.split('  ',1);p=root/rel;count+=1
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:fail.append(rel)
for name,col in [('SUPPLEMENTARY_TABLE_INDEX.csv','path'),('FIGURE_SOURCE_INDEX.csv','source_file')]:
    with (root/name).open() as f:
        for row in csv.DictReader(f):
            if not (root/row[col]).is_file():fail.append(row[col])
if fail:raise SystemExit('Missing or changed files: '+str(fail))
print(f'Verified {count} files and all indexed source paths.')
