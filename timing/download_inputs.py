"""Fetch public processed inputs with checksums; does not access private data."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parent
base='https://github.com/rasmushenningsson/SingleCellExampleData/releases/download/GSE164378_RNA_ADT_3P_v2/'
files=[('https://ndownloader.figshare.com/files/34464122','kang_2018.h5ad'),(base+'GSE164378_RNA_ADT_3P.csv.gz','hao_annotations.csv.gz')]
files += [(base+f'GSE164378_RNA_ADT_3P_{d}.h5',f'hao_{d}.h5') for d in ['P1','P2','P3','P4']]
for url,name in files:
    dest=ROOT/'data'/name
    if dest.exists():print('Already present',dest);continue
    subprocess.run([sys.executable,str(ROOT/'download_public.py'),url,str(dest)],check=True)
