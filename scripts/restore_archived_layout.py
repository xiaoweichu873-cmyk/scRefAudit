"""Optionally restore the four byte-identical copies expected by the full archived range driver."""
from pathlib import Path
from hashlib import sha256
from shutil import copyfile
import json

root=Path(__file__).resolve().parents[1]
review=root/'review'
checks=json.loads((review/'range_extension/input_sha256.json').read_text())
for rel,digest in checks.items():
    dest=review/rel
    if dest.name=='inputs.npz':
        shared=root/'timing/data/prepared'/dest.parent.name/'inputs.npz'
        if sha256(shared.read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Input checksum mismatch: {shared}')
        dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():copyfile(shared,dest)
    if sha256(dest.read_bytes()).hexdigest()!=digest:
        raise ValueError(f'Archived-layout checksum mismatch: {rel}')
print('Four original task archives and their metadata match the archived range protocol.')
