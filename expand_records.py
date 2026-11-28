"""Expand the archived research records and verify their SHA-256 checksums."""
import hashlib
from pathlib import Path
from zipfile import ZipFile

root = Path(__file__).resolve().parent
for volume in sorted((root / "records").glob("*.zip")):
    with ZipFile(volume) as archive:
        for item in archive.infolist():
            path = (root / item.filename).resolve()
            if not path.is_relative_to(root):
                raise ValueError(f"Invalid archive path {item.filename}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(archive.read(item))

manifest = (root / "checksums" / "expanded_records.sha256").read_text()
for line in manifest.splitlines():
    expected, name = line.split("  ", 1)
    if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
        raise ValueError(f"Checksum mismatch: {name}")
print(f"Verified {len(manifest.splitlines())} extracted records.")
