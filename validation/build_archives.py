"""Build the GitHub ZIP and R source archive without claiming an R CMD check.

Build a Python wheel first with:
  python -m pip wheel . --no-deps --no-build-isolation --wheel-dir ../dist
Run this script after validation. Output is in the sibling dist directory.
"""
import hashlib
import json
from pathlib import Path
import re
import tarfile
import tempfile
import zipfile

root = Path(__file__).resolve().parents[1]
dist = root.parent / "dist"
dist.mkdir(exist_ok=True)
version = "0.1.0"


def source_files():
    files = []
    for p in root.rglob("*"):
        rel = p.relative_to(root)
        if not p.is_file() or p.is_symlink():
            continue
        if any(part in {"__pycache__", ".git", "build", "dist"} or part.endswith(".egg-info")
               for part in rel.parts):
            continue
        if p.suffix in {".pyc", ".pyo"}:
            continue
        files.append(p)
    return sorted(files)


files = source_files()
archive = dist / f"scRefAudit_{version}_GitHub.zip"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6,
                     allowZip64=False) as z:
    for p in files:
        z.write(p, "scRefAudit/" + p.relative_to(root).as_posix())

ignore = [re.compile(line) for line in (root / ".Rbuildignore").read_text().splitlines() if line]
r_files = []
for p in files:
    rel = p.relative_to(root)
    names = [rel.as_posix()] + [parent.as_posix() for parent in rel.parents if str(parent) != "."]
    if any(regex.search(name) for regex in ignore for name in names):
        continue
    if rel.parts[0] in {".Rbuildignore", ".gitignore"}:
        continue
    r_files.append(p)

r_archive = dist / f"scRefAudit_{version}.tar.gz"
with tarfile.open(r_archive, "w:gz", compresslevel=6) as tar:
    for p in r_files:
        info = tar.gettarinfo(str(p), "scRefAudit/" + p.relative_to(root).as_posix())
        info.uid = info.gid = 0
        info.uname = info.gname = ""
        info.mode = 0o644
        with p.open("rb") as f:
            tar.addfile(info, f)

# Actual extraction and byte comparison, including the R archive.
with tempfile.TemporaryDirectory(prefix="screfaudit_check_") as tmp:
    tmp = Path(tmp)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        z.extractall(tmp / "zip")
    with tarfile.open(r_archive) as tar:
        tar.extractall(tmp / "r", filter="data")
    for kind, selected in [("zip", files), ("r", r_files)]:
        for p in selected:
            dest = tmp / kind / "scRefAudit" / p.relative_to(root)
            assert hashlib.sha256(dest.read_bytes()).digest() == hashlib.sha256(p.read_bytes()).digest()

records = []
for p in [archive, r_archive, dist / f"screfaudit-{version}-py3-none-any.whl"]:
    records.append(dict(file=p.name, size_bytes=p.stat().st_size,
                        sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
result = dict(files=records, zip_file_count=len(files), R_source_file_count=len(r_files),
              zip_extraction="passed", R_tar_extraction="passed", R_CMD_check="not run; R unavailable")
(dist / "release_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
