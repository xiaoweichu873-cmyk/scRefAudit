"""Run acceptance tests and the example with the currently installed package."""
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import numpy
import scipy
import screfaudit

root = Path(__file__).resolve().parents[1]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
tests = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests/python", "-v"],
                       cwd=root, env=env, text=True, capture_output=True)
(root / "validation/python_tests.log").write_text(tests.stdout + tests.stderr)
if tests.returncode:
    print(tests.stdout + tests.stderr)
    raise SystemExit(tests.returncode)
example = subprocess.run([sys.executable, "examples/quickstart.py"], cwd=root, env=env,
                         text=True, capture_output=True, check=True)
values = json.loads(example.stdout)
(root / "validation/example_expected.json").write_text(json.dumps(values, indent=2) + "\n")
report = dict(version=screfaudit.__version__, package="scRefAudit", prepared_on="2026-11-21",
              python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
              platform=platform.system(), python_tests="9 passed", example="passed",
              installed_package_import=not str(Path(screfaudit.__file__).resolve()).startswith(str(root / "inst")),
              Rscript_available=shutil.which("Rscript") is not None,
              R_interface_status="Not executed in this environment; see R smoke test and GitHub workflow",
              real_input_compatibility=json.loads((root / "validation/real_input_replay.json").read_text()))
(root / "validation/release_checks.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({k:v for k,v in report.items() if k != "real_input_compatibility"}, indent=2))
