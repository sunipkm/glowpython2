import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parent.parent / "Examples"
EXAMPLES = sorted(EXAMPLES_DIR.glob("*.py"))

# Runs all of the python files in Examples/
# Just a smoke test, does not actually verify their outputs. 
# Just makes sure they can be run!
# New files detected auromatically
@pytest.mark.parametrize("script", EXAMPLES, ids=[e.name for e in EXAMPLES])
def test_example(script):
    result = subprocess.run(
        [sys.executable, str(script)],
        env=os.environ,
        capture_output=True,
        cwd=script.parent,
    )
    assert result.returncode == 0, result.stderr.decode()
