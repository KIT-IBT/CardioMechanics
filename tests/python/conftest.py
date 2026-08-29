import sys
from pathlib import Path

# Make the tools/python scripts importable as modules for unit testing.
TOOLS_PYTHON = Path(__file__).resolve().parents[2] / "tools" / "python"
if str(TOOLS_PYTHON) not in sys.path:
    sys.path.insert(0, str(TOOLS_PYTHON))
