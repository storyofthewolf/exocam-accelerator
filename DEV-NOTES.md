# Development notes

Open bugs and loose ends, to deal with later.

## `accelerate.py` cannot import `trend_utils` outside the login environment

*Found 2026-11-12, hpc-connect live test D2. Open.*

Run through hpc-connect's `hpc-run` on Discover, `python accelerate.py --help`
dies before argparse runs:

    ModuleNotFoundError: No module named 'trend_utils'

`accelerate.py` does `import trend_utils as trend` (line 27, with a comment
pointing here), but `trend_utils` lives in the **exocam-trend** repo, not this
one. Run by hand, the import works only because the login environment or the
working directory supplies exocam-trend. An `hpc-run` payload runs in a
non-interactive shell that gets neither.

**Not the fix:** setting `PYTHONPATH` in the payload. hpc-connect then stops
treating the script as vetted (any loader variable set in the payload does
that), so it would need a token at a higher profile.

**Options:**

1. Locate exocam-trend relative to this file (both are checked out side by side
   under `~/pythonWorkspace/`) and add it to `sys.path` in `accelerate.py`
   before the import.
2. Install exocam-trend as a package in the GEOSpyD environment, so the import
   resolves from `$HPC_PYTHON_BIN/python` with no path handling.

Option 2 is cleaner, but it ties this repo to an install step on each machine;
option 1 needs no setup, but it assumes the side-by-side layout.
