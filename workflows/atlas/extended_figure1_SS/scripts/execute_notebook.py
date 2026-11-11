#!/usr/bin/env python3
"""Execute the delivered notebook cells with IPython, without a notebook server."""
import os
os.environ.setdefault("IPYTHONDIR", "/tmp/stromalatlas_signature_ipython")
from pathlib import Path
import time
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
root = Path(__file__).resolve().parents[1]
path = root / "notebooks/ExtendFig1_SignatureScore.ipynb"
nb = nbformat.read(path, as_version=4)
shell = InteractiveShell.instance()
count = 0
for cell in nb.cells:
    if cell.cell_type != "code":
        continue
    count += 1
    with capture_output() as captured:
        result = shell.run_cell(cell.source, store_history=False)
    outputs = []
    if captured.stdout:
        outputs.append(nbformat.v4.new_output("stream", name="stdout", text=captured.stdout))
    if captured.stderr:
        outputs.append(nbformat.v4.new_output("stream", name="stderr", text=captured.stderr))
    for output in captured.outputs:
        outputs.append(nbformat.v4.new_output("display_data", data=output.data, metadata=output.metadata))
    if result.error_before_exec or result.error_in_exec:
        raise RuntimeError(f"Notebook cell {count} failed") from (result.error_in_exec or result.error_before_exec)
    cell.outputs = outputs
    cell.execution_count = count
nb.metadata["execution"] = {"executor": "IPython InteractiveShell in-process; original notebook code cells", "executed_at": time.strftime("%Y-%m-%d %H:%M:%S"), "code_cells": count, "errors": 0}
nbformat.write(nb, path)
print(f"Executed {count} notebook code cells; no errors: {path}")
