#!/usr/bin/env python3
"""Phase entry point; all computation lives in signature_pipeline.py."""
import sys
from run_signature_score import main
if __name__ == "__main__":
    if "--stage" not in sys.argv: sys.argv.extend(["--stage", "audit"])
    main()
