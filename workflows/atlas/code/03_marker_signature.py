#!/usr/bin/env python
"""Run the markers phase; use the same --run-dir across phases."""
import sys
from stromal_atlas_pipeline import main

if __name__ == "__main__":
    main(["--stage", "markers", *sys.argv[1:]])
