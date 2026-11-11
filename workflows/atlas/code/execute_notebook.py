#!/usr/bin/env python
"""Execute the notebook against a completed or new run using the same phase driver."""
import argparse
import os
from pathlib import Path
import nbformat
from nbclient import NotebookClient

import stromal_atlas_pipeline as pipeline


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,default=pipeline.HERE.parent/'results/runs/atlas_metacell_20261031')
    p.add_argument('--config',type=Path,default=pipeline.HERE/'config.json')
    p.add_argument('--output',type=Path,default=pipeline.HERE/'stromal_atlas_pipeline.ipynb')
    args=p.parse_args()
    for key in ['IPYTHONDIR','JUPYTER_CONFIG_DIR','JUPYTER_DATA_DIR','JUPYTER_RUNTIME_DIR']:
        folder=pipeline.HERE.parent/'tmp'/key.lower();folder.mkdir(exist_ok=True)
        os.environ.setdefault(key,str(folder))
    nb=nbformat.read(pipeline.HERE/'stromal_atlas_pipeline.ipynb',as_version=4)
    setup=next(c for c in nb.cells if c.cell_type=='code' and 'RUN_DIR =' in c.source)
    setup.source='\n'.join('RUN_DIR = Path('+repr(str(args.run_dir.resolve()))+')' if line.startswith('RUN_DIR =')
                           else 'CONFIG = Path('+repr(str(args.config.resolve()))+')' if line.startswith('CONFIG =')
                           else line for line in setup.source.splitlines())
    client=NotebookClient(nb,timeout=7200,kernel_name='python3',resources={'metadata':{'path':str(pipeline.HERE)}})
    client.execute()
    nbformat.validate(nb)
    nbformat.write(nb,args.output)
    print('Executed notebook:',args.output)
    print('Code cells:',sum(c.cell_type=='code' for c in nb.cells))


if __name__=='__main__':main()
