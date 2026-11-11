"""Synthetic atlas audit, read-only input and relocated-config checks."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]


class AtlasReleaseTests(unittest.TestCase):
    def test_audit_from_external_working_directory_and_resume(self):
        with tempfile.TemporaryDirectory(prefix='tabulastromal_atlas_release_') as tmp:
            tmp = Path(tmp)
            cfg = json.loads((ROOT / 'workflows/atlas/code/config.json').read_text())
            cfg['inputs'] = {}
            cfg['metacell_inputs'] = {}
            hashes = {}
            rng = np.random.default_rng(12)
            for lin, states in [('fibro', ['qFibro', 'eFibro']), ('endo', ['ArtEC', 'CapEC'])]:
                state = [states[i % 2] for i in range(24)]
                obs = pd.DataFrame({
                    'Celltype_general': state,
                    'Celltype_sub': [s + '_marker' for s in state],
                    'Celltype': [s + '_marker_pan' for s in state],
                    'sample': [f's{i // 4}' for i in range(24)],
                    'GSE': ['GSE_SYNTHETIC'] * 24,
                    'tissue': ['synthetic_tissue'] * 24,
                    'Age': [20 + i // 4 for i in range(24)],
                }, index=[f'{lin}_fakecell{i}' for i in range(24)])
                a = ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(2, (24, 3))).astype('float32')),
                               obs=obs, var=pd.DataFrame(index=['G1', 'G2', 'G3']))
                a.obsm['X_umap'] = rng.normal(size=(24, 2)).astype('float32')
                a.uns['log1p'] = {}
                path = tmp / f'{lin}.h5ad'
                a.write_h5ad(path)
                cfg['inputs'][lin] = path.name  # resolves against config, not cwd
                hashes[path] = hashlib.sha256(path.read_bytes()).hexdigest()
            config = tmp / 'atlas.json'
            config.write_text(json.dumps(cfg))
            output = tmp / 'run'
            command = [sys.executable, str(ROOT / 'scripts/run_workflow.py'), 'atlas',
                       '--config', str(config), '--run-dir', str(output), '--stage', 'audit']
            first = subprocess.run(command, cwd=tmp, text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stderr[-2000:])
            markers = {}
            for lin in ['fibro', 'endo']:
                audit = json.loads((output / f'shared/data_audit_{lin}.json').read_text())
                self.assertEqual(audit['n_cells'], 24)
                self.assertEqual(audit['n_samples'], 6)
                self.assertIsNone(audit['n_donors'])
                path = output / f'shared/data_audit_{lin}.json'
                markers[path] = path.stat().st_mtime_ns
            again = subprocess.run(command + ['--resume'], cwd=tmp, text=True, capture_output=True)
            self.assertEqual(again.returncode, 0, again.stderr[-2000:])
            for path, stamp in markers.items():
                self.assertEqual(path.stat().st_mtime_ns, stamp)
            for path, expected in hashes.items():
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected)


if __name__ == '__main__':
    unittest.main()
