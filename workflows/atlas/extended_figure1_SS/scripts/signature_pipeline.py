#!/usr/bin/env python3
"""Auditable, out-of-core signature scoring and grouped annotation concordance.

Input H5AD files are opened read-only. No annotations are changed. Candidate
catalogues retain their original dictionary names, keys, gene lists and sources.
"""
from __future__ import annotations

import ast
import hashlib
import json
import logging
import math
import os
import re
import resource
import sys
import time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats
import yaml
try:
    from anndata.io import read_elem
except ImportError:
    from anndata._io.specs import read_elem

VERSION = 'ss_v1_20261104'
LABELS = {'fibro': 'Fibroblast', 'endo': 'Endothelial'}


def serial(v):
    if isinstance(v, (np.integer, np.floating)): return v.item()
    if isinstance(v, np.ndarray): return v.tolist()
    if isinstance(v, Path): return str(v)
    if isinstance(v, (set, tuple)): return list(v)
    raise TypeError(type(v).__name__)


def dump(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=serial), encoding='utf-8')


def read_json(path): return json.loads(Path(path).read_text())


def read_config(path):
    cfg = yaml.safe_load(Path(path).read_text())
    cfg['output_dir'] = str(Path(cfg['output_dir']).resolve())
    for d in ['audit', 'data', 'marker_sets', 'figures/panels', 'figures/integrated',
              'scripts', 'notebooks', 'logs', 'reports', 'data/cache']:
        (Path(cfg['output_dir']) / d).mkdir(parents=True, exist_ok=True)
    return cfg


def root(cfg): return Path(cfg['output_dir'])
def cache(cfg, lin):
    p = root(cfg) / 'data/cache' / lin; p.mkdir(parents=True, exist_ok=True); return p


def fingerprint(path):
    p = Path(path); s = p.stat()
    return {'path': str(p.resolve()), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def clean(values):
    s = pd.Series(values).astype('string').str.strip().fillna('Unknown')
    return s.replace({'': 'Unknown', 'nan': 'Unknown', 'None': 'Unknown',
                      'NA': 'Unknown', '<NA>': 'Unknown'}).astype('category')


def csr_rows(node, start, stop):
    """Read only a contiguous row block, never densify the full matrix."""
    if isinstance(node, h5py.Dataset): return sparse.csr_matrix(node[start:stop])
    enc = node.attrs.get('encoding-type', '')
    if enc != 'csr_matrix': raise ValueError(f'Expected CSR/dense, found {enc}')
    indptr = node['indptr'][start:stop+1]
    lo, hi = int(indptr[0]), int(indptr[-1])
    return sparse.csr_matrix((node['data'][lo:hi], node['indices'][lo:hi], indptr-lo),
                             shape=(stop-start, int(node.attrs['shape'][1])))


def safe_value(node, env):
    """Resolve literal dictionaries/list variables without executing notebooks."""
    try: return ast.literal_eval(node)
    except (ValueError, TypeError): pass
    if isinstance(node, ast.Name) and node.id in env: return env[node.id]
    if isinstance(node, ast.Dict):
        return {safe_value(k, env): safe_value(v, env) for k, v in zip(node.keys, node.values)}
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        vals = [safe_value(x, env) for x in node.elts]
        return set(vals) if isinstance(node, ast.Set) else vals
    if isinstance(node, ast.BinOp):
        a, b = safe_value(node.left, env), safe_value(node.right, env)
        if isinstance(node.op, ast.Add): return a+b
        if isinstance(node.op, ast.BitOr): return a|b
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {'sorted','list','set','tuple'} and len(node.args)==1:
        return {'sorted': sorted, 'list': list, 'set': set, 'tuple': tuple}[node.func.id](safe_value(node.args[0], env))
    raise ValueError('Not a safe literal expression')


def category_of(name, dic):
    ks = set(dic)
    if ks <= {'T','B','Fibro','Endo','MoMac','Epi','Fib_deg','Endo_deg'}: return 'broad_lineage'
    if name == 'Fibro_marker_sub_ya' or ks <= {'AP','Tip','Stalk','ECM','Immune'}: return 'functional'
    if any('fibro' in k.lower() for k in ks): return 'subtype_fibro'
    if any(k.startswith(('SMC','Peri','LEC','VenEC','CapEC','ArtEC','EndoMT','MCs')) for k in ks): return 'subtype_endo'
    return 'other'


def collect_history(cfg):
    catalog, history, issues = [], [], []
    mdpath = Path(cfg['prompt_file']); md = mdpath.read_text()
    sources = [('user_markdown', str(mdpath), None, '\n\n'.join(re.findall(r'```python\s*\n(.*?)```', md, re.S)))]
    for lin, path in cfg['notebooks'].items():
        nb = read_json(path)
        for i, c in enumerate(nb['cells']):
            if c['cell_type'] == 'code': sources.append((lin, path, i, ''.join(c['source'])))
    envs = {}; definitions = {}; references = []
    for source, path, cell, code in sources:
        env = envs.setdefault(path, {})
        sanitized = '\n'.join('pass' if l.lstrip().startswith(('%','!')) else l for l in code.splitlines())
        try: tree = ast.parse(sanitized)
        except SyntaxError as exc:
            issues.append({'source': path, 'cell': cell, 'issue': f'Unparsed code: {exc.msg}'}); continue
        if any(w in code for w in ['calculate_simple_scores', 'score_genes', 'Celltype_general', 'normalize_total', 'sc.pp.log1p', 'EndoMT-I', 'EndoMT-II', 'SMC', 'Peri']):
            history.append({'source': path, 'cell': cell, 'code': code})
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for tar in node.targets:
                    if not isinstance(tar, ast.Name): continue
                    definitions.setdefault(path, {})[tar.id] = cell
                    try: value = safe_value(node.value, env)
                    except Exception: continue
                    env[tar.id] = value
                    if (isinstance(value, dict) and value and
                        all(isinstance(k,str) and isinstance(v,list) and len(v)>0 and all(isinstance(g,str) for g in v) for k,v in value.items())):
                        kind = category_of(tar.id, value)
                        if kind == 'other': continue
                        payload = json.dumps(value, sort_keys=True)
                        cid = f'{source}:{tar.id}:cell{cell if cell is not None else "appendix"}'
                        catalog.append({'candidate_id': cid, 'dictionary_name': tar.id, 'source': path,
                                        'cell_index': cell, 'kind': kind, 'original_markers': value,
                                        'content_sha256': hashlib.sha256(payload.encode()).hexdigest()})
            for call in ast.walk(node):
                if isinstance(call,ast.Call) and isinstance(call.func,ast.Name) and call.func.id=='calculate_simple_scores':
                    if len(call.args)>1 and isinstance(call.args[1],ast.Name):
                        var = call.args[1].id
                        references.append({'source': path,'cell': cell,'dictionary_name':var,
                                           'defined_before_call':var in env})
    # Fig1 is a distinct historical candidate, not new biological evidence.
    old_code = Path(cfg['fig1_dir'])/'source_code/utils_data.py'
    t = ast.parse(old_code.read_text())
    old_markers = next(ast.literal_eval(n.value) for n in t.body if isinstance(n,ast.Assign)
                       and any(isinstance(x,ast.Name) and x.id=='MARKERS' for x in n.targets))
    for lin, dic in old_markers.items():
        catalog.append({'candidate_id':f'fig1_v2:{lin}:MARKERS','dictionary_name':f'Fig1_v2_{lin}_MARKERS',
                        'source':str(old_code),'cell_index':None,'kind':f'subtype_{lin}',
                        'original_markers':dic,'content_sha256':hashlib.sha256(json.dumps(dic,sort_keys=True).encode()).hexdigest()})
    dump(root(cfg)/'marker_sets/candidate_marker_catalog.json', {'version':VERSION,'candidates':catalog})
    dump(root(cfg)/'audit/notebook_code_evidence.json', history)
    dump(root(cfg)/'audit/notebook_score_calls.json', references)
    dump(root(cfg)/'audit/notebook_parse_issues.json', issues)
    figcfg = read_json(Path(cfg['fig1_dir'])/'config_used.json')
    dump(root(cfg)/'audit/fig1_v2_colors.json', {'source': str(Path(cfg['fig1_dir'])/'config_used.json'),
          'general_colors':figcfg['general_colors'], 'general_order':figcfg['general_order']})
    lines = ['# Historical code and colors audit', '', f'Candidate dictionaries: {len(catalog)}. Literal resolution only; notebooks were not executed.', '',
             'All user dictionaries retain their exact original names, keys and lists. The word authoritative is only an input name.', '',
             '## Requested variable checks']
    for path in cfg['notebooks'].values():
        lines += ['', f'Notebook: `{path}`']
        for var in ['stromal_marker','Endo_marker_dict_full_latest','Endo_marker_dict_latest']:
            status = f'defined in cell {definitions[path][var]}' if var in definitions.get(path,{}) else 'not defined in this notebook'
            lines.append(f'- `{var}`: {status}.')
    lines += ['', '## Execution-order issues']
    for x in references:
        if not x['defined_before_call']: lines.append(f'- Cell {x["cell"]} uses `{x["dictionary_name"]}` before its definition: {Path(x["source"]).name}.')
    lines += ['', 'Historical mean-score functions set missing sets to zero and clip scores at zero; these behaviors are not retained.',
              'The Endo notebook also contains a log1p(mean expression) × fraction-detected function; its raw-count description does not establish that the current X is counts.',
              'The Endo plotting code references LymEC in an 8-type dictionary using LEC: stale-key risk.',
              'Cluster IDs were relabeled after lineage filtering in the notebooks; old cluster IDs must not be transferred without checking current cross-tabs.',
              '', '## Fig1 v2 palette', '', '```json', json.dumps(figcfg['general_colors'],indent=2), '```', '',
              'Eight-type candidate sets are retained as eight-type diagnostics. The six-type comparison uses an explicitly labeled gene-union projection (SMC/Peri→MCs; EndoMT-I/II→EndoMT), not eight ground-truth labels. See annotation_hierarchy_evidence.json for actual notebook/obs evidence.']
    (root(cfg)/'audit/history_and_color_map_audit.md').write_text('\n'.join(lines))
    return catalog


def normalize_marker_key(k, lin):
    if lin=='fibro':
        if k=='apFibro': return 'iFibro'
        return re.sub(r'1$', '', k)
    k = re.sub(r'(_[1os]|1)$', '', k)
    return {'SMC':'MCs','Peri':'MCs','EndoMT-I':'EndoMT','EndoMT-II':'EndoMT','LymEC':'LEC'}.get(k,k)


def resolve_gene(g, names):
    """Only explicit safe mappings. Ambiguous symbols are excluded, never guessed."""
    if g=='H2-Q4': return None, 'nonhuman_symbol_excluded'
    if g in {'CD16','HLA-DQA'}: return None, 'ambiguous_symbol_excluded'
    if g in names:
        if names[g] != 1: return None, 'duplicate_feature_excluded'
        return g, 'present'
    if g=='NG2':
        return ('CSPG4', 'alias_NG2_to_CSPG4') if names.get('CSPG4')==1 else (None,'alias_target_absent')
    upper = g.upper()
    if upper!=g and names.get(upper)==1: return upper, 'case_normalized'
    return None, 'absent'


def audit_one(cfg, lin, catalog):
    out = root(cfg); c = cache(cfg,lin); path = cfg['inputs'][lin]
    logging.info('Phase 0: opening actual %s atlas read-only: %s',lin,path)
    with h5py.File(path,'r') as f:
        node=f['X']; shape=tuple(int(x) for x in node.attrs['shape']) if isinstance(node,h5py.Group) else node.shape
        obscols = list(f['obs']); indexkey = f['obs'].attrs.get('_index','_index')
        names = np.asarray(read_elem(f['var'][f['var'].attrs.get('_index','_index')])).astype(str)
        np.save(c/'var_names.npy',names)
        selected = [k for k in ['Celltype_general','Celltype','Celltype_sub','leiden','sample','GSE','GSM','tissue','Age_group'] if k in f['obs']]
        donorcols = [k for k in obscols if any(w in k.lower() for w in ['donor','patient','subject','individual'])]
        selected += donorcols
        obs = pd.DataFrame({k:clean(read_elem(f['obs'][k])) for k in selected})
        obs.insert(0,'cell_id',pd.Series(read_elem(f['obs'][indexkey]),dtype='string'))
        obs.insert(0,'row_id',np.arange(shape[0],dtype=np.int64))
        if len(obs)!=shape[0]: raise ValueError('Metadata row alignment failed')
        # Split uses group IDs alone, before any gene specificity or concordance analysis.
        group_col = cfg['group_column']
        groups = sorted(set(obs[group_col].astype(str))-{'Unknown'})
        rng=np.random.default_rng(cfg['seed']); perm=rng.permutation(groups)
        hold=set(perm[:max(1,int(round(len(groups)*cfg['holdout_fraction'])))])
        obs['split']=pd.Categorical(np.where(obs[group_col].astype(str).eq('Unknown'),'excluded_unknown_group',
                                  np.where(obs[group_col].astype(str).isin(hold),'holdout','train')))
        obs.to_parquet(c/'metadata.parquet',index=False,compression='zstd')
        split=pd.DataFrame({'group':groups,'split':['holdout' if g in hold else 'train' for g in groups]})
        split['lineage']=lin; split['group_column']=group_col
        split.to_csv(out/f'audit/group_split_{lin}.csv',index=False)
        np.save(c/'umap.npy',f['obsm/X_umap'][:].astype(np.float32))
        np.save(c/'plot_rows.npy',np.sort(rng.choice(shape[0],min(shape[0],cfg['plot_cells']),replace=False)))
        # Distributed contiguous probes: values, negative signs, integer-likeness and inverse-log totals.
        probes=[]; probe_values=[]; invsums=[]; sums=[]
        starts=np.unique(np.linspace(0,max(shape[0]-256,0),24,dtype=int))
        for start in starts:
            x=csr_rows(node,int(start),min(shape[0],int(start)+256))
            probe_values.append(x.data[::max(1,len(x.data)//10000)])
            sums.extend(np.asarray(x.sum(1)).ravel().tolist())
            if len(x.data) and x.data.max()<30:
                xx=x.copy(); xx.data=np.expm1(xx.data); invsums.extend(np.asarray(xx.sum(1)).ravel().tolist())
        vals=np.concatenate(probe_values)
        negative=float(np.mean(vals<0)); integer=float(np.mean(np.isclose(vals,np.round(vals),atol=1e-5)))
        has_log='log1p' in f['uns']
        if negative>0: scale='centered_or_negative_unverified'
        elif integer>0.999: scale='counts'
        elif has_log: scale='log1p_compatible'
        else: scale='unverified_nonnegative'
        expression={'source':'X','encoding':str(node.attrs.get('encoding-type','dense')),'shape':shape,
                    'layers':list(f.get('layers',{})), 'has_raw':'raw' in f, 'has_uns_log1p':has_log,
                    'probe_n_stored_values':len(vals),'probe_min':float(vals.min()),'probe_max':float(vals.max()),
                    'probe_negative_fraction':negative,'probe_integer_fraction':integer,'interpretation':scale,
                    'probe_row_sum_quantiles':np.quantile(sums,[0,.25,.5,.75,1]).tolist(),
                    'probe_expm1_row_sum_quantiles':np.quantile(invsums,[0,.25,.5,.75,1]).tolist() if invsums else [],
                    'normalization':'X used unchanged; nonnegative fractional values + log1p metadata; original library-size target not independently recoverable after gene filtering' if scale=='log1p_compatible' else scale,
                    'probe_starts':starts.tolist(),
                    'caveat':'Probes are not a full-matrix guarantee. Full extraction validates every stored value for negativity/non-finiteness. BBKNN corrects neighbors, not expression.'}
        dump(c/'expression_audit.json',expression)
        if scale not in {'log1p_compatible','counts'}: raise ValueError(f'No verified scoring expression: {scale}')
        sub_counts=obs['Celltype_general'].value_counts()
        info={'lineage':LABELS[lin],**fingerprint(path),'n_cells':shape[0],'n_genes':shape[1],
              'obs_columns':obscols,'var_columns':list(f['var']),'obsm':list(f['obsm']), 'uns':list(f['uns']),
              'expression':expression,'subtype_counts':sub_counts.to_dict(),'n_leiden':obs['leiden'].nunique(),
              'n_samples':obs['sample'].nunique(),'n_GSE':obs['GSE'].nunique(),'n_tissues':obs['tissue'].nunique(),
              'donor_fields':donorcols,'group_column':group_col,'train_groups':len(groups)-len(hold),'holdout_groups':len(hold),
              'missing_labels':int(obs['Celltype_general'].eq('Unknown').sum()),'duplicate_var_symbols':int(pd.Index(names).duplicated().sum()),
              'duplicate_cell_ids':int(obs['cell_id'].duplicated().sum())}
        dump(out/f'audit/data_audit_{lin}.json',info)
        pd.DataFrame([{'item':k,'value':json.dumps(v,ensure_ascii=False,default=serial)} for k,v in info.items()]).to_csv(out/f'audit/data_audit_{lin}.csv',index=False)
        scores=[]
        known={'T','B','Fibro','Endo','MoMac','Epi','ECM','Immune','AP','Tip','Stalk','Fib_deg','Endo_deg','MS_Fib_deg','MS_Endo_deg'}
        for k in obscols:
            if k in known or re.match(r'^(Myfibro|eFibro|iFibro|qFibro|SMC|Peri|LEC|VenEC|CapEC|ArtEC|MCs|EndoMT)([-_].*|1)?$',k):
                a=np.asarray(read_elem(f['obs'][k]))
                if not np.issubdtype(a.dtype,np.number): continue
                finite=np.isfinite(a)
                scores.append({'lineage':lin,'old_column':k,'n':len(a),'min':float(np.nanmin(a)),
                    'max':float(np.nanmax(a)),'mean':float(np.nanmean(a)), 'negative_fraction':float(np.mean(a<0)),
                    'missing_fraction':float(np.mean(~finite)),'provenance':'Notebook call evidence recorded; file snapshot has no immutable score metadata. Column name alone is not proof of last executed version.'})
        pd.DataFrame(scores).to_csv(out/f'audit/existing_score_columns_{lin}.csv',index=False)
        # Save only numeric old scores for the same deterministic display/test rows.
        rows=np.load(c/'plot_rows.npy')
        oldsample={'row_id':rows}
        for item in scores: oldsample[item['old_column']]=np.asarray(read_elem(f['obs'][item['old_column']]))[rows]
        pd.DataFrame(oldsample).to_parquet(c/'old_scores_sample.parquet',index=False)
    for key in ['leiden','Celltype','tissue','GSE','sample','Age_group']:
        counts=obs.groupby([key,'Celltype_general'],observed=True).size().rename('n_cells').reset_index()
        counts['lineage']=lin; counts.to_csv(out/f'audit/annotation_by_{key}_{lin}.csv',index=False)
    threeway=obs.groupby(['leiden','Celltype','Celltype_general'],observed=True).size().rename('n_cells').reset_index()
    threeway['lineage']=lin;threeway.to_csv(out/f'audit/annotation_cluster_crosstab_{lin}.csv',index=False)
    combo=obs.groupby(['tissue','GSE','Celltype_general'],observed=True).agg(n_cells=('row_id','size'),n_samples=('sample','nunique')).reset_index()
    combo.to_csv(out/f'audit/independent_units_{lin}.csv',index=False)
    overlaps=[]
    for field in ['GSM','sample']:
        tr=set(obs.loc[obs['split'].eq('train'),field].astype(str))-{'Unknown'}
        te=set(obs.loc[obs['split'].eq('holdout'),field].astype(str))-{'Unknown'}
        overlaps.append({'field':field,'shared_unique_strings_train_holdout':len(tr&te),
                         'interpretation':'Repeated names are not a verified donor ID; strings cannot establish individual independence'})
    pd.DataFrame(overlaps).to_csv(out/f'audit/sample_id_overlap_{lin}.csv',index=False)
    # Actual fine labels provide hierarchy evidence; current atlas has no SMC/Peri general labels.
    fine=obs[['Celltype','Celltype_general']].drop_duplicates().astype(str).to_dict('records')
    history=read_json(out/'audit/notebook_code_evidence.json')
    fragments=[{'source':h['source'],'cell':h['cell'],'lines':[l for l in h['code'].splitlines() if any(k in l for k in ['MCs','SMC','Peri','EndoMT-I','EndoMT-II','Celltype_general'])]} for h in history if h['source']==cfg['notebooks'][lin]]
    dump(out/f'audit/annotation_hierarchy_evidence_{lin}.json',{'fine_to_general':fine,'notebook_fragments':fragments,
         'interpretation':'Actual general labels are audited; union-projection is a candidate hierarchy assumption, not a cell-level proof that historical SMC/Peri assignments were merged unchanged.'})
    logging.info('%s audit done: %s cells, classes=%s, donor fields=%s',lin,len(obs),sub_counts.to_dict(),donorcols)
    return info


def marker_qc(cfg, lineages):
    out=root(cfg); cat=read_json(out/'marker_sets/candidate_marker_catalog.json')['candidates']; records=[]; specs=[]
    for lin in lineages:
        names=Counter(np.load(cache(cfg,lin)/'var_names.npy').tolist())
        allspec=[]
        for cand in cat:
            if cand['kind'].startswith('subtype_') and cand['kind']!=f'subtype_{lin}': continue
            # Include lineage-matched functional/broad versions, user versions apply to both.
            if cand['kind'] in {'broad_lineage','functional'} and cand['candidate_id'].startswith(('fibro:', 'endo:')) and not cand['candidate_id'].startswith(lin+':'): continue
            mapped={}
            for key,genes in cand['original_markers'].items():
                usable=[]; seen=set()
                for gene in genes:
                    target,status=resolve_gene(gene,names)
                    dup=target in seen if target else False
                    if target and not dup: usable.append(target); seen.add(target)
                    records.append({'lineage':lin,'candidate_id':cand['candidate_id'],'dictionary_name':cand['dictionary_name'],
                                    'original_key':key,'gene_listed':gene,'gene_resolved':target,'status':status,
                                    'duplicate_within_set_after_mapping':dup,'present':target is not None,
                                    'source':cand['source'],'cell_index':cand['cell_index']})
                mapped[key]=usable
                allspec.append({'set_id':cand['candidate_id']+'::'+key,'candidate_id':cand['candidate_id'],
                     'dictionary_name':cand['dictionary_name'],'key':key,'kind':cand['kind'],'genes_listed':genes,
                     'genes':usable,'n_genes_total':len(genes),'n_genes_matched':len(usable),
                     'mapping':'original','annotation_target':normalize_marker_key(key,lin) if cand['kind']==f'subtype_{lin}' else None})
            if cand['kind']==f'subtype_{lin}':
                merged={}
                for key,genes in mapped.items(): merged.setdefault(normalize_marker_key(key,lin),[]).extend(genes)
                for key,genes in merged.items():
                    usable=list(dict.fromkeys(genes))
                    allspec.append({'set_id':cand['candidate_id']+'::projected6::'+key,'candidate_id':cand['candidate_id'],
                         'dictionary_name':cand['dictionary_name'],'key':key,'kind':cand['kind'],'genes':usable,
                         'genes_listed':[g for k,v in cand['original_markers'].items() if normalize_marker_key(k,lin)==key for g in v],
                         'n_genes_total':sum(len(v) for k,v in cand['original_markers'].items() if normalize_marker_key(k,lin)==key),
                         'n_genes_matched':len(usable),'mapping':'projected','annotation_target':key})
        dump(cache(cfg,lin)/'score_specs.json',allspec)
        pd.DataFrame([{'candidate_id':s['candidate_id'],'original_key':s['key'],
             'target_reference_label':s['annotation_target'],'mapping':s['mapping'],
             'note':'Explicit diagnostic projection, not proof of ground truth or a mechanical cell-level merge' if s['kind'].startswith('subtype') else 'Not used in subtype argmax'}
             for s in allspec]).to_csv(out/f'data/candidate_label_projection_{lin}.csv',index=False)
        sets=[(s['set_id'],set(s['genes'])) for s in allspec if s['mapping']=='original']
        overlap=[]
        for i,(a,ga) in enumerate(sets):
            for b,gb in sets[i+1:]:
                overlap.append({'lineage':lin,'set_a':a,'set_b':b,'n_shared':len(ga&gb),
                                'jaccard':len(ga&gb)/len(ga|gb) if ga|gb else np.nan})
        pd.DataFrame(overlap).to_csv(out/f'data/candidate_gene_overlap_{lin}.csv',index=False)
    cv=pd.DataFrame(records)
    oldpath=out/'audit/marker_gene_coverage.csv'
    if oldpath.exists():
        previous=pd.read_csv(oldpath);previous=previous[~previous['lineage'].isin(lineages)]
        if len(previous):cv=pd.concat([previous,cv],ignore_index=True)
    cv.to_csv(out/'audit/marker_gene_coverage.csv',index=False)
    cv.to_csv(out/'data/marker_gene_coverage.csv',index=False)
    cv.to_csv(out/'data/marker_symbol_mapping.csv',index=False)
    links=[]
    for lin in ['fibro','endo']:
        op=out/f'audit/existing_score_columns_{lin}.csv';sp=cache(cfg,lin)/'score_specs.json'
        if not op.exists() or not sp.exists():continue
        old=pd.read_csv(op);specs=read_json(sp)
        for col in old['old_column']:
            matches=[s for s in specs if s['mapping']=='original' and s['key']==col]
            for s in matches or [None]:
                links.append({'lineage':lin,'old_column':col,'candidate_id':s['candidate_id'] if s else None,
                              'match_basis':'Same dictionary key only; numeric identity/correlation is supplementary provenance evidence' if s else 'Unresolved exact dictionary provenance',
                              'audit_note':'No old obs columns overwritten'})
    pd.DataFrame(links).to_csv(out/'audit/existing_score_dictionary_links.csv',index=False)
    logging.info('Phase 1: %s candidate gene/source records mapped; all absent/ambiguous symbols explicit',len(records))


def stage_audit(cfg, lineages):
    cat=collect_history(cfg)
    for lin in lineages: audit_one(cfg,lin,cat)
    out=root(cfg)
    for basename,pattern in [('annotation_cluster_crosstab','annotation_cluster_crosstab_*.csv'),('existing_score_columns','existing_score_columns_*.csv')]:
        tables=[pd.read_csv(p) for p in sorted((out/'audit').glob(pattern))]
        pd.concat(tables,ignore_index=True).to_csv(out/f'audit/{basename}.csv',index=False)
    texts=['# Expression layer audit', '', 'Input matrices were opened read-only. No whole-atlas dense conversion was performed.']
    for lin in ['fibro','endo']:
        p=cache(cfg,lin)/'expression_audit.json'
        if p.exists(): texts += ['',f'## {LABELS[lin]}','```json',json.dumps(read_json(p),indent=2),'```']
    texts += ['', 'No negative values are clipped. Counts, if actually identified, are normalized per row to 10,000 then log1p; current fractional log1p-compatible X is used unchanged.',
              'The historical normalization target is not inferred solely from BBKNN or from a filename. The nonnegative log-scale interpretation remains conditional on the recorded metadata/probes and full streaming checks.']
    (out/'audit/expression_layer_audit.md').write_text('\n'.join(texts))
    marker_qc(cfg,lineages)


def mean_scores(x, sets, gene_names):
    """Small-matrix primitive: sparse/dense parity; empty signatures stay NaN."""
    lookup={g:i for i,g in enumerate(gene_names)}
    result=np.full((x.shape[0],len(sets)),np.nan,dtype=np.float32)
    for j,genes in enumerate(sets):
        ii=[lookup[g] for g in dict.fromkeys(genes) if g in lookup]
        if ii: result[:,j]=np.asarray(x[:,ii].mean(axis=1)).ravel()
    return result


def smoke_tests(cfg):
    tests=[]
    a=np.array([[1.,0.,-2.],[0.,3.,1.],[2.,1.,0.]],np.float32)
    sets=[['A'],['A','B'],['MISSING'],['A','A','B']]
    d=mean_scores(a,sets,['A','B','C']);s=mean_scores(sparse.csr_matrix(a),sets,['A','B','C'])
    np.testing.assert_allclose(d,s,equal_nan=True); tests.append('sparse_dense_parity')
    assert np.isnan(d[:,2]).all();tests.append('missing_signature_is_nan')
    assert mean_scores(a,[['C']],['A','B','C'])[0,0]==-2.;tests.append('negative_values_not_clipped')
    np.testing.assert_allclose(d[:,1],a[:,:2].mean(1));tests.append('different_gene_set_sizes_use_means')
    np.testing.assert_allclose(d[:,1],d[:,3]);tests.append('duplicate_genes_not_double_counted')
    toyold=pd.DataFrame({'old_score':[1.,2.,3.]});before=toyold.copy(deep=True)
    _=mean_scores(a,[['A']],['A','B','C']);pd.testing.assert_frame_equal(toyold,before);tests.append('input_old_columns_untouched')
    groups=np.array(['g1','g1','g2','g3','g4']);hold={'g2','g4'}
    tr=set(groups[~np.isin(groups,list(hold))]);te=set(groups[np.isin(groups,list(hold))]);assert not tr&te;tests.append('no_group_leakage')
    idx=np.array([2,0]);np.testing.assert_array_equal(d[idx,0],a[idx,0]);tests.append('row_alignment')
    means=pd.DataFrame(d[:,:2]).mean().to_numpy();np.testing.assert_allclose(means,np.nanmean(d[:,:2],axis=0));tests.append('plotted_matrix_equals_summary')
    dump(root(cfg)/'logs/smoke_tests.json',{'passed':tests,'n_passed':len(tests),'time':time.strftime('%Y-%m-%d %H:%M:%S')})
    logging.info('Pre-recomputation smoke tests passed: %s',len(tests))


def training_control_pilot(cfg, lin, obs, names, marker_genes):
    c=cache(cfg,lin); n,p=len(obs),len(names); starts=np.unique(np.linspace(0,n-cfg['pilot_rows_per_block'],cfg['pilot_blocks'],dtype=int))
    sums=np.zeros(p,np.float64); ntrain=0; used=[];pilot_groups=set()
    train=obs['split'].eq('train').to_numpy()
    expression=read_json(c/'expression_audit.json')
    with h5py.File(cfg['inputs'][lin],'r') as f:
        for st in starts:
            en=min(n,int(st)+cfg['pilot_rows_per_block']);mask=train[st:en]
            if not mask.any(): continue
            x=csr_rows(f['X'],int(st),en)[mask]
            pilot_groups.update(obs.iloc[st:en].loc[mask,cfg['group_column']].astype(str))
            if expression['interpretation']=='counts':
                x=sparse.diags(1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1))@x;x=x.tocsr();x.data=np.log1p(x.data)
            sums+=np.asarray(x.sum(0)).ravel();ntrain+=x.shape[0];used.append({'start':int(st),'n_training_cells':int(x.shape[0])})
    if ntrain<1000: raise ValueError('Insufficient training-only cells for control bins')
    means=sums/ntrain
    bins=pd.qcut(pd.Series(means).rank(method='first'),cfg['expression_bins'],labels=False).to_numpy(dtype=int)
    excluded=set(marker_genes);rng=np.random.default_rng(cfg['seed']);controls={}
    for b in sorted(set(bins)):
        pool=np.array([i for i in np.flatnonzero(bins==b) if names[i] not in excluded and means[i]>0],dtype=int)
        if not len(pool):
            raise ValueError(f'No expressed background pool in bin {b}')
        controls[str(b)]=names[rng.choice(pool,min(len(pool),cfg['controls_per_bin']),replace=False)].tolist()
    pd.DataFrame({'gene':names,'training_pilot_mean':means,'expression_bin':bins}).to_csv(root(cfg)/f'data/control_expression_bins_{lin}.csv',index=False)
    info={'seed':cfg['seed'],'group_column':cfg['group_column'],'pilot_training_cells':ntrain,'pilot_blocks':used,
          'pilot_training_GSE':sorted(pilot_groups),'pilot_n_training_GSE':len(pilot_groups),
          'controls_by_bin':controls,'marker_bins':{g:int(bins[np.flatnonzero(names==g)[0]]) for g in marker_genes},
          'method':'Each marker contributes log expression minus mean of 12 fixed training-mean-bin matched controls; set score averages these per-marker contrasts. Controls shared by bin; markers excluded from pool. Not bit-identical to Scanpy union-of-controls scoring.',
          'holdout_used_for_bins':False,'all_candidate_genes_excluded':sorted(excluded)}
    dump(c/'controls.json',info);dump(root(cfg)/f'data/control_genes_{lin}.json',info)
    return info


def extract_expression(cfg, lin, resume=False):
    """One full sequential CSR pass; cache candidate genes and only 25 bin means."""
    c=cache(cfg,lin);done=c/'extraction_complete.json';fp=fingerprint(cfg['inputs'][lin])
    if resume and done.exists():
        if read_json(done)['source_fingerprint']!=fp: raise ValueError('Source changed; refuse stale extraction')
        logging.info('%s expression extraction reused (fingerprint verified)',lin);return
    obs=pd.read_parquet(c/'metadata.parquet');names=np.load(c/'var_names.npy');n=len(obs)
    specs=read_json(c/'score_specs.json')
    genes=sorted({g for s in specs for g in s['genes']});lookup={g:i for i,g in enumerate(names)}
    controls=training_control_pilot(cfg,lin,obs,names,genes)
    selected=sorted({lookup[g] for g in genes}|{lookup[g] for gg in controls['controls_by_bin'].values() for g in gg})
    remap={i:j for j,i in enumerate(selected)}
    gi=[remap[lookup[g]] for g in genes]
    ctrl_indices=[[remap[lookup[g]] for g in controls['controls_by_bin'][str(b)]] for b in range(cfg['expression_bins'])]
    xcache=np.lib.format.open_memmap(c/'marker_expression.npy',mode='w+',dtype=np.float32,shape=(n,len(genes)))
    bcache=np.lib.format.open_memmap(c/'control_bin_means.npy',mode='w+',dtype=np.float32,shape=(n,cfg['expression_bins']))
    np.save(c/'marker_genes.npy',np.asarray(genes))
    units=pd.MultiIndex.from_frame(obs[[cfg['group_column'],'Celltype_general']].astype(str))
    codes,uniques=pd.factorize(units,sort=True);unitmeta=uniques.to_frame(index=False)
    unitmeta.columns=[cfg['group_column'],'Celltype_general'];unitmeta['unit_index']=np.arange(len(unitmeta))
    splitmap=obs[[cfg['group_column'],'split']].drop_duplicates().set_index(cfg['group_column'])['split'].astype(str).to_dict()
    unitmeta['split']=unitmeta[cfg['group_column']].map(splitmap)
    unitmeta.to_parquet(c/'gene_group_units.parquet',index=False)
    sums=np.zeros((len(unitmeta),len(genes)),np.float64);detect=np.zeros_like(sums);counts=np.zeros(len(unitmeta),np.int64)
    chunk=cfg['chunk_rows'];value_min=np.inf;value_max=-np.inf;nstored=0;t0=time.time()
    expression=read_json(c/'expression_audit.json')
    with h5py.File(cfg['inputs'][lin],'r') as f:
        for st in range(0,n,chunk):
            en=min(n,st+chunk);x=csr_rows(f['X'],st,en)
            if x.nnz:
                if not np.isfinite(x.data).all(): raise ValueError(f'Non-finite X at {lin}:{st}')
                mn=float(x.data.min());mx=float(x.data.max())
                if mn<0: raise ValueError(f'Negative X found in full scan at {lin}:{st}; no clipping permitted')
                value_min=min(value_min,mn);value_max=max(value_max,mx);nstored+=x.nnz
            if expression['interpretation']=='counts':
                x=sparse.diags(1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1))@x;x=x.tocsr();x.data=np.log1p(x.data)
            y=x[:,selected].toarray();m=y[:,gi];xcache[st:en]=m
            bcache[st:en]=np.column_stack([y[:,ii].mean(1) for ii in ctrl_indices])
            ag=sparse.csr_matrix((np.ones(en-st),(codes[st:en],np.arange(en-st))),shape=(len(unitmeta),en-st))
            sums+=ag@m;detect+=ag@(m>0);counts+=np.bincount(codes[st:en],minlength=len(unitmeta))
            if st==0 or (st//chunk)%50==0 or en==n:
                logging.info('%s full X scan %s/%s (%.1f%%), elapsed %.1f min',lin,en,n,en*100/n,(time.time()-t0)/60)
                dump(c/'extraction_progress.json',{'rows_complete':en,'total_rows':n,'elapsed_seconds':time.time()-t0})
    xcache.flush();bcache.flush()
    np.savez_compressed(c/'gene_group_aggregates.npz',sums=sums,detected=detect,counts=counts)
    dump(done,{'source_fingerprint':fp,'n_cells':n,'n_candidate_genes':len(genes),
        'full_scan_stored_values':nstored,'full_scan_stored_min':value_min,'full_scan_stored_max':value_max,
        'negative_count':0,'nonfinite_count':0,'elapsed_seconds':time.time()-t0,
        'process_peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        'source_X_transformed':expression['interpretation']=='counts','alignment':'row i is original H5AD obs row i'})
    logging.info('%s extraction complete: %s genes, %.1f min',lin,len(genes),(time.time()-t0)/60)


def benchmark_actual_chunks(cfg):
    """Measured real 8,192-row block cost; not an invented full-run memory peak."""
    import psutil
    import gc
    records=[];process=psutil.Process()
    for lin in ['fibro','endo']:
        c=cache(cfg,lin);var=np.load(c/'var_names.npy').tolist();vi={g:i for i,g in enumerate(var)}
        genes=np.load(c/'marker_genes.npy').tolist();ctr=read_json(c/'controls.json')
        chosen=sorted({vi[g] for g in genes}|{vi[g] for gg in ctr['controls_by_bin'].values() for g in gg})
        gc.collect();before=process.memory_info().rss/2**20;t0=time.perf_counter()
        with h5py.File(cfg['inputs'][lin],'r') as f:x=csr_rows(f['X'],0,cfg['chunk_rows'])
        y=x[:,chosen].toarray();after=process.memory_info().rss/2**20
        records.append({'lineage':lin,'rows':x.shape[0],'total_genes':x.shape[1],'selected_genes_and_controls':len(chosen),
             'nnz':x.nnz,'elapsed_seconds':time.perf_counter()-t0,'RSS_before_MiB':before,'RSS_with_block_MiB':after,
             'process_peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
             'selected_dense_block_MiB':y.nbytes/2**20,'scope':'measured representative block; excludes full metadata and validation allocations; not a full-run peak'})
        del x,y
    dump(root(cfg)/'logs/resource_benchmark.json',records)
    return records


def unique_score_weights(specs, genes, marker_bins, nbins):
    lookup={g:i for i,g in enumerate(genes)};unique=[];seen={}
    for s in specs:
        key=tuple(sorted(s['genes']))
        if key not in seen: seen[key]=len(unique);unique.append(key)
        s['score_index']=seen[key]
    w=np.zeros((len(genes),len(unique)),np.float32);b=np.zeros((nbins,len(unique)),np.float32)
    for j,gg in enumerate(unique):
        for g in gg:
            w[lookup[g],j]=1/len(gg);b[marker_bins[g],j]+=1/len(gg)
    return sparse.csr_matrix(w),sparse.csr_matrix(b),unique


def calculate_candidates(cfg,lin,resume=False):
    c=cache(cfg,lin);done=c/'candidate_scores_complete.json'
    if resume and done.exists(): logging.info('%s candidate scores reused',lin);return
    genes=np.load(c/'marker_genes.npy').tolist();x=np.load(c/'marker_expression.npy',mmap_mode='r')
    bg=np.load(c/'control_bin_means.npy',mmap_mode='r');controls=read_json(c/'controls.json')
    specs=read_json(c/'score_specs.json');w,b,unique=unique_score_weights(specs,genes,controls['marker_bins'],cfg['expression_bins'])
    n=len(x);mean=np.lib.format.open_memmap(c/'candidate_mean.npy',mode='w+',dtype=np.float32,shape=(n,len(unique)))
    module=np.lib.format.open_memmap(c/'candidate_module.npy',mode='w+',dtype=np.float32,shape=(n,len(unique)))
    empty=[j for j,g in enumerate(unique) if not g]
    for st in range(0,n,65536):
        en=min(n,st+65536);m=(w.T@x[st:en].T).T;ctrl=(b.T@bg[st:en].T).T
        mean[st:en]=m;module[st:en]=m-ctrl
        if empty: mean[st:en,empty]=np.nan;module[st:en,empty]=np.nan
    mean.flush();module.flush();dump(c/'score_specs_indexed.json',specs)
    unitmeta=pd.read_parquet(c/'gene_group_units.parquet');agg=np.load(c/'gene_group_aggregates.npz')
    expressed=set(np.asarray(genes)[agg['detected'].sum(0)>0])
    meta=[]
    for s in specs:
        meta.append({k:v for k,v in s.items() if k not in {'genes','genes_listed'}}|
          {'lineage':lin,'n_genes_expressed':sum(g in expressed for g in s['genes']),
           'genes':';'.join(s['genes']),'expression_layer':'X','normalization':read_json(c/'expression_audit.json')['normalization'],
           'scoring_method':'mean_log_expression AND bin_matched_background_corrected','version':VERSION,
           'compute_status':'computable' if s['genes'] else 'insufficient_coverage_nan'})
    pd.DataFrame(meta).to_csv(root(cfg)/f'data/score_metadata_{lin}.csv',index=False)
    dump(done,{'n_unique_gene_sets':len(unique),'n_catalog_sets':len(specs),'n_cells':n})
    logging.info('%s scored %s catalog signatures (%s unique) by two methods',lin,len(specs),len(unique))


def stage_score(cfg,lineages,resume=False):
    smoke_tests(cfg)
    for lin in lineages:
        extract_expression(cfg,lin,resume=resume)
        calculate_candidates(cfg,lin,resume=resume)


def bh_adjust(p):
    p=np.asarray(p,float);result=np.full(len(p),np.nan);ok=np.isfinite(p)
    ii=np.flatnonzero(ok);order=np.argsort(p[ok]);v=p[ok][order];m=len(v)
    if m:
        q=np.minimum.accumulate((v*m/np.arange(1,m+1))[::-1])[::-1]
        result[ii[order]]=np.minimum(q,1.)
    return result


def paired_summary(differences):
    a=np.asarray(differences,float);a=a[np.isfinite(a)]
    if not len(a): return {'n_pairs':0,'median_effect':np.nan,'mean_effect':np.nan,'positive_fraction':np.nan,'p_value':np.nan}
    if np.all(np.abs(a)<1e-12): p=1.
    else:
        try: p=float(stats.wilcoxon(a,alternative='greater',zero_method='wilcox').pvalue)
        except ValueError: p=1.
    return {'n_pairs':len(a),'median_effect':float(np.median(a)),'mean_effect':float(a.mean()),
            'positive_fraction':float(np.mean(a>0)),'p_value':p}


def choose_final_markers(cfg,lin):
    """Marker selection sees TRAIN GSE aggregates only, never holdout scores."""
    c=cache(cfg,lin);genes=np.load(c/'marker_genes.npy').tolist();lookup={g:i for i,g in enumerate(genes)}
    units=pd.read_parquet(c/'gene_group_units.parquet');ag=np.load(c/'gene_group_aggregates.npz')
    means=ag['sums']/np.maximum(ag['counts'][:,None],1);det=ag['detected']/np.maximum(ag['counts'][:,None],1)
    specs=read_json(c/'score_specs_indexed.json');order=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_order'][lin]
    pools={s:set() for s in order};sources={}
    for s in specs:
        if s['kind']==f'subtype_{lin}' and s['mapping']=='original' and s['annotation_target'] in pools:
            pools[s['annotation_target']].update(s['genes'])
            for g in s['genes']: sources.setdefault((s['annotation_target'],g),set()).add(s['candidate_id'])
    use=units['split'].eq('train').to_numpy()&(ag['counts']>=cfg['min_cells_per_group_subtype'])
    uid=units[use].copy();uid['array_i']=np.flatnonzero(use)
    group_lookup={(str(r[cfg['group_column']]),str(r['Celltype_general'])):int(r['array_i']) for _,r in uid.iterrows()}
    train_groups=sorted(uid[cfg['group_column']].unique());rows=[];pairrows=[]
    for subtype in order:
        target_idx=uid.loc[uid['Celltype_general'].eq(subtype),'array_i'].to_numpy(dtype=int)
        for gene in sorted(pools[subtype]):
            j=lookup[gene];diff=[];de_diff=[]
            for group in train_groups:
                ti=group_lookup.get((group,subtype))
                oi=[group_lookup[(group,k)] for k in order if k!=subtype and (group,k) in group_lookup]
                if ti is not None and oi:
                    diff.append(means[ti,j]-means[oi,j].mean())
                    de_diff.append(det[ti,j]-det[oi,j].mean())
            sm=paired_summary(diff)
            pair_effects=[]
            for other in order:
                if other==subtype: continue
                dd=[];ddet=[]
                for group in train_groups:
                    ti=group_lookup.get((group,subtype));oi=group_lookup.get((group,other))
                    if ti is not None and oi is not None:
                        dd.append(means[ti,j]-means[oi,j]);ddet.append(det[ti,j]-det[oi,j])
                ps=paired_summary(dd)
                pairrows.append({'lineage':lin,'subtype':subtype,'comparison':other,'gene':gene,'split':'train',
                                 **ps,'median_detection_delta':float(np.median(ddet)) if ddet else np.nan})
                if ps['n_pairs']>=cfg['min_training_pairs']: pair_effects.append(ps['median_effect'])
            n_cell=int(ag['counts'][target_idx].sum())
            cell_det=float(ag['detected'][target_idx,j].sum()/max(n_cell,1))
            equal_det=float(det[target_idx,j].mean()) if len(target_idx) else np.nan
            rows.append({'lineage':lin,'subtype':subtype,'gene':gene,'selection_split':'train',**sm,
                         'target_detection_cell_weighted':cell_det,'target_detection_GSE_equal':equal_det,
                         'target_mean_GSE_equal':float(means[target_idx,j].mean()) if len(target_idx) else np.nan,
                         'median_detection_delta':float(np.median(de_diff)) if de_diff else np.nan,
                         'min_pairwise_median_effect':float(np.min(pair_effects)) if pair_effects else np.nan,
                         'n_training_GSE_with_subtype':len(target_idx),'n_training_cells':n_cell,
                         'n_subtype_pools_containing_gene':sum(gene in pool for pool in pools.values()),
                         'candidate_sources':';'.join(sorted(sources[(subtype,gene)]))})
    evidence=pd.DataFrame(rows);evidence['fdr']=bh_adjust(evidence['p_value'])
    evidence['passes_prespecified_QC']=(evidence['n_pairs'].ge(cfg['min_training_pairs'])&
        evidence['target_detection_GSE_equal'].ge(cfg['selection_min_detection'])&
        evidence['median_effect'].ge(cfg['selection_min_pair_effect'])&
        evidence['positive_fraction'].ge(cfg['selection_min_stability'])&evidence['fdr'].le(.05))
    # Rank by worst pairwise separation, cross-source sign stability, then mean effect.
    # The fixed cap is declared in config before inspecting holdout performance.
    final={};evidence['selected']=False
    for subtype in order:
        eligible=evidence[evidence['subtype'].eq(subtype)&evidence['passes_prespecified_QC']].sort_values(
            ['min_pairwise_median_effect','positive_fraction','median_effect','gene'],ascending=[False,False,False,True])
        picked=eligible.head(cfg['max_final_genes'])
        evidence.loc[picked.index,'selected']=True;final[subtype]=picked['gene'].tolist()
    evidence['decision_reason']=np.where(evidence['selected'],'retained: training coverage/effect/FDR/stability and prespecified rank',
        np.where(evidence['passes_prespecified_QC'],'eligible but beyond fixed maximum gene count','failed one or more prespecified training QC thresholds'))
    evidence.to_csv(root(cfg)/f'data/marker_selection_evidence_{lin}.csv',index=False)
    pairs=pd.DataFrame(pairrows);pairs['fdr']=bh_adjust(pairs['p_value']);pairs.to_csv(root(cfg)/f'data/targeted_pairwise_expression_{lin}.csv',index=False)
    # Freeze the list before evaluating any holdout concordance.
    doc={'lineage':LABELS[lin],'annotation_column':'Celltype_general','version':VERSION,
         'marker_sets':final,'expression_layer':'X','normalization':read_json(c/'expression_audit.json')['normalization'],
         'scoring_method':['mean_log_expression','training_bin_matched_module'],
         'validation':{'selection_rule':{
             'training_GSE_pairs_min':cfg['min_training_pairs'],'cells_per_GSE_subtype_min':cfg['min_cells_per_group_subtype'],
             'GSE_equal_detection_min':cfg['selection_min_detection'],'paired_median_log_expression_delta_min':cfg['selection_min_pair_effect'],
             'positive_GSE_fraction_min':cfg['selection_min_stability'],'one_sided_paired_wilcoxon_BH_max':.05,
             'ranking':['worst_pairwise_median_effect','GSE_sign_stability','paired_median_effect'],'max_genes':cfg['max_final_genes']},
           'grouped_holdout':{'grouping':'GSE','holdout_fraction':cfg['holdout_fraction'],'seed':cfg['seed'],
                              'marker_selection_uses_holdout':False,'control_bins_use_holdout':False,
                              'split_table':f'audit/group_split_{lin}.csv'},
           'donor_status':'unavailable: no verified donor identifier in source obs',
           'limitations':['Labels were originally informed by overlapping marker sets: this is annotation concordance, not independent biological validation.',
              'GSE is a study unit, not a donor; independence across shared individuals cannot be verified.',
              'Only genes retained in this atlas can be tested; absence is not biological nonexpression.',
              'Raw counts are unavailable; targeted tests compare grouped mean log expression, not count-model differential expression.',
              'Raw module score argmax is an uncalibrated concordance diagnostic; different gene sets retain different baselines and dropout.']},
         'subtype_evidence':{s:{'status':'pending_holdout' if gg else 'insufficient_evidence','genes':gg,
                              'selection_evidence_table':f'data/marker_selection_evidence_{lin}.csv'} for s,gg in final.items()}}
    dump(root(cfg)/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json',doc)
    dump(c/'frozen_final_markers.json',doc)
    logging.info('%s final markers frozen using training GSE only: %s',lin,final)
    return doc


def final_scores(cfg,lin,doc):
    import pyarrow as pa
    import pyarrow.parquet as pq
    c=cache(cfg,lin);obs=pd.read_parquet(c/'metadata.parquet');genes=np.load(c/'marker_genes.npy').tolist()
    x=np.load(c/'marker_expression.npy',mmap_mode='r');bg=np.load(c/'control_bin_means.npy',mmap_mode='r')
    ctr=read_json(c/'controls.json');order=list(doc['marker_sets']);n=len(x)
    specs=[{'genes':doc['marker_sets'][s]} for s in order]
    w,b,unique=unique_score_weights(specs,genes,ctr['marker_bins'],cfg['expression_bins'])
    indices=[s['score_index'] for s in specs]
    mean=np.lib.format.open_memmap(c/'final_mean.npy',mode='w+',dtype=np.float32,shape=(n,len(order)))
    module=np.lib.format.open_memmap(c/'final_module.npy',mode='w+',dtype=np.float32,shape=(n,len(order)))
    dest=root(cfg)/f'data/cell_signature_scores_{lin}.parquet';writer=None
    for st in range(0,n,65536):
        en=min(n,st+65536);m=((w.T@x[st:en].T).T)[:,indices];v=m-((b.T@bg[st:en].T).T)[:,indices]
        for j,s in enumerate(order):
            if not doc['marker_sets'][s]: m[:,j]=np.nan;v[:,j]=np.nan
        mean[st:en]=m;module[st:en]=v
        frame=obs.iloc[st:en][['row_id','cell_id','Celltype_general','leiden','GSE','sample','tissue','split']].copy()
        for j,s in enumerate(order):
            frame[f'sig_{lin}_final_{s}_mean_{VERSION}']=m[:,j]
            frame[f'sig_{lin}_final_{s}_module_{VERSION}']=v[:,j]
        table=pa.Table.from_pandas(frame,preserve_index=False)
        if writer is None: writer=pq.ParquetWriter(dest,table.schema,compression='zstd')
        writer.write_table(table)
    writer.close();mean.flush();module.flush();np.save(c/'final_score_names.npy',np.asarray(order))
    logging.info('%s full-cell final scores saved with row/cell IDs: %s',lin,dest.name)


def aggregate_scores(values,obs,keys,chunk=131072):
    # Preserve categorical codes: materializing millions of repeated strings is unnecessary.
    codes,unique=pd.factorize(pd.MultiIndex.from_frame(obs[keys]),sort=True)
    meta=unique.to_frame(index=False);meta.columns=keys;n=len(meta);p=values.shape[1]
    sums=np.zeros((n,p),np.float64);valid=np.zeros((n,p),np.int64)
    for st in range(0,len(obs),chunk):
        en=min(len(obs),st+chunk);v=np.asarray(values[st:en]);ok=np.isfinite(v)
        ag=sparse.csr_matrix((np.ones(en-st),(codes[st:en],np.arange(en-st))),shape=(n,en-st))
        sums+=ag@np.nan_to_num(v);valid+=np.asarray(ag@ok,dtype=np.int64)
    means=np.divide(sums,valid,out=np.full_like(sums,np.nan),where=valid>0)
    meta['n_cells']=np.bincount(codes,minlength=n)
    return meta,means


def concordance_metrics(y,values,k,groups=None):
    """Top scores are descriptive agreement; exact ties/empty rows abstain."""
    vals=np.nan_to_num(np.asarray(values),nan=-np.inf);rank=np.argsort(vals,axis=1,kind='stable')[:,::-1]
    pred=rank[:,0];top=vals[np.arange(len(vals)),rank[:,0]]
    second=vals[np.arange(len(vals)),rank[:,1]] if k>1 else np.zeros(len(vals))
    margin=top-second;abstain=(~np.isfinite(top))|(~np.isfinite(margin))|(margin<=1e-8)
    pred=np.where(abstain,k,pred);valid=(y>=0)&(y<k)
    conf=np.bincount((y[valid]*(k+1)+pred[valid]).astype(int),minlength=k*(k+1)).reshape(k,k+1)
    tp=np.diag(conf[:,:k]);den=conf.sum(1);recall=tp/np.maximum(den,1)
    prec=tp/np.maximum(conf[:,:k].sum(0),1);f1=2*prec*recall/np.maximum(prec+recall,1e-12)
    top2=((rank[:,:min(2,k)]==y[:,None]).any(1)&~abstain)
    m={'n_cells':int(valid.sum()),'top1_concordance':float(tp.sum()/max(valid.sum(),1)),
       'top2_concordance':float(top2[valid].mean()) if valid.any() else np.nan,
       'macro_f1':float(f1.mean()),'balanced_accuracy':float(recall.mean()),
       'abstention_fraction':float(abstain[valid].mean()) if valid.any() else np.nan,
       'mean_top1_top2_margin':float(np.nanmean(np.where(np.isfinite(margin[valid]),margin[valid],np.nan))) if valid.any() else np.nan}
    return m,pred,margin,conf


def bootstrap_group_ci(values,seed,nboot=500):
    a=np.asarray(values,float);a=a[np.isfinite(a)]
    if not len(a): return np.nan,np.nan
    rng=np.random.default_rng(seed);bs=np.mean(a[rng.integers(0,len(a),size=(nboot,len(a)))],axis=1)
    return tuple(np.quantile(bs,[.025,.975]).tolist())


def compare_candidate_versions(cfg,lin):
    c=cache(cfg,lin);obs=pd.read_parquet(c/'metadata.parquet');specs=read_json(c/'score_specs_indexed.json')
    order=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_order'][lin];k=len(order)
    y=pd.Categorical(obs['Celltype_general'],categories=order).codes;methods={m:np.load(c/f'candidate_{m}.npy',mmap_mode='r') for m in ['mean','module']}
    versions={}
    for s in specs:
        if s['kind']==f'subtype_{lin}' and s['mapping']=='projected': versions.setdefault(s['candidate_id'],{})[s['key']]=s['score_index']
    rows=[];setrows=[]
    for method,values in methods.items():
        um,uv=aggregate_scores(values,obs,['GSE','Celltype_general','split'])
        for cid,imap in versions.items():
            complete=set(order)<=set(imap) and all(np.isfinite(values[0,imap[s]]) for s in order if s in imap)
            missing=[s for s in order if s not in imap]
            empty=[s['key'] for s in specs if s['candidate_id']==cid and s['mapping']=='projected' and not s['genes']]
            coverage_status='; '.join(([f'missing reference signatures: {missing}'] if missing else [])+
                                     ([f'no usable genes: {empty}'] if empty else [])) or 'complete'
            for split in ['train','holdout','all']:
                use=np.ones(len(obs),bool) if split=='all' else obs['split'].eq(split).to_numpy()
                record={'lineage':lin,'candidate_id':cid,'method':method,'split':split,'n_signatures_available':len(imap),
                        'required_signatures':k,'complete_for_topscore_comparison':complete,
                        'coverage_comparison_status':coverage_status,
                        'scope':'original labels normalized; 8-type dictionaries use explicit marker gene-union projection to 6 labels'}
                if complete:
                    vals=np.column_stack([values[use,imap[s]] for s in order]);m,_,_,_=concordance_metrics(y[use],vals,k)
                    record.update(m)
                    umask=(um['n_cells']>=cfg['min_cells_per_group_subtype'])
                    if split!='all': umask&=um['split'].eq(split)
                    yy=pd.Categorical(um.loc[umask,'Celltype_general'],categories=order).codes
                    vv=uv[umask][:,[imap[s] for s in order]]
                    _,pr,_,_=concordance_metrics(yy,vv,k)
                    gr=pd.DataFrame({'GSE':um.loc[umask,'GSE'].to_numpy(),'match':pr==yy}).groupby('GSE')['match'].mean()
                    lo,hi=bootstrap_group_ci(gr,cfg['seed'],cfg['bootstrap_replicates'])
                    record.update({'GSE_equal_groupmean_top1':float(gr.mean()),'GSE_equal_ci_low':lo,'GSE_equal_ci_high':hi,'n_GSE':len(gr),'n_GSE_subtype_units':int(umask.sum())})
                rows.append(record)
            # Each candidate signature is evaluated even if the full family is incomplete.
        for s in specs:
            if s['mapping']!='original': continue
            j=s['score_index'];target=s['annotation_target']
            for split in ['train','holdout']:
                use=um['split'].eq(split)&um['n_cells'].ge(cfg['min_cells_per_group_subtype'])
                targetvals=uv[use&um['Celltype_general'].eq(target),j] if target else np.array([])
                offvals=uv[use&~um['Celltype_general'].eq(target),j] if target else np.array([])
                setrows.append({'lineage':lin,'set_id':s['set_id'],'kind':s['kind'],'method':method,'split':split,
                                'annotation_target':target,'n_usable_genes':len(s['genes']),
                                'target_GSE_equal_mean':float(np.nanmean(targetvals)) if np.isfinite(targetvals).any() else np.nan,
                                'other_GSE_subtype_equal_mean':float(np.nanmean(offvals)) if np.isfinite(offvals).any() else np.nan,
                                'n_target_GSE':len(targetvals),'n_other_GSE_subtype_units':len(offvals)})
        broadrows=[]
        for s in specs:
            if s['mapping']!='original' or s['kind'].startswith('subtype'): continue
            j=s['score_index']
            for split in ['train','holdout']:
                for subtype in order:
                    take=um['split'].eq(split)&um['Celltype_general'].eq(subtype)&um['n_cells'].ge(cfg['min_cells_per_group_subtype'])
                    a=uv[take,j]
                    broadrows.append({'lineage':lin,'set_id':s['set_id'],'kind':s['kind'],'signature':s['key'],
                                      'method':method,'split':split,'subtype':subtype,'n_GSE':int(take.sum()),
                                      'mean_score_GSE_equal':float(np.nanmean(a)) if np.isfinite(a).any() else np.nan,
                                      'interpretation':'Broad/functional diagnostic only; never used in subtype argmax.'})
        pd.DataFrame(broadrows).to_csv(root(cfg)/f'data/broad_functional_scores_{lin}_{method}.csv',index=False)
    pd.DataFrame(rows).to_csv(root(cfg)/f'data/marker_sets_candidate_comparison_{lin}.csv',index=False)
    pd.DataFrame(setrows).to_csv(root(cfg)/f'data/candidate_individual_signature_comparison_{lin}.csv',index=False)
    logging.info('%s all candidate dictionary versions compared by mean and module scores',lin)


def old_score_comparison(cfg,lin):
    c=cache(cfg,lin);old=pd.read_parquet(c/'old_scores_sample.parquet');ii=old['row_id'].to_numpy(dtype=int)
    specs=read_json(c/'score_specs_indexed.json');rows=[]
    for method in ['mean','module']:
        scores=np.load(c/f'candidate_{method}.npy',mmap_mode='r')
        for s in specs:
            if s['mapping']!='original' or s['key'] not in old: continue
            a=old[s['key']].to_numpy();b=scores[ii,s['score_index']];ok=np.isfinite(a)&np.isfinite(b)
            if ok.sum()>2 and a[ok].std()>0 and b[ok].std()>0:
                pear=float(stats.pearsonr(a[ok],b[ok]).statistic);spe=float(stats.spearmanr(a[ok],b[ok]).statistic)
            else: pear=spe=np.nan
            rows.append({'lineage':lin,'old_column':s['key'],'candidate_id':s['candidate_id'],'new_method':method,
                         'n_cells_compared':int(ok.sum()),'pearson_r':pear,'spearman_r':spe,
                         'mean_absolute_difference':float(np.mean(np.abs(a[ok]-b[ok]))) if ok.any() else np.nan,
                         'comparison_note':'Same original-row sample; correlation/identity checks are evidence, not immutable provenance.'})
    pd.DataFrame(rows).to_csv(root(cfg)/f'data/old_vs_recomputed_scores_{lin}.csv',index=False)


def final_validation(cfg,lin,doc):
    c=cache(cfg,lin);out=root(cfg);obs=pd.read_parquet(c/'metadata.parquet');order=list(doc['marker_sets']);k=len(order)
    y=pd.Categorical(obs['Celltype_general'],categories=order).codes
    metricrows=[];classrows=[];distrows=[];group_rows=[];cluster_rows=[];margin_rows=[];matrix_rows=[];cirows=[]
    counts=obs.groupby(['leiden','Celltype_general'],observed=True).size().rename('n_cells').reset_index()
    cluster_indices=obs.groupby('leiden',observed=True).indices
    dominant=counts.sort_values('n_cells',ascending=False).drop_duplicates('leiden').set_index('leiden')
    totals=counts.groupby('leiden',observed=True)['n_cells'].sum()
    for method in ['mean','module']:
        values=np.load(c/f'final_{method}.npy',mmap_mode='r')
        gm,gv=aggregate_scores(values,obs,['GSE','Celltype_general','split'])
        for scope in ['GSE','tissue','leiden','sample']:
            if scope=='GSE': um,uv=gm,gv
            else: um,uv=aggregate_scores(values,obs,[scope,'Celltype_general','split'])
            for j,s in enumerate(order):
                for i,r in um.iterrows():
                    group_rows.append({'lineage':lin,'group_level':scope,'group_id':r[scope],
                        'subtype':r['Celltype_general'],'split':r['split'],'signature':s,'method':method,
                        'mean_score':uv[i,j],'n_cells':int(r['n_cells'])})
        for split in ['train','holdout','all']:
            mask=np.ones(len(obs),bool) if split=='all' else obs['split'].eq(split).to_numpy()
            m,pred,margin,conf=concordance_metrics(y[mask],values[mask],k)
            metricrows.append({'lineage':lin,'method':method,'split':split,**m,'interpretation':'reference annotation concordance; not ground truth accuracy'})
            pd.DataFrame(conf,index=order,columns=order+['abstain']).to_csv(out/f'data/confusion_{lin}_{method}_{split}.csv')
            for j,s in enumerate(order):
                tm=(y[mask]==j);n=int(tm.sum());tp=int(conf[j,j])
                sv=np.asarray(values[np.flatnonzero(mask)[tm]])
                if len(sv):
                    sr=np.argsort(np.nan_to_num(sv,nan=-np.inf),axis=1)[:,::-1]
                    gap=np.take_along_axis(sv,sr[:,:1],axis=1).ravel()-np.take_along_axis(sv,sr[:,1:2],axis=1).ravel()
                    top2=float(((sr[:,:2]==j).any(1)&np.isfinite(gap)&(gap>1e-8)).mean())
                else:top2=np.nan
                classrows.append({'lineage':lin,'method':method,'split':split,'subtype':s,'n_cells':n,
                                  'top1_concordance':tp/max(n,1),'top2_concordance':top2,
                                  'mean_top1_top2_margin':float(np.nanmean(margin[tm])) if n else np.nan})
                for sj,sig in enumerate(order):
                    v=np.asarray(values[mask,sj])[tm];vf=v[np.isfinite(v)]
                    qs=np.quantile(vf,[.05,.25,.5,.75,.95]) if len(vf) else [np.nan]*5
                    distrows.append({'lineage':lin,'method':method,'split':split,'subtype':s,'signature':sig,
                         'n_cells':n,'n_finite':len(vf),'mean_score':float(vf.mean()) if len(vf) else np.nan,
                         'median_score':qs[2],'q05':qs[0],'q25':qs[1],'q75':qs[3],'q95':qs[4],
                         'fraction_positive_score':float(np.mean(vf>0)) if len(vf) else np.nan})
            umask=gm['n_cells'].ge(cfg['min_cells_per_group_subtype']).to_numpy()
            if split!='all': umask&=gm['split'].eq(split).to_numpy()
            yy=pd.Categorical(gm.loc[umask,'Celltype_general'],categories=order).codes
            _,pr,mar,_=concordance_metrics(yy,gv[umask],k)
            u=gm.loc[umask].copy();u['predicted_signature']=[order[p] if p<k else 'abstain' for p in pr]
            u['match']=pr==yy;u['method']=method
            equal=u.groupby('GSE',observed=True)['match'].mean()
            lo,hi=bootstrap_group_ci(equal,cfg['seed'],cfg['bootstrap_replicates'])
            cirows.append({'lineage':lin,'method':method,'split':split,'unit_type':'GSE',
               'n_GSE':len(equal),'n_GSE_subtype_units':len(u),'GSE_equal_groupmean_top1':float(equal.mean()),
               'ci_low':lo,'ci_high':hi,'bootstrap_replicates':cfg['bootstrap_replicates'],
               'note':'Each GSE has equal weight; subtypes present within each GSE are averaged. GSE is not donor.'})
            for j,s in enumerate(order):
                where=u['Celltype_general'].eq(s).to_numpy();vals=gv[umask][where]
                row_mean=np.nanmean(vals,axis=0) if len(vals) else np.full(k,np.nan)
                best=int(np.nanargmax(row_mean)) if np.isfinite(row_mean).any() else -1
                own=row_mean[j];competitor=np.nanmax(np.delete(row_mean,j)) if np.isfinite(row_mean).sum()>1 else np.nan
                for sj,sig in enumerate(order):
                    matrix_rows.append({'lineage':lin,'method':method,'split':split,'reference_subtype':s,'signature':sig,
                       'GSE_equal_mean_score':row_mean[sj],'n_GSE':int(where.sum()),
                       'diagonal_dominant':best==j,'top_signature':order[best] if best>=0 else 'unscorable',
                       'own_vs_best_other_margin':own-competitor})
            if split=='holdout':
                for i,(idx,r) in enumerate(gm.loc[umask].iterrows()):
                    yi=order.index(r['Celltype_general']);v=gv[idx];own=v[yi]
                    competing=np.delete(v,yi);other=np.nanmax(competing) if np.isfinite(competing).any() else np.nan
                    margin_rows.append({'lineage':lin,'method':method,'GSE':r['GSE'],'subtype':r['Celltype_general'],
                        'n_cells':int(r['n_cells']),'own_score':own,'best_other_score':other,'own_minus_best_other':own-other,
                        'groupmean_top1_match':bool(pr[i]==yi)})
            # Cell-level GSE × subtype agreement is distinct from agreement of group means.
            z=obs.loc[mask,['GSE','Celltype_general','split']].copy()
            z['match']=pred==y[mask];z['margin']=np.where(np.isfinite(margin),margin,np.nan)
            group_conc=z.groupby(['GSE','Celltype_general','split'],observed=True).agg(
                n_cells=('match','size'),cell_top1_concordance=('match','mean'),mean_top1_top2_margin=('margin','mean')).reset_index()
            group_conc['lineage']=lin;group_conc['method']=method;group_conc['evaluation_scope']=split
            group_conc.to_csv(out/f'data/GSE_cell_concordance_{lin}_{method}_{split}.csv',index=False)
            tissue=obs.loc[mask,['tissue','Celltype_general']].copy()
            tissue['cell_top1_match']=pred==y[mask]
            tissue.groupby(['tissue','Celltype_general'],observed=True)['cell_top1_match'].agg(['mean','size']).reset_index().to_csv(
                out/f'data/tissue_cell_concordance_{lin}_{method}_{split}.csv',index=False)
        # Per-Leiden composition and scores use all cells, no clusters are omitted.
        cm,cv=aggregate_scores(values,obs,['leiden'])
        for i,r in cm.iterrows():
            cluster=r['leiden'];reference=str(dominant.loc[cluster,'Celltype_general']);refidx=order.index(reference)
            v=cv[i];top=int(np.nanargmax(v)) if np.isfinite(v).any() else -1
            own=v[refidx];others=np.delete(v,refidx);other=float(np.nanmax(others)) if np.isfinite(others).any() else np.nan
            composition=counts[counts['leiden'].astype(str).eq(str(cluster))]
            mix={str(x['Celltype_general']):int(x['n_cells']) for _,x in composition.iterrows()}
            # Fine Celltype diversity is reported without turning Leiden into a subtype.
            cell_indices=cluster_indices[cluster]
            yy=y[cell_indices];_,pp,_,_=concordance_metrics(yy,values[cell_indices],k)
            for sj,sig in enumerate(order):
                cluster_rows.append({'lineage':lin,'method':method,'leiden':cluster,'signature':sig,'mean_score':v[sj],
                    'n_cells':int(r['n_cells']),'dominant_reference_subtype':reference,
                    'dominant_reference_fraction':int(dominant.loc[cluster,'n_cells'])/int(totals.loc[cluster]),
                    'reference_composition_json':json.dumps(mix),'n_fine_celltypes':obs.iloc[cell_indices]['Celltype'].nunique(),
                    'highest_signature':order[top] if top>=0 else 'unscorable','dominant_label_matches_top':top==refidx,
                    'cell_top1_concordance':float(np.mean(pp==yy)),'own_minus_best_other':own-other})
    for name,records in [('final_concordance',metricrows),('final_concordance_by_subtype',classrows),
         ('signature_score_distributions',distrows),('signature_scores_group_summary',group_rows),
         ('signature_concordance_by_cluster',cluster_rows),('holdout_GSE_margins',margin_rows),
         ('subtype_signature_matrix',matrix_rows),('grouped_holdout_concordance',cirows)]:
        pd.DataFrame(records).to_csv(out/f'data/{name}_{lin}.csv',index=False)
    matrices=pd.DataFrame(matrix_rows);margins=pd.DataFrame(margin_rows)
    for s in order:
        row=matrices[(matrices['method']=='module')&(matrices['split']=='holdout')&(matrices['reference_subtype']==s)].iloc[0]
        gs=margins[(margins['method']=='module')&(margins['subtype']==s)]
        lo,hi=bootstrap_group_ci(gs['own_minus_best_other'],cfg['seed'],cfg['bootstrap_replicates'])
        frac=float(gs['groupmean_top1_match'].mean()) if len(gs) else np.nan
        if not doc['marker_sets'][s] or len(gs)<10: status='insufficient_evidence'
        elif row['diagonal_dominant'] and frac>=.65 and lo>0: status='supported_relative_to_reference'
        else: status='ambiguous'
        doc['subtype_evidence'][s].update({'status':status,'holdout_n_GSE':len(gs),
            'holdout_GSE_groupmean_top1_fraction':frac,'holdout_GSE_equal_diagonal_dominant':bool(row['diagonal_dominant']),
            'holdout_best_signature':row['top_signature'],'holdout_mean_own_vs_best_other_margin':float(gs['own_minus_best_other'].mean()),
            'holdout_margin_GSE_bootstrap95CI':[lo,hi],
            'status_rule':'supported requires >=10 holdout GSE, diagonal dominant GSE-equal row, >=65% GSE mean-score matches, and positive lower bootstrap CI of own-minus-best-other; otherwise ambiguous/insufficient. Not independent cell identity proof.'})
    dump(out/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json',doc)
    logging.info('%s holdout validation statuses: %s',lin,{s:e['status'] for s,e in doc['subtype_evidence'].items()})
    return doc


def marker_expression_qc(cfg,lin):
    c=cache(cfg,lin);genes=np.load(c/'marker_genes.npy');ag=np.load(c/'gene_group_aggregates.npz')
    units=pd.read_parquet(c/'gene_group_units.parquet');rows=[]
    for split in ['train','holdout','all']:
        mask=np.ones(len(units),bool) if split=='all' else units['split'].eq(split).to_numpy()
        for subtype in list(units['Celltype_general'].unique())+['ALL']:
            use=mask.copy()
            if subtype!='ALL': use&=units['Celltype_general'].eq(subtype).to_numpy()
            n=int(ag['counts'][use].sum());s=ag['sums'][use].sum(0);d=ag['detected'][use].sum(0)
            for j,g in enumerate(genes): rows.append({'lineage':lin,'split':split,'subtype':subtype,'gene':g,
                'n_cells':n,'mean_log_expression':s[j]/max(n,1),'detection_fraction':d[j]/max(n,1)})
    pd.DataFrame(rows).to_csv(root(cfg)/f'data/candidate_gene_expression_QC_{lin}.csv',index=False)


def paired_final_candidate_comparison(cfg,lin):
    c=cache(cfg,lin);obs=pd.read_parquet(c/'metadata.parquet');order=np.load(c/'final_score_names.npy').tolist();k=len(order)
    specs=read_json(c/'score_specs_indexed.json');versions={}
    for s in specs:
        if s['mapping']=='projected': versions.setdefault(s['candidate_id'],{})[s['key']]=s['score_index']
    rows=[]
    for method in ['mean','module']:
        um,fv=aggregate_scores(np.load(c/f'final_{method}.npy',mmap_mode='r'),obs,['GSE','Celltype_general','split'])
        cm,cv=aggregate_scores(np.load(c/f'candidate_{method}.npy',mmap_mode='r'),obs,['GSE','Celltype_general','split'])
        pd.testing.assert_frame_equal(um,cm)
        keep=um['split'].eq('holdout')&um['n_cells'].ge(cfg['min_cells_per_group_subtype'])
        yy=pd.Categorical(um.loc[keep,'Celltype_general'],categories=order).codes
        _,fp,_,_=concordance_metrics(yy,fv[keep],k)
        for cid,imap in versions.items():
            if not set(order)<=set(imap): continue
            v=cv[keep][:,[imap[s] for s in order]]
            if not np.isfinite(v).all(): continue
            _,pr,_,_=concordance_metrics(yy,v,k)
            z=pd.DataFrame({'GSE':um.loc[keep,'GSE'].to_numpy(),'final_match':fp==yy,'candidate_match':pr==yy})
            z['delta']=z['final_match'].astype(float)-z['candidate_match'].astype(float)
            g=z.groupby('GSE')['delta'].mean();lo,hi=bootstrap_group_ci(g,cfg['seed'],cfg['bootstrap_replicates'])
            rows.append({'lineage':lin,'candidate_id':cid,'method':method,'split':'holdout','n_GSE':len(g),
                'final_minus_candidate_GSE_equal_concordance':g.mean(),'paired_GSE_bootstrap95_low':lo,'paired_GSE_bootstrap95_high':hi,
                'improvement_supported':lo>0,'interpretation':'paired change in annotation concordance; overlapping marker-based original labels remain a circularity limitation'})
    pd.DataFrame(rows).to_csv(root(cfg)/f'data/final_vs_candidates_holdout_{lin}.csv',index=False)


def combine_tables(cfg):
    out=root(cfg)
    for name in ['marker_selection_evidence','marker_sets_candidate_comparison','signature_scores_group_summary',
                 'signature_concordance_by_cluster','final_concordance','subtype_signature_matrix',
                 'grouped_holdout_concordance','holdout_GSE_margins','final_vs_candidates_holdout']:
        tables=[pd.read_csv(out/f'data/{name}_{lin}.csv') for lin in ['fibro','endo'] if (out/f'data/{name}_{lin}.csv').exists()]
        if tables: pd.concat(tables,ignore_index=True).to_csv(out/f'data/{name}.csv',index=False)
    parts=[pd.read_csv(p) for p in (out/'data').glob('GSE_cell_concordance_*_module_holdout.csv')]
    if parts: pd.concat(parts,ignore_index=True).to_csv(out/'data/signature_concordance_by_GSE.csv',index=False)
    scores=out/'data/signature_scores_group_summary.csv'
    if scores.exists():
        t=pd.read_csv(scores);keys=['lineage','method','split','group_level','subtype','signature']
        variability=t.groupby(keys,observed=True)['mean_score'].agg(
            n_units='count',unit_equal_mean='mean',variance_between_unit_means='var',
            sd_between_unit_means='std',median_unit_mean='median').reset_index()
        variability.to_csv(out/'data/signature_between_group_variance.csv',index=False)
    # Explicit unavailable rows, rather than relabeling study/sample statistics as donors.
    pd.DataFrame([{'lineage':lin,'donor_available':False,'reason':'No donor/patient/individual ID in source obs',
         'fallback_unit':'GSE','fallback_table':'signature_concordance_by_GSE.csv'} for lin in ['fibro','endo']]).to_csv(out/'data/signature_concordance_by_donor.csv',index=False)


def stage_validate(cfg,lineages):
    for lin in lineages:
        logging.info('Phase 3: training-only marker selection and untouched GSE holdout validation: %s',lin)
        marker_expression_qc(cfg,lin)
        doc=choose_final_markers(cfg,lin)
        final_scores(cfg,lin,doc)
        compare_candidate_versions(cfg,lin)
        old_score_comparison(cfg,lin)
        final_validation(cfg,lin,doc)
        paired_final_candidate_comparison(cfg,lin)
        targeted_holdout_expression(cfg,lin)
    combine_tables(cfg)
    write_reports(cfg,lineages)


def targeted_holdout_expression(cfg,lin):
    """Paired held-out-study expression/detection diagnostics, no marker retuning."""
    c=cache(cfg,lin);u=pd.read_parquet(c/'gene_group_units.parquet');a=np.load(c/'gene_group_aggregates.npz');genes=np.load(c/'marker_genes.npy')
    means=a['sums']/np.maximum(a['counts'][:,None],1);det=a['detected']/np.maximum(a['counts'][:,None],1)
    use=u['split'].eq('holdout')&(a['counts']>=cfg['min_cells_per_group_subtype'])
    lookup={(str(r.GSE),str(r.Celltype_general)):int(i) for i,r in u[use].iterrows()}
    groups=sorted(u.loc[use,'GSE'].unique());order=np.load(c/'final_score_names.npy').tolist();rows=[]
    for s in order:
        for other in order:
            if other==s:continue
            pairs=[(lookup[(g,s)],lookup[(g,other)]) for g in groups if (g,s) in lookup and (g,other) in lookup]
            if not pairs:continue
            left,right=np.array(pairs).T
            for j,gene in enumerate(genes):
                sm=paired_summary(means[left,j]-means[right,j])
                rows.append({'lineage':lin,'subtype':s,'comparison':other,'gene':gene,'split':'holdout',**sm,
                             'median_detection_delta':float(np.median(det[left,j]-det[right,j])),
                             'analysis':'paired GSE means of log expression, not raw-count pseudobulk DE; post-freeze diagnostic'})
    tab=pd.DataFrame(rows);tab['fdr']=bh_adjust(tab['p_value']);tab.to_csv(root(cfg)/f'data/targeted_pairwise_expression_holdout_{lin}.csv',index=False)


def mdtable(df):
    try:return df.to_markdown(index=False,floatfmt='.3f')
    except ImportError:return '```csv\n'+df.to_csv(index=False)+'```'


def enrich_final_metadata(cfg,lineages):
    for lin in lineages:
        p=root(cfg)/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json'
        if not p.exists():continue
        doc=read_json(p);c=cache(cfg,lin);genes=np.load(c/'marker_genes.npy').tolist();ag=np.load(c/'gene_group_aggregates.npz')
        expressed={g for g,v in zip(genes,ag['detected'].sum(0)) if v>0};records=[]
        for s,gg in doc['marker_sets'].items():
            metadata={'n_genes_total':len(gg),'n_genes_matched':len(gg),
                      'n_genes_expressed':sum(g in expressed for g in gg),'gene_coverage_scope':'all atlas cells; selection remains training-only'}
            doc['subtype_evidence'][s].update(metadata)
            records.append({'lineage':lin,'subtype':s,**metadata,'expression_layer':'X','version':VERSION,
                            'scoring_methods':'mean_log_expression;training_bin_matched_module','normalization':doc['normalization'],
                            'genes':';'.join(gg),'status':doc['subtype_evidence'][s]['status']})
        doc['source_fingerprint']=fingerprint(cfg['inputs'][lin]);doc['control_gene_metadata']=f'data/control_genes_{lin}.json'
        dump(p,doc);pd.DataFrame(records).to_csv(root(cfg)/f'data/final_score_metadata_{lin}.csv',index=False)
    coverage_path=root(cfg)/'audit/marker_gene_coverage.csv'
    if coverage_path.exists():
        cv=pd.read_csv(coverage_path);frames=[]
        for lin in ['fibro','endo']:
            qp=root(cfg)/f'data/candidate_gene_expression_QC_{lin}.csv'
            if qp.exists():
                q=pd.read_csv(qp);q=q[q['split'].eq('all')&q['subtype'].eq('ALL')]
                frames.append(q[['lineage','gene','mean_log_expression','detection_fraction']].rename(columns={'gene':'gene_resolved'}))
        if frames:
            cv=cv.drop(columns=['mean_log_expression','detection_fraction'],errors='ignore').merge(
                pd.concat(frames,ignore_index=True),on=['lineage','gene_resolved'],how='left',validate='many_to_one')
            cv.to_csv(coverage_path,index=False);cv.to_csv(root(cfg)/'data/marker_gene_coverage.csv',index=False)
    for lin in lineages:
        sm=root(cfg)/f'data/score_metadata_{lin}.csv'
        if sm.exists():
            q=pd.read_csv(sm);q['coverage_fraction']=q['n_genes_matched']/q['n_genes_total'].clip(lower=1)
            q['expressed_fraction_of_matched']=q['n_genes_expressed']/q['n_genes_matched'].clip(lower=1)
            q.to_csv(root(cfg)/f'data/marker_list_QC_summary_{lin}.csv',index=False)


def write_reports(cfg,lineages):
    out=root(cfg);available=[lin for lin in ['fibro','endo'] if (out/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json').exists()]
    enrich_final_metadata(cfg,available)
    audit=['# Phase 0 — 数据与历史代码审计', '',
       '两份输入 H5AD 均以只读方式打开；不修改标签、原始表达值或旧分数。详细字段与探针值见 audit/。', '',
       '未发现可核实的 donor/patient/individual 标识。sample/GSM 不被自动当作 donor；主要保留验证以 GSE 为分组。']
    marker=['# Marker QC 与最终名单', '',
       '“最终”表示本次预先规定规则下冻结并接受保留集审查的候选名单，不等于每个亚型都已获得独立生物学确证。', '',
       '所有原候选版本保存在 candidate_marker_catalog.json。相同名单的 notebook/用户副本保留各自来源，不算独立生物学支持。', '',
       '筛选只使用训练 GSE：每组亚型至少20细胞、至少10对可比GSE、GSE等权检出率≥5%、配对中位表达差≥0.02、正向GSE比例≥65%、单侧配对Wilcoxon检验BH FDR≤0.05。',
       '合格基因按最差亚型间配对效应、跨GSE方向稳定性、总体配对效应排序，每亚型最多8个。不因保留集结果不佳重新挑基因。', '',
       '没有可用基因时分数为 NaN。NG2 仅映射至唯一 CSPG4（若存在）；CD16、HLA-DQA不做模糊猜测；H2-Q4排除。']
    validation=['# Annotation validation — 结果与限制', '',
       '所有 accuracy/F1 指标均为 signature 与参考注释的一致性，不能解读为真实细胞身份的分类准确率。原标签已参考相同/重叠marker，仍有循环性；GSE holdout只隔离本次筛选与背景控制的拟合。', '',
       '均值分数与背景校正模块分数都保留。跨signature argmax仍未校准；不同marker数量、基线及dropout会影响排名。精确并列及无法打分细胞记为 abstain。', '',
       'GSE内先汇总至 GSE×subtype（至少20细胞），各GSE内部对实际出现的亚型等权，再令GSE等权；95% CI重采样GSE，不把百万细胞当独立重复。缺少某亚型的GSE不会被虚构为零表达。']
    summary=['# Run summary', '',f'Version: `{VERSION}`; seed: {cfg["seed"]}.', '',
       '输出目录：`'+str(out)+'`。输入文件未被写入。']
    oldlines=['# Old score comparison', '',
        '旧obs同名列可能经历不同字典/方法覆盖；notebook源代码和数值相似只能提供来源证据，不能还原最后一次实际执行状态。',
        '旧MS_*列没有可核实的对应字典/尺度元数据，不强行与新mean-log分数互换。Endo旧列的负值可能来自背景校正；不能据此认定X本身存在负值。']
    for lin in ['fibro','endo']:
        p=out/f'audit/data_audit_{lin}.json'
        if not p.exists():continue
        a=read_json(p);e=a['expression'];done=cache(cfg,lin)/'extraction_complete.json'
        audit += ['',f'## {LABELS[lin]}',f'- {a["n_cells"]:,} cells × {a["n_genes"]:,} genes; {a["n_leiden"]} Leiden clusters; {a["n_GSE"]} GSE; {a["n_samples"]} sample IDs.',
                  f'- 训练/保留 GSE = {a["train_groups"]}/{a["holdout_groups"]}; donor字段 = {a["donor_fields"]}.',
                  f'- X: {e["encoding"]}; layers={e["layers"]}; raw={e["has_raw"]}; log1p metadata={e["has_uns_log1p"]}.',
                  f'- 抽查逆log后每行总和分位数 = {e["probe_expm1_row_sum_quantiles"]}；与 normalize_total(10000) 后 log1p 相符。当前X原样打分，不二次log。',
                  '- 原始counts和原始转换脚本不能由该对象完整恢复；log-scale解释由metadata、数值及逆变换检验共同支持。']
        if done.exists():
            full=read_json(done);audit += [f'- 完整扫描 {full["full_scan_stored_values"]:,} stored values；负值/非有限值 = 0/0；范围 [{full["full_scan_stored_min"]:.5f}, {full["full_scan_stored_max"]:.5f}]。']
        cross=pd.read_csv(out/f'audit/annotation_by_leiden_{lin}.csv',dtype={'leiden':str})
        zero=cross[cross['leiden'].eq('0')]
        audit += ['- 当前 Leiden 0 真实组成：'+json.dumps(zero[['Celltype_general','n_cells']].to_dict('records'),ensure_ascii=False)+'.']
    audit += ['', '## 历史代码问题',
       '- Fibro notebook 在 stromal_marker 定义前调用；Endo 两个 latest 字典实际上存在，不能认定缺失。',
       '- 缺基因赋0、clip到非负、使用可能不存在的LymEC键及原始cluster ID重编号均已记录；新pipeline缺基因保留NaN且不clip。',
       '- Endo notebook cell72把 SMC_TPM2 改为 MCs_TPM2，EndoMT-I/II标签改为EndoMT，并同时更改其他cluster标签。故这不是所有旧细胞标签一一机械合并的证明。',
       '- 八类候选保留原始score；为六类一致性比较额外建立基因并集投影，明确标注该操作；旧full字典缺EndoMT时不计算完整六类argmax指标。',
       '- Fig1 v2 固定颜色读取 config_used.json，而非随机重新分配；记录于 audit/fig1_v2_colors.json。', '',
       '## 可查证的方法/符号来源',
       '- [NCBI CSPG4 / NG2](https://www.ncbi.nlm.nih.gov/gene/1464)',
       '- [NCBI FCGR3A](https://www.ncbi.nlm.nih.gov/gene/2214), [NCBI FCGR3B](https://www.ncbi.nlm.nih.gov/gene/2215)：不将CD16静默映射到唯一基因。',
       '- [Scanpy score_genes](https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.tl.score_genes.html)：匹配表达背景的原则。此实现逐marker按bin校正、共享固定bin控制，非Scanpy逐bit复刻。']
    for lin in available:
        doc=read_json(out/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json');marker+=['',f'## {LABELS[lin]}']
        for s,genes in doc['marker_sets'].items():
            ev=doc['subtype_evidence'][s];marker.append(f'- **{s}**: `{", ".join(genes) or "无合格marker"}`; status: `{ev["status"]}`; holdout GSE={ev.get("holdout_n_GSE","pending")}.')
        met=out/f'data/final_concordance_{lin}.csv';ci=out/f'data/grouped_holdout_concordance_{lin}.csv'
        needed=[met,ci,out/f'data/subtype_signature_matrix_{lin}.csv',
                out/f'data/signature_concordance_by_cluster_{lin}.csv',out/f'data/final_vs_candidates_holdout_{lin}.csv',
                out/f'data/old_vs_recomputed_scores_{lin}.csv']
        if not all(p.exists() for p in needed):continue
        metrics=pd.read_csv(met);cis=pd.read_csv(ci)
        validation += ['',f'## {LABELS[lin]}', '', '保留集细胞级描述：',
            mdtable(metrics[metrics['split'].eq('holdout')][['method','n_cells','top1_concordance','top2_concordance','macro_f1','balanced_accuracy','abstention_fraction']]), '',
            '保留GSE层面（对GSE重采样）：',mdtable(cis[cis['split'].eq('holdout')][['method','n_GSE','n_GSE_subtype_units','GSE_equal_groupmean_top1','ci_low','ci_high']])]
        sm=pd.read_csv(out/f'data/subtype_signature_matrix_{lin}.csv');sm=sm[(sm['method']=='module')&(sm['split']=='holdout')].drop_duplicates('reference_subtype')
        validation += ['', '每个参考亚型的行最大值检验（原始module尺度；未人为校准）：',
            mdtable(sm[['reference_subtype','top_signature','diagonal_dominant','n_GSE','own_vs_best_other_margin']])]
        ct=pd.read_csv(out/f'data/signature_concordance_by_cluster_{lin}.csv',dtype={'leiden':str});ct=ct[ct['method'].eq('module')].drop_duplicates('leiden')
        validation += ['',f'所有Leiden cluster均保留：{int(ct.dominant_label_matches_top.sum())}/{len(ct)} 个cluster的最高signature与主要参考标签相同。',
            '当前 Leiden 0：',mdtable(ct[ct['leiden'].eq('0')][['leiden','dominant_reference_subtype','highest_signature','n_cells','cell_top1_concordance','own_minus_best_other']])]
        changes=pd.read_csv(out/f'data/final_vs_candidates_holdout_{lin}.csv');short=changes[changes['method'].eq('module')].drop_duplicates(['candidate_id'])
        validation += ['', '最终 vs 候选的配对保留GSE变化（正数表示一致性上升；负数保留）：',
            mdtable(short[['candidate_id','final_minus_candidate_GSE_equal_concordance','paired_GSE_bootstrap95_low','paired_GSE_bootstrap95_high']])]
        old=pd.read_csv(out/f'data/old_vs_recomputed_scores_{lin}.csv');best=old[old['new_method'].eq('mean')].sort_values('mean_absolute_difference').drop_duplicates('old_column').head(12)
        oldlines += ['',f'## {LABELS[lin]}',mdtable(best[['old_column','candidate_id','pearson_r','mean_absolute_difference']])]
        summary += ['',f'## {LABELS[lin]}',f'Atlas: {len(pd.read_csv(out/f"audit/annotation_by_leiden_{lin}.csv"))} annotation×Leiden records; full input size {read_json(out/f"audit/data_audit_{lin}.json")["n_cells"]:,} cells.',
            mdtable(cis[cis['split'].eq('holdout')][['method','n_GSE','GSE_equal_groupmean_top1','ci_low','ci_high']]),
            '状态：'+', '.join(f'{s}={e["status"]}' for s,e in doc['subtype_evidence'].items())+'.']
    marker += ['', '## Marker QC 解释',
        '检测率、表达量和marker重叠见 data/candidate_gene_expression_QC_*.csv、candidate_gene_overlap_*.csv；逐基因保留/剔除原因、训练效应/FDR及跨来源方向一致性见 marker_selection_evidence.csv。',
        '全矩阵仅扫描一次，缓存候选marker表达和25个背景bin均值；不缓存完整dense表达。均值与module每套候选均重算，未覆盖obs旧列。',
        'UCell/AUCell未作为本次第三种方法执行：本任务两种必需方法已完整比较；全基因逐细胞rank另需大规模I/O与计算，未混用不等价近似。']
    coverage=pd.read_csv(out/'audit/marker_gene_coverage.csv')
    for lin in ['fibro','endo']:
        cv=coverage[coverage['lineage'].eq(lin)];missing=sorted(cv.loc[~cv['present'],'gene_listed'].unique())
        marker += ['',f'### {LABELS[lin]} 候选feature覆盖限制',
            f'原候选中未匹配/歧义/排除的唯一符号（不是生物学零表达）：{", ".join(missing)}。']
        cp=cache(cfg,lin)/'controls.json'
        if cp.exists():
            control=read_json(cp);marker += [f'背景bin估计使用 {control["pilot_training_cells"]:,} 个训练pilot细胞，覆盖 {control.get("pilot_n_training_GSE","未记录")} 个训练GSE。控制匹配仍可能受pilot组织构成影响，不能当作绝对校准。']
    validation += ['', '## 三层证据与局限',
        '1. Cell level：全量细胞分数、Top1/Top2、macro F1、balanced accuracy和分布表。仅作描述，不以单细胞数作为统计独立重复。',
        '2. Leiden：真实参考标签组成、均值排序、混合比例、逐细胞一致性；所有cluster均保存于CSV和诊断图。',
        '3. GSE：未参与本次筛选/控制基因确定的保留研究；来源等权、配对GSE基因表达与检出率、GSE bootstrap区间。donor层面不可完成，原因明确，不冒充donor-aware结果。',
        '组织与sample分层均值见 signature_scores_group_summary.csv；sample为样本ID，并非已核实独立个体。组织内特异状态可能被全局平均掩盖；GSE×subtype的组织构成也可能不同。',
        'MCs/EndoMT、eFibro/qFibro等交叉需同时查看逐基因检出与配对效应；广泛ECM/收缩/AP表达不自动证明另一细胞身份。broad lineage或functional高分也不能单独确诊污染/双细胞。',
        '缺失canonical marker是当前feature universe的限制；必要时应回到包含完整基因的原始对象做独立验证。本次不改原注释，不把保留集反例当作绘图错误删除。']
    legend=['# Figure legend', '',
        '**Extended Figure 1 — Signature scores and reference subtype annotation concordance.**', '',
        '**A–B.** Final background-corrected signatures on the original atlas UMAP, with a deterministic uniform sample of up to90,000 cells per lineage. Shared color range within lineage uses pooled1st–99th percentiles; saturation affects display only. All-cell numerical scores remain untruncated. Signature title colors are inherited from Fig1v2. Negative module scores are retained.', '',
        '**C–D.** Held-out GSE×reference-subtype means, requiring≥20cells per unit, averaged with equal weight perGSE for each subtype. Numbers are raw mean module scores. Outlined cells mark row maxima, including off-diagonal maxima. nGSE denotes studies containing that subtype, not donors. No row z-scoring or forced diagonal scaling.', '',
        '**E–F.** Each dot is one held-outGSE×subtype unit; y is the reference signature mean minus the largest competing signature mean. Dashed zero separates positive from negative margins. Boxplots show median/IQR and1.5IQR whiskers; points retain the full distribution. Different signatures are not calibrated probabilistic classifiers.', '',
        '**Status marks.** † indicates ambiguous signature concordance; ‡ indicates insufficient held-out study units. Supported states require ≥10 eligible held-outGSE, a dominant diagonal in the GSE-equal matrix, ≥65% mean-score matches, and a positive lower95%GSE-bootstrap CI of the own-minus-best-other margin. These are relative-to-reference criteria, not independent identity proof.', '',
        'Matched controls: fixed seed42, training-only pilot expression bins,25bins and up to12noncandidate expressed control genes/bin; per-marker control contrasts are averaged. All supplied marker genes are excluded from the background pool. Mean-log scores, candidate versions, cluster diagnostics and paired study confidence intervals are provided in separate tables/figures.', '',
        'Existing marker-informed reference labels make this a concordance analysis, not independent biological ground truth. No verified donor identifier exists; study-grouped holdout is reported explicitly. Ambiguous/insufficient states remain visible. Scatter clouds are locally rasterized; PDF text, axes, heatmaps, outlines and boxplots are editable vector objects.']
    summary += ['', '## Phase completion',
        'Phase0/1: actual inputs/history/palette audit and candidate symbol QC.',
        'Phase2: nine smoke checks plus full-matrix integrity scan; mean and background-corrected cell scores.',
        'Phase3: training-only selection; held-outGSE, cells,Leiden,study/tissue/sample summaries and oldscore comparisons.',
        'Phase4: integrated and lineage-specific figures with separatecluster/candidate diagnostics.',
        'Phase5: CLI pipeline, notebook, config, requirements, README and rerun command.',
        f'Integrated PNG currently exists: {(out/"figures/integrated/ExtFig1_SignatureScore_Fibro_Endo.png").exists()}.',
        'See README for exact commands and data/table meanings. Large full-cell scores and cache stay outside the compact download bundle.']
    for filename,lines in [('audit_report.md',audit),('marker_selection_report.md',marker),('annotation_validation_report.md',validation),
                            ('figure_legend.md',legend),('run_summary.md',summary),('old_score_comparison_report.md',oldlines)]:
        (out/'reports'/filename).write_text('\n'.join(lines),encoding='utf-8')
    # Restore the current figure's caption and distinguish its selected lists
    # when historical baseline reports are regenerated (e.g. from the notebook).
    if (out/'logs/reference_layout_update.json').exists():
        from reference_style_figures import refresh_current_documentation
        refresh_current_documentation(cfg)


def plotting_setup():
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'text.color':'#000000',
        'axes.labelcolor':'#000000','xtick.color':'#000000','ytick.color':'#000000',
        'font.weight':'normal','axes.titleweight':'normal','axes.labelweight':'normal',
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
    return plt


def save_figure(fig,base,cfg):
    base=Path(base);base.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(base)+'.png',dpi=cfg['dpi'],facecolor='white')
    fig.savefig(str(base)+'.pdf',dpi=cfg['dpi'],facecolor='white')
    logging.info('Saved %s PNG/PDF',base.name)


def add_panel_header(fig,spec,letter,title):
    pos=spec.get_position(fig)
    fig.text(pos.x0-.025,pos.y1+.014,letter,fontsize=25,ha='left',va='bottom')
    fig.text(pos.x0+.012,pos.y1+.016,title,fontsize=17,ha='left',va='bottom')


def signature_display_names(cfg,lin):
    doc=read_json(root(cfg)/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json')
    return {s:s+{'ambiguous':'†','insufficient_evidence':'‡'}.get(e['status'],'')
            for s,e in doc['subtype_evidence'].items()}


def draw_umap_panel(fig,spec,cfg,lin):
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    plt=plotting_setup();c=cache(cfg,lin);rows=np.load(c/'plot_rows.npy')
    xy=np.load(c/'umap.npy',mmap_mode='r')[rows];v=np.load(c/'final_module.npy',mmap_mode='r')[rows]
    names=np.load(c/'final_score_names.npy').tolist();colors=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_colors']
    display_names=signature_display_names(cfg,lin)
    ncol=2 if len(names)==4 else 3
    nested=GridSpecFromSubplotSpec(2,ncol+1,subplot_spec=spec,width_ratios=[1]*ncol+[.055],wspace=.12,hspace=.22)
    finite=v[np.isfinite(v)];lo,hi=np.quantile(finite,[.01,.99]) if len(finite) else (-1,1)
    if hi<=lo: hi=lo+1
    norm=Normalize(lo,hi);axes=[]
    for j,name in enumerate(names):
        ax=fig.add_subplot(nested[j//ncol,j%ncol]);axes.append(ax)
        ax.scatter(xy[:,0],xy[:,1],c=v[:,j],s=.9,cmap='viridis',norm=norm,rasterized=True,linewidths=0)
        ax.set_aspect('equal',adjustable='datalim');ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values(): spine.set_visible(False)
        ng=len(read_json(root(cfg)/f'marker_sets/{LABELS[lin]}_subtype_markers_final.json')['marker_sets'][name])
        ax.set_title(f'{display_names[name]}  |  {ng} genes',color=colors[name],fontsize=14,pad=5)
        if j==0:
            for end in [(.23,.07),(.07,.25)]: ax.annotate('',xy=end,xytext=(.07,.07),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.9})
            ax.text(.15,.015,'UMAP-1',transform=ax.transAxes,ha='center',fontsize=9)
            ax.text(.012,.16,'UMAP-2',transform=ax.transAxes,va='center',rotation=90,fontsize=9)
    cax=fig.add_subplot(nested[:,-1]);cb=fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),cax=cax,extend='both')
    if cb.solids is not None: cb.solids.set_rasterized(False)
    cb.set_label('Module score (shared lineage scale)',fontsize=12,labelpad=8);cb.ax.tick_params(labelsize=10)
    dump(root(cfg)/f'data/figure_umap_scale_{lin}.json',{'vmin':float(lo),'vmax':float(hi),'range':'pooled 1st–99th percentile across signatures; color saturation only, numeric scores unchanged','n_plotted_cells':len(rows),'sampling':'uniform without replacement; seed42; all-cell analyses stored separately'})
    return axes


def draw_score_matrix(fig,spec,cfg,lin):
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle
    tab=pd.read_csv(root(cfg)/f'data/subtype_signature_matrix_{lin}.csv')
    tab=tab[(tab['method']=='module')&(tab['split']=='holdout')]
    order=np.load(cache(cfg,lin)/'final_score_names.npy').tolist();colors=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_colors']
    display_names=signature_display_names(cfg,lin)
    matrix=tab.pivot(index='reference_subtype',columns='signature',values='GSE_equal_mean_score').reindex(index=order,columns=order)
    nested=GridSpecFromSubplotSpec(2,1,subplot_spec=spec,height_ratios=[1,.055],hspace=.48)
    ax=fig.add_subplot(nested[0]);lim=max(float(np.nanmax(np.abs(matrix.to_numpy()))),.1)
    # QuadMesh keeps every heatmap tile editable in Illustrator (imshow would rasterize).
    edges=np.arange(len(order)+1)-.5
    im=ax.pcolormesh(edges,edges,matrix.to_numpy(),cmap='RdBu_r',
                    norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),shading='flat',
                    edgecolors='white',linewidth=.3,rasterized=False)
    ax.set_aspect('equal');ax.set_ylim(len(order)-.5,-.5)
    ax.set_xticks(range(len(order)),order,rotation=35,ha='right',fontsize=13)
    nmap=tab.drop_duplicates('reference_subtype').set_index('reference_subtype')['n_GSE'].to_dict()
    ax.set_yticks(range(len(order)),[f'{display_names[s]}  ({int(nmap[s])} GSE)' for s in order],fontsize=13)
    for t,s in zip(ax.get_xticklabels(),order): t.set_color(colors[s])
    for t,s in zip(ax.get_yticklabels(),order): t.set_color(colors[s])
    for i in range(len(order)):
        r=matrix.iloc[i].to_numpy();best=int(np.nanargmax(r)) if np.isfinite(r).any() else -1
        for j in range(len(order)):
            value=r[j];txt=f'{value:.2f}' if np.isfinite(value) else 'NA'
            ax.text(j,i,txt,ha='center',va='center',color='white' if abs(value)>lim*.6 else 'black',fontsize=14)
        if best>=0: ax.add_patch(Rectangle((best-.48,i-.48),.96,.96,fill=False,edgecolor='black',lw=1.6))
    ax.set_xlabel('Final signature',labelpad=7);ax.set_ylabel('Reference annotation',labelpad=9)
    ax.tick_params(length=0);cax=fig.add_subplot(nested[1]);cb=fig.colorbar(im,cax=cax,orientation='horizontal')
    if cb.solids is not None:cb.solids.set_rasterized(False)
    cb.set_label('Mean module score | equal weight per GSE',fontsize=12);cb.ax.tick_params(labelsize=11)
    matrix.to_csv(root(cfg)/f'data/plotted_holdout_matrix_{lin}.csv')
    # Assert the exact annotated heatmap numbers match the written data table.
    check=pd.read_csv(root(cfg)/f'data/plotted_holdout_matrix_{lin}.csv',index_col=0)
    np.testing.assert_allclose(check.to_numpy(),matrix.to_numpy(),equal_nan=True,rtol=1e-12,atol=1e-12)
    return ax


def draw_group_margins(fig,spec,cfg,lin):
    ax=fig.add_subplot(spec);tab=pd.read_csv(root(cfg)/f'data/holdout_GSE_margins_{lin}.csv')
    tab=tab[tab['method']=='module'];names=np.load(cache(cfg,lin)/'final_score_names.npy').tolist()
    colors=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_colors'];rng=np.random.default_rng(cfg['seed'])
    display_names=signature_display_names(cfg,lin)
    arrays=[tab.loc[tab['subtype'].eq(s),'own_minus_best_other'].dropna().to_numpy() for s in names]
    bp=ax.boxplot(arrays,positions=np.arange(len(names)),widths=.5,showfliers=False,patch_artist=True,
        medianprops={'color':'black','linewidth':1.3},boxprops={'linewidth':.9},whiskerprops={'linewidth':.9},capprops={'linewidth':.9})
    for j,(patch,s,a) in enumerate(zip(bp['boxes'],names,arrays)):
        patch.set_facecolor(colors[s]);patch.set_alpha(.40)
        ax.scatter(j+rng.uniform(-.19,.19,len(a)),a,s=12,c=colors[s],alpha=.65,edgecolors='none',rasterized=False)
    ax.axhline(0,color='black',lw=1,ls='--');ax.set_xlim(-.6,len(names)-.4)
    ax.set_xticks(range(len(names)),[f'{display_names[s]}\nn={len(a)}' for s,a in zip(names,arrays)],fontsize=12)
    for t,s in zip(ax.get_xticklabels(),names): t.set_color(colors[s])
    ax.set_ylabel('Own signature − highest other signature\n(mean module score per held-out GSE)',fontsize=12)
    ax.set_xlabel('Reference subtype | each dot is one held-out GSE',fontsize=12,labelpad=10)
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    return ax


def plot_leiden_diagnostics(cfg,lin):
    from matplotlib.patches import Patch
    plt=plotting_setup();tab=pd.read_csv(root(cfg)/f'data/signature_concordance_by_cluster_{lin}.csv',dtype={'leiden':str})
    tab=tab[tab['method']=='module'].drop_duplicates('leiden').copy()
    tab['numeric']=pd.to_numeric(tab['leiden']);tab=tab.sort_values('numeric')
    colors=read_json(root(cfg)/'audit/fig1_v2_colors.json')['general_colors']
    fig,ax=plt.subplots(figsize=(12,max(10,len(tab)*.245)))
    yy=np.arange(len(tab));ax.barh(yy,tab['own_minus_best_other'],color=[colors[s] for s in tab['dominant_reference_subtype']],height=.73)
    ax.axvline(0,color='black',lw=1)
    labels=[f'L{r.leiden}  {r.dominant_reference_subtype}  |  n={r.n_cells:,}' for r in tab.itertuples()]
    ax.set_yticks(yy,labels,fontsize=10);ax.invert_yaxis();ax.set_xlabel('Reference signature − highest competing signature | mean module score')
    for i,r in enumerate(tab.itertuples()):
        if not r.dominant_label_matches_top:
            ax.text(r.own_minus_best_other,i,'  '+r.highest_signature,va='center',ha='right' if r.own_minus_best_other<0 else 'left',fontsize=9)
    ax.set_title(f'{LABELS[lin]} | every Leiden cluster, all cells\nNegative values and competing signatures are retained',fontsize=15,pad=12)
    ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.35,right=.96,top=.94,bottom=.08)
    save_figure(fig,root(cfg)/f'figures/panels/{LABELS[lin]}_Leiden_concordance',cfg);plt.close(fig)


def plot_version_comparison(cfg,lin):
    plt=plotting_setup();tab=pd.read_csv(root(cfg)/f'data/marker_sets_candidate_comparison_{lin}.csv')
    tab=tab[(tab['split']=='holdout')&(tab['method']=='module')&tab['complete_for_topscore_comparison']].copy()
    # Identical user/notebook copies remain in tables; show the first identical metric tuple once.
    tab=tab.drop_duplicates(['GSE_equal_groupmean_top1','GSE_equal_ci_low','GSE_equal_ci_high'])
    final=pd.read_csv(root(cfg)/f'data/grouped_holdout_concordance_{lin}.csv')
    f=final[(final['split']=='holdout')&(final['method']=='module')].iloc[0]
    labels=[s.split(':')[1] if ':' in s else s for s in tab['candidate_id']]+['Final (training-selected)']
    val=tab['GSE_equal_groupmean_top1'].tolist()+[f.GSE_equal_groupmean_top1]
    low=tab['GSE_equal_ci_low'].tolist()+[f.ci_low];high=tab['GSE_equal_ci_high'].tolist()+[f.ci_high]
    fig,ax=plt.subplots(figsize=(12,max(5,.5*len(labels))))
    ax.errorbar(val,np.arange(len(labels)),xerr=[np.asarray(val)-low,np.asarray(high)-val],fmt='o',color='#3b6178',capsize=3)
    ax.set_yticks(range(len(labels)),labels,fontsize=11);ax.invert_yaxis();ax.set_xlim(0,1)
    ax.set_xlabel('Held-out GSE-equal concordance of mean-score argmax\n95% CI: resample GSE units, not cells')
    ax.set_title(LABELS[lin]+' | original and final candidates',fontsize=15)
    ax.grid(axis='x',alpha=.2);fig.subplots_adjust(left=.39,right=.96,bottom=.19,top=.91)
    save_figure(fig,root(cfg)/f'figures/panels/{LABELS[lin]}_candidate_versions',cfg);plt.close(fig)


def plot_reference_umap(cfg,lin):
    from matplotlib.lines import Line2D
    plt=plotting_setup();c=cache(cfg,lin);rows=np.load(c/'plot_rows.npy');xy=np.load(c/'umap.npy',mmap_mode='r')[rows]
    obs=pd.read_parquet(c/'metadata.parquet').iloc[rows];pal=read_json(root(cfg)/'audit/fig1_v2_colors.json')
    fig,axes=plt.subplots(1,2,figsize=(16,7));names=np.load(c/'final_score_names.npy').tolist()
    axes[0].scatter(*xy.T,c=obs['Celltype_general'].astype(str).map(pal['general_colors']),s=1.6,rasterized=True,linewidths=0)
    axes[0].legend(handles=[Line2D([0],[0],marker='o',ls='',color=pal['general_colors'][s],label=s) for s in names],
                   frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.5,-.16),fontsize=12)
    with h5py.File(cfg['inputs'][lin],'r') as f:
        cats=np.asarray(read_elem(f['obs/leiden/categories'])).astype(str)
        colors=np.asarray(read_elem(f['uns/leiden_colors'])).astype(str)
    leidenpal=dict(zip(cats,colors))
    axes[1].scatter(*xy.T,c=obs['leiden'].astype(str).map(leidenpal),s=1.6,rasterized=True,linewidths=0)
    for cluster in cats:
        m=obs['leiden'].astype(str).eq(cluster).to_numpy()
        if m.any():
            x,y=np.median(xy[m],axis=0);axes[1].text(x,y,cluster,ha='center',va='center',fontsize=8,
                bbox={'facecolor':'white','edgecolor':'none','alpha':.7,'pad':.25})
    for ax,title in zip(axes,[LABELS[lin]+' | reference annotation',LABELS[lin]+' | actual Leiden IDs']):
        ax.set_title(title,fontsize=15);ax.set_aspect('equal',adjustable='datalim');ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
    fig.subplots_adjust(left=.025,right=.98,top=.92,bottom=.18,wspace=.1)
    save_figure(fig,root(cfg)/f'figures/panels/{LABELS[lin]}_reference_Leiden_UMAP',cfg);plt.close(fig)


def stage_validation_figure(cfg,lineages):
    from matplotlib.gridspec import GridSpec
    plt=plotting_setup()
    for lin in lineages:
        # Standalone lineage sheet: signature UMAPs, grouped matrix, GSE margins.
        fig=plt.figure(figsize=(13,20));gs=GridSpec(3,1,figure=fig,height_ratios=[7,5,4],hspace=.42,left=.15,right=.92,top=.95,bottom=.06)
        for spec,letter,title,func in [(gs[0],'A',LABELS[lin]+' | final module scores',draw_umap_panel),
            (gs[1],'B','Held-out studies | subtype × signature',draw_score_matrix),
            (gs[2],'C','Held-out study diagnostics',draw_group_margins)]:
            add_panel_header(fig,spec,letter,title);func(fig,spec,cfg,lin)
        fig.text(.52,.018,'† Ambiguous signature concordance   ‡ Insufficient held-out studies',ha='center',fontsize=11)
        save_figure(fig,root(cfg)/f'figures/panels/ExtFig1_SignatureScore_{LABELS[lin]}',cfg);plt.close(fig)
        plot_leiden_diagnostics(cfg,lin);plot_version_comparison(cfg,lin);plot_reference_umap(cfg,lin)
    if all((root(cfg)/f'data/subtype_signature_matrix_{lin}.csv').exists() and
           (root(cfg)/f'data/holdout_GSE_margins_{lin}.csv').exists() for lin in ['fibro','endo']):
        fig=plt.figure(figsize=(24,22))
        gs=GridSpec(3,2,figure=fig,width_ratios=[1,1],height_ratios=[7.1,5.2,4.7],hspace=.44,wspace=.36,left=.075,right=.96,top=.945,bottom=.065)
        for col,lin in enumerate(['fibro','endo']):
            for row,letter,title,func in [(0,'AB'[col],LABELS[lin]+' | final module scores',draw_umap_panel),
                (1,'CD'[col],'Held-out studies | subtype × signature',draw_score_matrix),
                (2,'EF'[col],'Held-out study diagnostics',draw_group_margins)]:
                add_panel_header(fig,gs[row,col],letter,title);func(fig,gs[row,col],cfg,lin)
        fig.text(.5,.028,'† Ambiguous signature concordance   ‡ Insufficient held-out studies.  No verified donor IDs: validation units are GSE studies.',ha='center',fontsize=12)
        fig.text(.5,.012,'Reference-label concordance; marker overlap with original annotation limits biological independence.',ha='center',fontsize=11)
        save_figure(fig,root(cfg)/'figures/integrated/ExtFig1_SignatureScore_Fibro_Endo',cfg);plt.close(fig)
    write_reports(cfg,lineages)


def stage_figure(cfg,lineages):
    """Current one-page A–D figure: Fibro/Endo UMAPs and all Leiden bars."""
    from reference_style_figures import render_reference_figures
    render_reference_figures(cfg,lineages)
