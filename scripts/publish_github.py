#!/usr/bin/env python3
"""Create an empty GitHub repository and push the reviewed, committed source."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(command, env=None):
    return subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True)


def checked(command, env=None):
    p = run(command, env)
    if p.returncode:
        raise RuntimeError(p.stderr.strip() or p.stdout.strip() or 'Command failed')
    return p.stdout.strip()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', default='ytwang21/TabulaStromal')
    p.add_argument('--visibility', choices=['public', 'private'], default='public')
    p.add_argument('--gh', default=shutil.which('gh'))
    a = p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', a.repo):
        p.error('Expected owner/name for --repo.')
    if not a.gh or not Path(a.gh).is_file():
        p.error('Install gh and log in: gh auth login --scopes repo,workflow')
    gh = str(Path(a.gh).resolve())
    try:
        checked([sys.executable, 'scripts/validate_repository.py', '--verify-checksums'])
        if checked(['git', 'status', '--porcelain']):
            raise RuntimeError('Commit the reviewed source tree before publishing.')
        if checked(['git', 'branch', '--show-current']) != 'main':
            raise RuntimeError('Publish from main.')
        owner, name = a.repo.split('/')
        login = checked([gh, 'api', 'user', '--jq', '.login'])
        if login.casefold() != owner.casefold():
            raise RuntimeError(f'Authenticated owner {login} differs from {owner}.')
        lookup = run([gh, 'repo', 'view', a.repo, '--json', 'nameWithOwner,isEmpty,isPrivate'])
        if lookup.returncode:
            if 'Could not resolve to a Repository' not in lookup.stderr and '404' not in lookup.stderr:
                raise RuntimeError('Cannot verify remote existence: ' + lookup.stderr)
            checked([gh, 'api', '--method', 'POST', 'user/repos', '-f', f'name={name}',
                     '-f', 'description=Human stromal atlas, age prediction and CellAge analysis workflows',
                     '-F', f'private={str(a.visibility == "private").lower()}', '-F', 'auto_init=false'])
        else:
            info = json.loads(lookup.stdout)
            if not info['isEmpty']:
                raise RuntimeError('Remote is nonempty; refusing to overwrite its history.')
            if info['isPrivate'] != (a.visibility == 'private'):
                raise RuntimeError('Remote visibility differs from the requested visibility.')
        url = f'https://github.com/{a.repo}.git'
        origin = run(['git', 'remote', 'get-url', 'origin'])
        if origin.returncode:
            checked(['git', 'remote', 'add', 'origin', url])
        elif origin.stdout.strip() != url:
            raise RuntimeError('Origin differs; it was not changed.')
        # This helper applies only to this process, never global Git settings.
        env = {**os.environ, 'GIT_CONFIG_COUNT': '1',
               'GIT_CONFIG_KEY_0': 'credential.https://github.com.helper',
               'GIT_CONFIG_VALUE_0': '!' + gh + ' auth git-credential',
               'GIT_TERMINAL_PROMPT': '0'}
        checked(['git', 'push', '--set-upstream', 'origin', 'main'], env)
        local = checked(['git', 'rev-parse', 'HEAD'])
        remote = checked([gh, 'api', f'repos/{a.repo}/commits/main', '--jq', '.sha'])
        if local != remote:
            raise RuntimeError('Remote HEAD differs from the local reviewed commit.')
        print(f'Published https://github.com/{a.repo} at {remote}')
        return 0
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
