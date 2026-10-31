#!/usr/bin/env python3
"""Verify recorded source hashes against the worktree or read-only git history."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import subprocess


def resolve_source_hashes(expected, source_root=None):
    """Find exact bytes, without executing old code or replacing recorded hashes.

    A newer checker can replay the mathematics while this separate check
    identifies the historical implementation that generated the record.
    The lookup is local and read-only; it never fetches or changes git state.
    """
    source_root = Path(source_root or Path(__file__).resolve().parent)
    repo = source_root.parent
    resolved = {}
    unresolved = []
    for name, digest in expected.items():
        path = source_root/name
        if path.is_file() and sha256(path.read_bytes()).hexdigest() == digest:
            resolved[name] = {'location': 'worktree', 'sha256': digest}
            continue
        relative = str(path.relative_to(repo))
        history = subprocess.run(['git', 'log', '--all', '--format=%H', '--', relative],
                                 cwd=repo, capture_output=True, text=True)
        found = False
        if history.returncode == 0:
            for revision in history.stdout.splitlines():
                blob = subprocess.run(['git', 'show', revision+':'+relative], cwd=repo, capture_output=True)
                if blob.returncode == 0 and sha256(blob.stdout).hexdigest() == digest:
                    resolved[name] = {'location': 'git', 'revision': revision, 'path': relative, 'sha256': digest}
                    found = True
                    break
        if not found:
            unresolved.append(name)
    return resolved, unresolved
