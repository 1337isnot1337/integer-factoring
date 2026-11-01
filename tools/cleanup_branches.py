#!/usr/bin/env python3
"""Delete old remote branches fully contained in the default branch.

Preview by default; --apply enables deletion. Requires Git and authenticated gh.
Local branches and worktrees are never deleted.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def run(*command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def api(repository, endpoint, paginate=False):
    command = ["gh", "api", f"repos/{repository}{endpoint}"]
    if paginate:
        command += ["--paginate", "--slurp"]
    result = json.loads(run(*command))
    return [item for page in result for item in page] if paginate else result


def open_pr_branches(repository):
    pulls = api(repository, "/pulls?state=open&per_page=100", paginate=True)
    return {p[side]["ref"] for p in pulls for side in ("head", "base")}


def integrated(sha, default_sha):
    result = subprocess.run(
        ["git", "merge-base", "--is-ancestor", sha, default_sha], cwd=ROOT,
        check=False,
    )
    if result.returncode not in (0, 1):
        raise RuntimeError(f"Cannot inspect ancestry of {sha}")
    return result.returncode == 0


def cleanup(repository, apply=False):
    default = api(repository, "")["default_branch"]
    run("git", "fetch", "--prune", "origin")
    default_sha = run("git", "rev-parse", f"refs/remotes/origin/{default}")
    in_use = open_pr_branches(repository)
    cutoff = time.time() - 24 * 60 * 60
    candidates = []
    for branch in api(repository, "/branches?per_page=100", paginate=True):
        name, sha = branch["name"], branch["commit"]["sha"]
        if name == default or branch["protected"] or name in in_use:
            continue
        if not integrated(sha, default_sha):
            continue
        if int(run("git", "show", "-s", "--format=%ct", sha)) > cutoff:
            continue
        candidates.append((name, sha))
        print(f"Integrated for cleanup: {name} at {sha}")

    if not apply or not candidates:
        print(f"{'Preview' if not apply else 'Cleanup'}: {len(candidates)} eligible branches")
        return

    # Refresh the default tip and PR references immediately before deletion.
    if api(repository, "")["default_branch"] != default:
        raise RuntimeError("Default branch changed; retry cleanup")
    run("git", "fetch", "--prune", "origin")
    default_sha = run("git", "rev-parse", f"refs/remotes/origin/{default}")
    branches = {b["name"]: b for b in api(
        repository, "/branches?per_page=100", paginate=True,
    )}
    in_use = open_pr_branches(repository)
    for name, sha in candidates:
        current = branches.get(name)
        if (current is None or current["commit"]["sha"] != sha
                or current["protected"] or name in in_use
                or not integrated(sha, default_sha)):
            raise RuntimeError(f"Branch changed or is in use: {name}; retry cleanup")

    # Exact leases also reject pushes made after the API checks. Atomic push
    # ensures a changed branch aborts the whole batch without partial deletion.
    run("git", "push", "--atomic",
        *(f"--force-with-lease=refs/heads/{name}:{sha}" for name, sha in candidates),
        "origin", *(f":refs/heads/{name}" for name, sha in candidates))
    print(f"Deleted {len(candidates)} integrated remote branches")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    repository = os.environ.get("GITHUB_REPOSITORY") or run(
        "gh", "repo", "view", "--json", "nameWithOwner", "--jq", ".nameWithOwner",
    )
    cleanup(repository, apply=args.apply)


if __name__ == "__main__":
    main()
