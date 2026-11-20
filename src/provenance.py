#!/usr/bin/env python3
"""
provenance.py
-------------
Records what produced a result, so any output file can be traced back to the
exact inputs, code and software versions that made it.

The problem this solves: six months from now you open `out/annotation.json`
and cannot tell which version of the input files it came from, which version of
Biopython was installed, or whether the numbers in your write-up match the file
in front of you. Results without provenance are anecdotes.

Every script calls `stamp()` and puts the result in its JSON output under
"provenance". That block contains:

  timestamp     when the run happened (UTC, ISO 8601)
  git_commit    the exact code state, if the repo is clean enough to say
  git_dirty     whether there were uncommitted changes (a dirty run is not
                reproducible, and saying so is better than hiding it)
  inputs        every input file with its size and SHA-256 checksum
  packages      versions of the libraries whose numbers affect the result
  python        interpreter version

The checksums are the important part. If an input file changes by one byte,
its checksum changes completely, so you can always prove whether a result
corresponds to the data you think it does.
"""

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Libraries whose version can change a number in the output. Packages that only
# affect formatting are deliberately not listed.
TRACKED_PACKAGES = ["biopython", "ViennaRNA", "pandas", "numpy"]


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _git(args, root):
    try:
        out = subprocess.run(["git", "-C", str(root)] + args,
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def _package_versions():
    versions = {}
    try:
        from importlib.metadata import version, PackageNotFoundError
    except ImportError:
        return versions
    for name in TRACKED_PACKAGES:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            pass
    return versions


def stamp(script, inputs, root=None, notes=None):
    """Build a provenance record.

    script : the script producing the output, normally __file__
    inputs : list of Path objects the run depended on
    notes  : optional dict of anything else worth recording
    """
    script = Path(script).resolve()
    root = Path(root) if root else script.parents[1]

    recorded = []
    for p in inputs:
        p = Path(p)
        if p.exists():
            recorded.append(dict(
                path=str(p.relative_to(root)) if root in p.parents else str(p),
                bytes=p.stat().st_size,
                sha256=sha256(p),
            ))
        else:
            recorded.append(dict(path=str(p), missing=True))

    commit = _git(["rev-parse", "HEAD"], root)
    status = _git(["status", "--porcelain"], root)

    rec = dict(
        script=script.name,
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        git_commit=commit,
        git_dirty=bool(status) if status is not None else None,
        inputs=recorded,
        packages=_package_versions(),
        python=platform.python_version(),
        platform=platform.platform(),
    )
    if notes:
        rec["notes"] = notes
    return rec


def summary_line(rec):
    """One line for the top of a text report."""
    commit = (rec["git_commit"] or "no-git")[:8]
    dirty = " (uncommitted changes)" if rec.get("git_dirty") else ""
    return f"run {rec['timestamp_utc']}  commit {commit}{dirty}"


if __name__ == "__main__":
    # Running this file directly prints a stamp for the repo, which is a quick
    # way to check the environment is what you think it is.
    root = Path(__file__).resolve().parents[1]
    print(json.dumps(stamp(__file__, [], root=root), indent=2))
