#!/usr/bin/env python3
"""Run all article figure scripts.

This script discovers every Python file in ``examples/article_figures`` and
executes them in sorted order from the repository root. Each figure script keeps
its own output behavior, so generated files are written wherever the individual
script sends them, usually the VLab4Mic output directory.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
FIGURE_SCRIPTS_DIR = REPO_ROOT / "examples" / "article_figures"


def discover_scripts(pattern: str = "*.py") -> list[Path]:
    """Return all article figure scripts in deterministic order."""
    return sorted(FIGURE_SCRIPTS_DIR.glob(pattern))


def run_script(script: Path, *, env: dict[str, str]) -> int:
    """Run one figure script and return its process exit code."""
    relative_script = script.relative_to(REPO_ROOT)
    print(f"\n=== Running {relative_script} ===", flush=True)
    completed = subprocess.run(
        [sys.executable, str(relative_script)],
        cwd=REPO_ROOT,
        env=env,
        check=False,
    )
    if completed.returncode == 0:
        print(f"=== Finished {relative_script} ===", flush=True)
    else:
        print(
            f"=== Failed {relative_script} with exit code {completed.returncode} ===",
            flush=True,
        )
    return completed.returncode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all scripts in examples/article_figures."
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue running remaining scripts after a script fails.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the scripts that would run without executing them.",
    )
    parser.add_argument(
        "--pattern",
        default="*.py",
        help="Glob pattern to select scripts inside examples/article_figures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    scripts = discover_scripts(args.pattern)

    if not FIGURE_SCRIPTS_DIR.is_dir():
        print(f"Figure script directory not found: {FIGURE_SCRIPTS_DIR}", file=sys.stderr)
        return 2
    if not scripts:
        print(f"No scripts found in {FIGURE_SCRIPTS_DIR} matching {args.pattern}")
        return 0

    print("Article figure scripts:")
    for script in scripts:
        print(f" - {script.relative_to(REPO_ROOT)}")

    if args.dry_run:
        return 0

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    failures: list[tuple[Path, int]] = []
    passed = 0
    ran = 0

    for script in scripts:
        ran += 1
        returncode = run_script(script, env=env)
        if returncode != 0:
            failures.append((script, returncode))
            if not args.keep_going:
                break
        else:
            passed += 1

    print("\n=== Summary ===")
    print(f"Scripts discovered: {len(scripts)}")
    print(f"Scripts run: {ran}")
    print(f"Scripts completed: {passed}")
    print(f"Scripts failed: {len(failures)}")

    if failures:
        for script, returncode in failures:
            print(f" - {script.relative_to(REPO_ROOT)}: exit code {returncode}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
