#!/usr/bin/env python3
"""Launch the native Nextflow workflow using a JSON configuration.

The optional JSON 'environment' object selects executables and a MAJIQ license
for this run. It does not orchestrate sample processing.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def configure_environment(config, base_dir):
    settings = config.get("environment", {})
    if not isinstance(settings, dict):
        raise ValueError("JSON 'environment' must be an object")

    env = os.environ.copy()
    conda = settings.get("conda_env")
    if conda:
        prefix = Path(conda).expanduser()
        if not prefix.is_absolute():
            prefix = base_dir / prefix
        prefix = prefix.resolve()
        bin_dir = prefix / "bin"
        if not bin_dir.is_dir():
            raise FileNotFoundError(f"Conda environment bin directory not found: {bin_dir}")
        env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")

    license_name = settings.get("majiq_license_file")
    if license_name:
        license_file = Path(license_name).expanduser()
        if not license_file.is_absolute():
            license_file = base_dir / license_file
        license_file = license_file.resolve()
        if not license_file.is_file():
            raise FileNotFoundError(f"MAJIQ license file not found: {license_file}")
        env["MAJIQ_LICENSE_FILE"] = str(license_file)

    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params-file", required=True, help="Nextflow JSON parameters")
    parser.add_argument("-c", "--config", default=str(HERE / "splicescape-v2.config"))
    parser.add_argument("-profile", "--profile", default="local")
    parser.add_argument("--nextflow", default="nextflow")
    parser.add_argument("--resume", action="store_true")
    args, extra = parser.parse_known_args()

    params_path = Path(args.params_file).expanduser().resolve()
    config_path = Path(args.config).expanduser().resolve()
    if not params_path.is_file():
        parser.error(f"JSON file not found: {params_path}")
    if not config_path.is_file():
        parser.error(f"Nextflow config not found: {config_path}")

    try:
        settings = json.loads(params_path.read_text())
        if not isinstance(settings, dict):
            raise ValueError("Top-level JSON must be an object")
        env = configure_environment(settings, params_path.parent)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    if not shutil.which(args.nextflow, path=env.get("PATH", "")):
        parser.error(f"Nextflow executable not found: {args.nextflow}")
    if settings.get("environment", {}).get("conda_env") and not shutil.which("STAR", path=env["PATH"]):
        parser.error("STAR not found in configured PATH")

    db_value = settings.get("db")
    if not isinstance(db_value, str) or not db_value.strip():
        parser.error("JSON must specify a non-empty 'db' path for run-wide locking")
    db_path = Path(db_value).expanduser()
    if not db_path.is_absolute():
        db_path = params_path.parent / db_path
    db_path = db_path.resolve()
    if not db_path.parent.is_dir():
        parser.error(f"Database parent directory does not exist: {db_path.parent}")

    cmd = [
        args.nextflow, "run", str(HERE / "splicescape.nf"),
        "-c", str(config_path), "-profile", args.profile,
        "-params-file", str(params_path),
    ]
    if args.resume:
        cmd.append("-resume")
    cmd.extend(extra)
    # The whole preflight -> analysis -> import -> replacement lifecycle is
    # protected, not just the final SQLite write. The lock file must persist.
    # Advisory: launchers bypassing run_splicescape.py do not honor the lock.
    lock_path = db_path.with_name(db_path.name + ".splicescape.lock")
    try:
        with lock_path.open("a+") as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print(f"Database busy: another SpliceScape launcher holds {lock_path}",
                      file=sys.stderr, flush=True)
                return 73
            print(f"Database run lock acquired: {lock_path}", flush=True)
            print("Launching Nextflow with JSON-specified environment", flush=True)
            return subprocess.run(cmd, cwd=HERE, env=env, check=False).returncode
    except OSError as exc:
        print(f"Cannot acquire database run lock {lock_path}: {exc}",
              file=sys.stderr, flush=True)
        return 73


if __name__ == "__main__":
    sys.exit(main())
