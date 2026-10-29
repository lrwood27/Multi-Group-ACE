"""Fetch checksum-pinned authentic neutron ACE fixtures into ignored storage.

Run: python tests/fetch_ce_data.py
No data are downloaded during ordinary test discovery.
"""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).with_name("ce_data_manifest.json")


def matches(path, entry):
    return (
        path.is_file()
        and path.stat().st_size == entry["size"]
        and hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    )


def fetch(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for entry in json.loads(MANIFEST.read_text())["files"]:
        target = output_dir / entry["name"]
        if matches(target, entry):
            print(f"Verified: {target}")
            continue
        staged = None
        try:
            with NamedTemporaryFile(dir=output_dir, delete=False) as output:
                staged = Path(output.name)
                with urlopen(entry["url"], timeout=120) as source:
                    while chunk := source.read(1024 * 1024):
                        output.write(chunk)
            if not matches(staged, entry):
                raise ValueError(f"Fixture checksum/size mismatch for {entry['name']}")
            staged.replace(target)
            print(f"Downloaded and verified: {target}")
        finally:
            if staged is not None:
                staged.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / ".test-data" / "ce")
    fetch(parser.parse_args().outdir)
