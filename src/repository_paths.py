"""Resolve relocated files when reading historical provenance records."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def resolve_recorded_path(path):
    """Map only explicitly relocated files; never replace an existing file."""
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
    if path.exists():
        return path
    try:
        relative = path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path
    moves = json.loads((ROOT / 'docs/provenance/relocated_files.json').read_text())
    return ROOT / moves[relative] if relative in moves else path
