"""Write gold JSON files back in each file's own original layout, so a correction diff
shows only the changed values (2026-10-20)."""
import json
from pathlib import Path

STYLE = {  # bank -> json.dumps kwargs reproducing the committed file byte-for-byte
    "groundzero": dict(indent=2, ensure_ascii=False),
    "hard50": dict(indent=2, ensure_ascii=False),
    "extrahard": dict(indent=1, sort_keys=True, ensure_ascii=False),
}


def dump(bank: str, obj) -> str:
    return json.dumps(obj, **STYLE[bank])


def write(path: Path, bank: str, obj) -> None:
    path.write_text(dump(bank, obj), encoding="utf-8")
