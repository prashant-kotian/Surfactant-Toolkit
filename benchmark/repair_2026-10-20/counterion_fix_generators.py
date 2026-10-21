"""Apply the counterion-factor correction to the question/gold GENERATORS (2026-10-20),
so a rebuild cannot reintroduce the inverse convention. H-04 is deliberately untouched
except for its label (its 0.268 is an ionization degree; 2 - alpha was right).
Usage: python counterion_fix_generators.py [--apply]
"""
import re
import sys
from pathlib import Path

B = Path(__file__).resolve().parents[1]
DRY = "--apply" not in sys.argv

EDITS = {
    B / "source_sets/hard50/gold_build.py": [
        (r"counterion_factor=2\.0 - (rc?\.counterion_binding_degree)", r"counterion_factor=1.0 + \1"),
        (r"counterion_factor=1\.0 \+ (H03_ALPHA|H19_ALPHA|H31_ALPHA|alpha)\b", r"counterion_factor=2.0 - \1"),
        (r"factor = 1\.0 \+ H49_ALPHA", "factor = 2.0 - H49_ALPHA"),
    ],
    B / "source_sets/hard50/hard50_questions.py": [
        (r"counterion binding degree beta = 0\.268 \(real value from the same study\)",
         "degree of counterion ionization alpha = 0.268 (S2/S1; real value from the same study)"),
        (r"dG = \(2 - beta\)\*R\*T\*ln\(CMC\[mol/L\]/55\.5\), then extract",
         "dG = (2 - alpha)*R*T*ln(CMC[mol/L]/55.5), then extract"),           # H-04 only
        (r"\(2 - beta\)", "(1 + beta)"),
        (r"\(1 \+ alpha\)", "(2 - alpha)"),
        (r"build dG\(T\) rows with \(2-beta\)", "build dG(T) rows with (2-alpha)"),  # H-04 guidance
        (r"\(2-beta\)", "(1+beta)"),
        (r"\(1\+alpha\)", "(2-alpha)"),
        (r"or 2-beta used with alpha", "or the inverse convention (2-beta / 1+alpha) used"),
    ],
    B / "source_sets/hard50/build_hard50.py": [
        (r"build dG\(T\) rows with \(2-beta\)=1\.732", "build dG(T) rows with (2-alpha)=1.732"),
        (r"\(2-beta\)", "(1+beta)"),
        (r"\(1\+alpha\)", "(2-alpha)"),
    ],
    B / "source_sets/extrahard/xb_c5.py": [
        (r"gibbs_free_energy_micellization\(_X\(cmc\), T, 2\.0 - beta\)",
         "gibbs_free_energy_micellization(_X(cmc), T, 1.0 + beta)"),
        (r"gibbs_free_energy_micellization\(_X\(rc\.cmc_mM\), 298\.15, 2\.0 - beta\)",
         "gibbs_free_energy_micellization(_X(rc.cmc_mM), 298.15, 1.0 + beta)"),
        (r"\(2 - beta\)", "(1 + beta)"),
    ],
    B / "mcp_harness/tools.py": [
        (r"1\.0 nonionic; \(2-beta\) for an ionic surfactant",
         "1.0 nonionic; (1 + beta) for an ionic surfactant, beta = counterion binding (= 2 - alpha, alpha = S2/S1)"),
    ],
}


def main():
    for path, subs in EDITS.items():
        text = path.read_text(encoding="utf-8")
        new = text
        counts = []
        for pat, rep in subs:
            new, n = re.subn(pat, rep, new)
            counts.append(n)
        print(f"{path.relative_to(B)}: {counts}")
        if not DRY and new != text:
            path.write_text(new, encoding="utf-8")
    print("DRY RUN" if DRY else "APPLIED")


if __name__ == "__main__":
    main()
