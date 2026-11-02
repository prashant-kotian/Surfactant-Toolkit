# research/ — pre-publication research code

This directory is **not** part of the `surfactantkit` package. Nothing here is
imported by `src/`, shipped in the distribution, exposed through the MCP server,
or covered by the validation tiers in `src/surfactantkit/validation_tiers.py`.

It exists for one purpose: to let a claim be **tested and, if necessary, killed**
before it becomes a promise in a paper. Code promoted out of here must first earn
a validation tier and a place in `tests/`.

## Live thread

| Thread | Question | Status |
|---|---|---|
| `laod_tau2/` | Can the slow micellar relaxation time τ₂ be recovered from the harmonic content of a large-amplitude oscillatory dilatational (LAOD) response? | Kill-criterion study **complete**; viable in a quantified bandwidth window. See `laod_tau2/FINDINGS.md` |

## Running it

```bash
pip install -e ".[calc,dev]"          # numpy, scipy, pytest, rdkit
PYTHONPATH=src:. python -m pytest research/laod_tau2/tests -q
PYTHONPATH=src:. python -m research.laod_tau2.study_bandwidth_rule
```

The bandwidth study takes ~25 minutes; the anchors take seconds. **Run the
anchors first** — if they fail, every downstream number is meaningless.

## House rules

1. Analytic anchor before numerical result. A model that cannot reproduce
   Lucassen–van den Tempel in the small-amplitude limit is wrong, regardless of
   how plausible its output looks.
2. State the kill criterion before running the study that could fail it.
3. Negative results get written down (`FINDINGS.md` §6), including the ones that
   were our own mistakes.
4. Parameter regimes must be justified from independent data before a sensitivity
   study is run in them.
