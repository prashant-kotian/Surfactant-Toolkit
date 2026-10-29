# Reproducibility and hallucination sets

## R-01..R-18 (internal GZ-51..GZ-68), five answers each, unaugmented

Questions: `benchmark/source_sets/groundzero/reproducibility_bank_x5.txt`. These 18 questions are
not part of SurfBench-200.

    python extract_runs.py   # transcripts/ -> runs.json (each answer assigned by its unique answer keys)
    python score_runs.py     # reference rules + hallucination codes -> scored.json

`reproducibility_answers.csv` (one row per answer) and `repro_summary.json` are the tables built
from `scored.json` for the paper.

| Model | Protocol | Transcripts (`transcripts/`) |
|---|---|---|
| Claude | one chat per run per block of questions | `unagumented 18th sept/question repeat type/` (GZ-51..60), `unaggumented output 19th sept/` (GZ-61..65; `New folder/` GZ-66..68) |
| ChatGPT | new session per prompt | `unaugmented_runs_20260921/chatgpt_GZ51-68_repro_x5_unaugmented.txt` |
| Gemini | one chat per run | `new output/New folder/reproducibility type/gemini/` (the file named "GZ-69-78" is run 3 of GZ-51..68) |
| Qwen | one chat per run | `new output/New folder/reproducibility type/qween/` |

Each model ran in the same version and setting as in the SurfBench-200 runs.

`unaugmented_runs_20260921/gemini_GZ51-68_repro_x5_unaugmented.txt` came from a signed-out Gemini
session (Flash-Lite), not the benchmark's Gemini setting; it is extracted into runs.json for the
record and never scored. File and folder names are those of the original exports.

Outcomes: correct / acceptable (judgement questions R-02, R-05, R-06 intermediate, R-10 in
1.5-4 mM, R-14) / error / hallucination / missing. Hallucination types: F value for a withheld
input, C stated condition overridden, U unsupported claim (R-17 coded from the full text, table
`CODES_67` in `score_runs.py`), N value with no data to read.

## Eleven SurfBench-200 questions, five answers each (`sb_extension/`)

SB-011, SB-012, SB-017, SB-038, SB-042, SB-061, SB-065, SB-069, SB-070, SB-073 and SB-169, asked
five times to each web model, unaugmented, every prompt in a new chat. Question texts are in
`sb_extension/questions/`; answers in `sb_extension/runs/r<N>/<model>_unaug/<QID>.txt` (the
FINAL_JSON block) with `<QID>.meta.json` (model, setting and tool state verified on the page before
sending). Answers collected by hand are in `manual_answers/` (prompts in `manual_prompts/`);
hand-collected runs beyond the five per question are kept in `manual_answers/extra/` and not scored.

    python sb_extension/score_ext.py      # -> scored_ext.json
    python -m pytest sb_extension/test_score_ext.py

`benchmark_repeats_answers.csv` and `ext_summary.json` are the tables built from `scored_ext.json`.
