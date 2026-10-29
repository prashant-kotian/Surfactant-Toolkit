# SB-151..SB-200 (orig_hard50): first-release transcripts

Raw answers of the project's first benchmark release for the fourth authoring set, copied
unchanged from the development repository (`benchmark/blind_run/`, `benchmark/augmented_run/`):

- `blind_run/<model>_blind_Q01-Q60.txt` -- unaugmented answers (wording:
  `benchmark/source_sets/orig_hard50/question_bank_unaugmented.txt`), Q01..Q50 = SB-151..SB-200;
  Q51..Q60 were not part of SurfBench-200.
- `blind_run/rectification_reruns/` -- answers to Q17, Q31, Q34, Q46 and Q59 after their wording was
  clarified on 2026-10-04 (the bracketed notes in the question files).
- `augmented_run/` -- augmented answers and session logs (wording: `orig_hard50/Q*.md`).
- `grade_blind_hard60.py`, `grade_augmented_hard60.py`, `answers_data.py`,
  `augmented_answers.py` -- the first-release grading scripts and parsed answers.

The verdicts used in the results are the frozen grading files in
`benchmark/source_sets/orig_hard50/*_grade_results.json`; `benchmark/repair_2026-10-20/
orig_hard50_grade.py` reproduces all 250 of them from the same rule.
