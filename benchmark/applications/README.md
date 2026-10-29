# Application problems P1..P7: statements, transcripts and toolkit reference calculations

- `problem_statements.md` -- the seven problems as authored, copied unchanged. Each model received
  only the shared preamble and the PROBLEM block; the "ANALYST APPENDIX" at the end was never
  shown to a model. File names it mentions refer to the development repository.
- `transcripts/` -- the model answers, exported from the chat interfaces and copied unchanged:

  | File | Model | Condition | Problems |
  |---|---|---|---|
  | `chatgpt_augmented_P1-P7.txt` | ChatGPT 5.6 | toolkit MCP server attached | P1..P7 |
  | `chatgpt_unaugmented_P1-P7.txt` | ChatGPT 5.6 | no toolkit | P1..P7 |
  | `claude_augmented_P1-P7.txt` | Claude Opus 4.8 | toolkit MCP server attached | P1..P7 |
  | `claude_unaugmented_P1-P5.txt` | Claude Opus 4.8 | no toolkit | P1..P5 |

- `reference_2026-10-05_results_B_C_5_6_7.json` -- the toolkit reference calculations made when the
  problems were written.
- `compute_p2_p6.py`, `compute_p3.py`, `compute_p4.py` (outputs `values_p2_p6.json`, `values_p3.json`, `values_p4.json`) -- the toolkit values
  quoted in the paper for P2/P6, P3 and P4, recomputed on this release:

      python benchmark/applications/compute_p2_p6.py
      python benchmark/applications/compute_p3.py
      python benchmark/applications/compute_p4.py

  Each script uses this repository's `src/` (or a toolkit path given as the first argument).
