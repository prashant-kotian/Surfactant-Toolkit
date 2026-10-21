# SurfBench-200 — grand results

Single benchmark, 200 questions (the four authoring sets are internal provenance only).
`n` = gradable questions after the declared exclusion policy (see DATA_NOTES.md);
`Aug` = toolkit/MCP-augmented run, `Una` = unaugmented baseline. Lift = Aug − Una (pp).

| Model | Unaugmented | Augmented | Lift | Gradable n /200 | Flags |
|---|---|---|---|---|---|
| Claude | 127/198 (64.1%) | 178/185 (96.2%) | +32.1 pp | 198 | extrahard: count-only grading (both arms); orig_hard50: 1 unanswered (ungradeable, excluded from n) |
| ChatGPT | 127/195 (65.1%) | 186/199 (93.5%) | +28.4 pp | 195 | hard50: count-only grading; extrahard: count-only grading (both arms); orig_hard50: 4 unanswered (ungradeable, excluded from n) |
| DeepSeek | 120/198 (60.6%) | 181/196 (92.3%) | +31.7 pp | 198 | groundzero: partial transcript (both arms); hard50: count-only grading; extrahard: count-only grading (both arms); hard50: SI-only grading |
| Qwen | 105/199 (52.8%) | — | — | 199 | hard50: count-only grading; extrahard: count-only grading |
| Gemini | 95/186 (51.1%) | — | — | 186 | groundzero: partial transcript; hard50: count-only grading; extrahard: count-only grading |
| Gemini-CLI | 119/146 (81.5%) | 194/196 (99.0%) | +17.5 pp | 146 | — |

## Per-authoring-set breakdown (provenance; not separate benchmarks)

| Model | groundzero U/A | hard50 U/A | extrahard U/A | orig_hard50 U/A |
|---|---|---|---|---|---|
| Claude | 33/49 · 42/47 | 36/50 · 38/40 | 33/50 · 48/48 | 25/49 · 50/50 |
| ChatGPT | 31/49 · 39/49 | 36/50 · 47/50 | 33/50 · 50/50 | 27/46 · 50/50 |
| DeepSeek | 32/48 · 39/46 | 34/50 · 46/50 | 28/50 · 47/50 | 26/50 · 49/50 |
| Qwen | 31/49 · — | 33/50 · — | 21/50 · — | 20/50 · — |
| Gemini | 15/40 · — | 31/47 · — | 26/49 · — | 23/50 · — |
| Gemini-CLI | 48/49 · 47/49 | 48/50 · 50/50 | 23/47 · 47/47 | — · 50/50 |
