# SurfBench-200 — grand results

Single benchmark, 200 questions (the four authoring sets are internal provenance only).
`n` = gradable questions after the declared exclusion policy (see DATA_NOTES.md);
`Aug` = toolkit/MCP-augmented run, `Una` = unaugmented baseline. Lift = Aug − Una (pp).

| Model | Unaugmented | Augmented | Lift | Gradable n /200 | Flags |
|---|---|---|---|---|---|
| Claude | 129/200 (64.5%) | 192/200 (96.0%) | +31.5 pp | 200 | extrahard: count-only grading (both arms) |
| ChatGPT | 128/200 (64.0%) | 187/200 (93.5%) | +29.5 pp | 200 | hard50: count-only grading; extrahard: count-only grading (both arms) |
| DeepSeek | 122/200 (61.0%) | 184/200 (92.0%) | +31.0 pp | 200 | hard50: count-only grading; extrahard: count-only grading (both arms); groundzero: SI-only grading; hard50: SI-only grading |
| Gemini | 101/200 (50.5%) | 198/200 (99.0%) | +48.5 pp | 200 | unaugmented = Gemini web; hard50: count-only grading; extrahard: count-only grading; augmented = Gemini CLI (Antigravity) |
| Qwen | 106/200 (53.0%) | — | — | 200 | hard50: count-only grading; extrahard: count-only grading |

## Per-authoring-set breakdown (provenance; not separate benchmarks)

| Model | groundzero U/A | hard50 U/A | extrahard U/A | orig_hard50 U/A |
|---|---|---|---|---|---|
| Claude | 34/50 · 44/50 | 36/50 · 48/50 | 33/50 · 50/50 | 26/50 · 50/50 |
| ChatGPT | 32/50 · 40/50 | 36/50 · 47/50 | 33/50 · 50/50 | 27/50 · 50/50 |
| DeepSeek | 34/50 · 42/50 | 34/50 · 46/50 | 28/50 · 47/50 | 26/50 · 49/50 |
| Gemini | 19/50 · 48/50 | 32/50 · 50/50 | 27/50 · 50/50 | 23/50 · 50/50 |
| Qwen | 32/50 · — | 33/50 · — | 21/50 · — | 20/50 · — |
