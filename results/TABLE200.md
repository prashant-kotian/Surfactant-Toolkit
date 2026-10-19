# SurfBench-200 — grand results

Single benchmark, 200 questions (the four authoring sets are internal provenance only).
`n` = gradable questions after the declared exclusion policy (see DATA_NOTES.md);
`Aug` = toolkit/MCP-augmented run, `Una` = unaugmented baseline. Lift = Aug − Una (pp).

| Model | Unaugmented | Augmented | Lift | Gradable n /200 | Flags |
|---|---|---|---|---|---|
| Claude | 118/188 (62.8%) | 181/190 (95.3%) | +32.5 pp | 188 | extrahard: count-only grading (both arms); orig_hard50: 1 unanswered (ungradeable, excluded from n) |
| ChatGPT | 120/185 (64.9%) | 175/190 (92.1%) | +27.2 pp | 185 | hard50: count-only grading; extrahard: count-only grading (both arms); orig_hard50: 4 unanswered (ungradeable, excluded from n) |
| DeepSeek | 114/188 (60.6%) | 173/187 (92.5%) | +31.9 pp | 188 | groundzero: partial transcript (both arms); hard50: count-only grading; extrahard: count-only grading (both arms); hard50: SI-only grading |
| Qwen | 98/189 (51.9%) | — | — | 189 | hard50: count-only grading; extrahard: count-only grading |
| Gemini | 89/179 (49.7%) | — | — | 179 | groundzero: partial transcript; hard50: count-only grading; extrahard: count-only grading |
| Gemini-CLI | 113/139 (81.3%) | 187/189 (98.9%) | +17.6 pp | 139 | — |

## Per-authoring-set breakdown (provenance; not separate benchmarks)

| Model | groundzero U/A | hard50 U/A | extrahard U/A | orig_hard50 U/A |
|---|---|---|---|---|---|
| Claude | 27/43 · 36/43 | 33/47 · 45/47 | 33/49 · 50/50 | 25/49 · 50/50 |
| ChatGPT | 27/43 · 31/43 | 34/47 · 44/47 | 32/49 · 50/50 | 27/46 · 50/50 |
| DeepSeek | 29/42 · 32/40 | 31/47 · 44/47 | 28/49 · 48/50 | 26/50 · 49/50 |
| Qwen | 27/43 · — | 30/47 · — | 21/49 · — | 20/50 · — |
| Gemini | 10/33 · — | 29/47 · — | 27/49 · — | 23/50 · — |
| Gemini-CLI | 42/43 · 41/43 | 46/47 · 47/47 | 25/49 · 49/49 | — · 50/50 |
