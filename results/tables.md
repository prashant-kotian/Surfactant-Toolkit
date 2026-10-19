# Per-source-set results (internal diagnostic view)

> These are the authoring groups behind the single SurfBench-200 benchmark;
> the headline table for reporting is `results/TABLE200.md`. (generated — rerun
> benchmark/build_tables.py)

## Unaugmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 5 | 10 | 50.0 | partial coverage (n=10) |
| GroundZero-50 | ChatGPT | 22 | 43 | 51.2 |  |
| GroundZero-50 | DeepSeek | 29 | 42 | 69.0 | partial coverage (n=42) |
| GroundZero-50 | Gemini (web) | 10 | 33 | 30.3 | partial coverage (n=33) |
| GroundZero-50 | Qwen | 27 | 43 | 62.8 |  |
| GroundZero-50 | Antigravity CLI | 42 | 43 | 97.7 |  |
| Hard50 (new) | Claude | 33 | 47 | 70.2 |  |
| Hard50 (new) | ChatGPT | 34 | 47 | 72.3 | count-only |
| Hard50 (new) | DeepSeek | 31 | 47 | 66.0 | count-only |
| Hard50 (new) | Gemini (web) | 29 | 47 | 61.7 | count-only |
| Hard50 (new) | Qwen | 30 | 47 | 63.8 | count-only |
| Hard50 (new) | Antigravity CLI | 46 | 47 | 97.9 |  |
| ExtraHard-50 | Claude | 33 | 49 | 67.3 | count-only |
| ExtraHard-50 | ChatGPT | 32 | 49 | 65.3 | count-only |
| ExtraHard-50 | DeepSeek | 28 | 49 | 57.1 | count-only |
| ExtraHard-50 | Gemini (web) | 27 | 49 | 55.1 | count-only |
| ExtraHard-50 | Qwen | 21 | 49 | 42.9 | count-only |
| ExtraHard-50 | Antigravity CLI | 25 | 49 | 51.0 | pro-high tier (retired) — disclosed, final |

## Augmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 36 | 43 | 83.7 |  |
| GroundZero-50 | ChatGPT | 31 | 43 | 72.1 |  |
| GroundZero-50 | DeepSeek | 32 | 40 | 80.0 | SI-only; partial coverage (n=40) |
| GroundZero-50 | Antigravity CLI | 41 | 43 | 95.3 |  |
| Hard50 (new) | Claude | 45 | 47 | 95.7 |  |
| Hard50 (new) | ChatGPT | 44 | 47 | 93.6 |  |
| Hard50 (new) | DeepSeek | 44 | 47 | 93.6 | SI-only |
| Hard50 (new) | Antigravity CLI | 47 | 47 | 100.0 |  |
| ExtraHard-50 | Claude | 50 | 50 | 100.0 | count-only |
| ExtraHard-50 | ChatGPT | 50 | 50 | 100.0 | count-only |
| ExtraHard-50 | DeepSeek | 48 | 50 | 96.0 | count-only |
| ExtraHard-50 | Antigravity CLI | 49 | 49 | 100.0 | pro-high tier (retired) — disclosed, final |

## Augmentation lift

| Source set | Model | unaug % | aug % | Δ pts | comparability |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 50.0 | 83.7 | +33.7 | ok |
| GroundZero-50 | ChatGPT | 51.2 | 72.1 | +20.9 | ok |
| GroundZero-50 | DeepSeek | 69.0 | 80.0 | +11.0 | ok |
| GroundZero-50 | Antigravity CLI | 97.7 | 95.3 | -2.4 | ok |
| Hard50 (new) | Claude | 70.2 | 95.7 | +25.5 | ok |
| Hard50 (new) | ChatGPT | 72.3 | 93.6 | +21.3 | harness differs |
| Hard50 (new) | DeepSeek | 66.0 | 93.6 | +27.6 | harness differs |
| Hard50 (new) | Antigravity CLI | 97.9 | 100.0 | +2.1 | ok |
| ExtraHard-50 | Claude | 67.3 | 100.0 | +32.7 | ok |
| ExtraHard-50 | ChatGPT | 65.3 | 100.0 | +34.7 | ok |
| ExtraHard-50 | DeepSeek | 57.1 | 96.0 | +38.9 | ok |
| ExtraHard-50 | Antigravity CLI | 51.0 | 100.0 | +49.0 | pro-high tier |

## Declared exclusions

| Source set | n | n cells | excluded | reason |
|---|---|---|---|---|
| GroundZero-50 | 50 | 10/33/40/42/43 | GZ-15, GZ-26, GZ-36, GZ-37, GZ-38, GZ-39, GZ-40 | question requires image/chart input the harnesses cannot deliver |
| Hard50 (new) | 50 | 47 | H-36, H-37, H-38 | image/figure-based questions excluded from all cells |
| ExtraHard-50 | 50 | 49/50 | X-19 | question body carried a kappa*a typo (1.45 vs 14.5) at run time for the pre-fix packets; fixed in-bank 2026-10-16, later-run cells keep it |

## Provenance ledger

| cell | passed/n | % | harness | source | sha | gold | notes |
|---|---|---|---|---|---|---|---|
| extrahard / ChatGPT / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json | 009e8d7f512c | post-fix | count-only |
| extrahard / ChatGPT / unaug | 32/49 | 65.3 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | a23be00595f4 | post-fix | count-only |
| extrahard / Claude / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | f834e22a4171 | post-fix | count-only |
| extrahard / Claude / unaug | 33/49 | 67.3 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | a23be00595f4 | post-fix | count-only |
| extrahard / DeepSeek / aug | 48/50 | 96.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | f834e22a4171 | post-fix | count-only |
| extrahard / DeepSeek / unaug | 28/49 | 57.1 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | a23be00595f4 | post-fix | count-only |
| extrahard / Gemini-CLI / aug | 49/49 | 100.0 | agy-cli | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json | 056d793768c9 | post-fix | pro-high tier (retired) — disclosed, final |
| extrahard / Gemini-CLI / unaug | 25/49 | 51.0 | agy-cli | benchmark/source_sets/extrahard/antigravity_unaugmented_grade_results.json | 48fb70a0410c | post-fix | pro-high tier (retired) — disclosed, final |
| extrahard / Gemini / unaug | 27/49 | 55.1 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | a23be00595f4 | post-fix | count-only |
| extrahard / Qwen / unaug | 21/49 | 42.9 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | a23be00595f4 | post-fix | count-only |
| groundzero / ChatGPT / aug | 31/43 | 72.1 | mcp | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json | acecfd25766d | post-fix |  |
| groundzero / ChatGPT / unaug | 22/43 | 51.2 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | 54beb1e5502c | post-fix |  |
| groundzero / Claude / aug | 36/43 | 83.7 | mcp | benchmark/source_sets/groundzero/claude_augmented_grade_results.json | 0e391e177120 | post-fix |  |
| groundzero / Claude / unaug | 5/10 | 50.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | 54beb1e5502c | post-fix | partial coverage (n=10) |
| groundzero / DeepSeek / aug | 32/40 | 80.0 | mcp | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json | 8537cdd2b24d | post-fix | SI-only; partial coverage (n=40) |
| groundzero / DeepSeek / unaug | 29/42 | 69.0 | mcp | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json | f28606937ad9 | post-fix | partial coverage (n=42) |
| groundzero / Gemini-CLI / aug | 41/43 | 95.3 | agy-cli | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json | c08a4aae10a4 | post-fix |  |
| groundzero / Gemini-CLI / unaug | 42/43 | 97.7 | agy-cli | benchmark/source_sets/groundzero/antigravity_unaugmented_grade_results.json | 3637904c8a02 | post-fix |  |
| groundzero / Gemini / unaug | 10/33 | 30.3 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | 54beb1e5502c | post-fix | partial coverage (n=33) |
| groundzero / Qwen / unaug | 27/43 | 62.8 | mcp | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json | 5b42285f44a8 | post-fix |  |
| hard50 / ChatGPT / aug | 44/47 | 93.6 | mcp | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json | 4e6f78c94eb5 | post-fix |  |
| hard50 / ChatGPT / unaug | 34/47 | 72.3 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 88591d6d66b4 | post-fix | count-only |
| hard50 / Claude / aug | 45/47 | 95.7 | mcp | benchmark/source_sets/hard50/claude_augmented_grade_results.json | d86f784506b9 | post-fix |  |
| hard50 / Claude / unaug | 33/47 | 70.2 | mcp | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json | 42c45e2da3cd | post-fix |  |
| hard50 / DeepSeek / aug | 44/47 | 93.6 | mcp | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json | 36041a918ab4 | post-fix | SI-only |
| hard50 / DeepSeek / unaug | 31/47 | 66.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 88591d6d66b4 | post-fix | count-only |
| hard50 / Gemini-CLI / aug | 47/47 | 100.0 | agy-cli | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json | 09e72354d484 | post-fix |  |
| hard50 / Gemini-CLI / unaug | 46/47 | 97.9 | agy-cli | benchmark/source_sets/hard50/antigravity_unaugmented_grade_results.json | 198dda619bdd | post-fix |  |
| hard50 / Gemini / unaug | 29/47 | 61.7 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 88591d6d66b4 | post-fix | count-only |
| hard50 / Qwen / unaug | 30/47 | 63.8 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 88591d6d66b4 | post-fix | count-only |
