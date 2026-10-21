# Per-source-set results (internal diagnostic view)

> These are the authoring groups behind the single SurfBench-200 benchmark;
> the headline table for reporting is `results/TABLE200.md`. (generated — rerun
> benchmark/build_tables.py)

## Unaugmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 33 | 49 | 67.3 |  |
| GroundZero-50 | ChatGPT | 31 | 49 | 63.3 |  |
| GroundZero-50 | DeepSeek | 32 | 48 | 66.7 | partial coverage (n=48) |
| GroundZero-50 | Gemini (web) | 15 | 40 | 37.5 | partial coverage (n=40) |
| GroundZero-50 | Qwen | 31 | 49 | 63.3 |  |
| GroundZero-50 | Antigravity CLI | 48 | 49 | 98.0 |  |
| Hard50 (new) | Claude | 36 | 50 | 72.0 |  |
| Hard50 (new) | ChatGPT | 36 | 50 | 72.0 | count-only |
| Hard50 (new) | DeepSeek | 34 | 50 | 68.0 | count-only |
| Hard50 (new) | Gemini (web) | 31 | 47 | 66.0 | count-only |
| Hard50 (new) | Qwen | 33 | 50 | 66.0 | count-only |
| Hard50 (new) | Antigravity CLI | 48 | 50 | 96.0 |  |
| ExtraHard-50 | Claude | 33 | 50 | 66.0 | count-only |
| ExtraHard-50 | ChatGPT | 33 | 50 | 66.0 | count-only |
| ExtraHard-50 | DeepSeek | 28 | 50 | 56.0 | count-only |
| ExtraHard-50 | Gemini (web) | 26 | 49 | 53.1 | count-only |
| ExtraHard-50 | Qwen | 21 | 50 | 42.0 | count-only |
| ExtraHard-50 | Antigravity CLI | 23 | 47 | 48.9 | pro-high tier (retired) — disclosed, final |

## Augmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 42 | 47 | 89.4 |  |
| GroundZero-50 | ChatGPT | 39 | 49 | 79.6 |  |
| GroundZero-50 | DeepSeek | 39 | 46 | 84.8 | SI-only; partial coverage (n=46) |
| GroundZero-50 | Antigravity CLI | 47 | 49 | 95.9 |  |
| Hard50 (new) | Claude | 38 | 40 | 95.0 |  |
| Hard50 (new) | ChatGPT | 47 | 50 | 94.0 |  |
| Hard50 (new) | DeepSeek | 46 | 50 | 92.0 | SI-only |
| Hard50 (new) | Antigravity CLI | 50 | 50 | 100.0 |  |
| ExtraHard-50 | Claude | 48 | 48 | 100.0 | count-only |
| ExtraHard-50 | ChatGPT | 50 | 50 | 100.0 | count-only |
| ExtraHard-50 | DeepSeek | 47 | 50 | 94.0 | count-only |
| ExtraHard-50 | Antigravity CLI | 47 | 47 | 100.0 | pro-high tier (retired) — disclosed, final |

## Augmentation lift

| Source set | Model | unaug % | aug % | Δ pts | comparability |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 67.3 | 89.4 | +22.1 | ok |
| GroundZero-50 | ChatGPT | 63.3 | 79.6 | +16.3 | ok |
| GroundZero-50 | DeepSeek | 66.7 | 84.8 | +18.1 | ok |
| GroundZero-50 | Antigravity CLI | 98.0 | 95.9 | -2.1 | ok |
| Hard50 (new) | Claude | 72.0 | 95.0 | +23.0 | ok |
| Hard50 (new) | ChatGPT | 72.0 | 94.0 | +22.0 | harness differs |
| Hard50 (new) | DeepSeek | 68.0 | 92.0 | +24.0 | harness differs |
| Hard50 (new) | Antigravity CLI | 96.0 | 100.0 | +4.0 | ok |
| ExtraHard-50 | Claude | 66.0 | 100.0 | +34.0 | ok |
| ExtraHard-50 | ChatGPT | 66.0 | 100.0 | +34.0 | ok |
| ExtraHard-50 | DeepSeek | 56.0 | 94.0 | +38.0 | ok |
| ExtraHard-50 | Antigravity CLI | 48.9 | 100.0 | +51.1 | pro-high tier |

## Declared exclusions

| Source set | n | n cells | excluded | reason |
|---|---|---|---|---|
| GroundZero-50 | 50 | 40/46/47/48/49 | GZ-16, GZ-20, GZ-26 | question requires image/chart input the harnesses cannot deliver |
| Hard50 (new) | 50 | 40/47/50 | H-02, H-03, H-08, H-16, H-19, H-31, H-38, H-46, H-47, H-49 | image/figure-based questions excluded from all cells |
| ExtraHard-50 | 50 | 47/48/49/50 | X-19, X-41, X-47 | question body carried a kappa*a typo (1.45 vs 14.5) at run time for the pre-fix packets (fixed in-bank 2026-10-16); excluded from every arm so unaugmented/augmented pairs share one question set |

## Provenance ledger

| cell | passed/n | % | harness | source | sha | gold | notes |
|---|---|---|---|---|---|---|---|
| extrahard / ChatGPT / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json | bb19fd79144a | post-fix | count-only |
| extrahard / ChatGPT / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Claude / aug | 48/48 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 03c5689f29cd | post-fix | count-only |
| extrahard / Claude / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / DeepSeek / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 03c5689f29cd | post-fix | count-only |
| extrahard / DeepSeek / unaug | 28/50 | 56.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Gemini-CLI / aug | 47/47 | 100.0 | agy-cli | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json | 6043b0793db0 | post-fix | pro-high tier (retired) — disclosed, final |
| extrahard / Gemini-CLI / unaug | 23/47 | 48.9 | agy-cli | benchmark/source_sets/extrahard/antigravity_unaugmented_grade_results.json | a22ae6e7d157 | post-fix | pro-high tier (retired) — disclosed, final |
| extrahard / Gemini / unaug | 26/49 | 53.1 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Qwen / unaug | 21/50 | 42.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| groundzero / ChatGPT / aug | 39/49 | 79.6 | mcp | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json | 26bf68e59aa6 | post-fix |  |
| groundzero / ChatGPT / unaug | 31/49 | 63.3 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix |  |
| groundzero / Claude / aug | 42/47 | 89.4 | mcp | benchmark/source_sets/groundzero/claude_augmented_grade_results.json | a0bc7ccb7625 | post-fix |  |
| groundzero / Claude / unaug | 33/49 | 67.3 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix |  |
| groundzero / DeepSeek / aug | 39/46 | 84.8 | mcp | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json | 42ac653966c4 | post-fix | SI-only; partial coverage (n=46) |
| groundzero / DeepSeek / unaug | 32/48 | 66.7 | mcp | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json | 2e8f8288e360 | post-fix | partial coverage (n=48) |
| groundzero / Gemini-CLI / aug | 47/49 | 95.9 | agy-cli | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json | 5a753572989a | post-fix |  |
| groundzero / Gemini-CLI / unaug | 48/49 | 98.0 | agy-cli | benchmark/source_sets/groundzero/antigravity_unaugmented_grade_results.json | c58703abce24 | post-fix |  |
| groundzero / Gemini / unaug | 15/40 | 37.5 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix | partial coverage (n=40) |
| groundzero / Qwen / unaug | 31/49 | 63.3 | mcp | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json | f8b8615eb935 | post-fix |  |
| hard50 / ChatGPT / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json | 33cd934a4128 | post-fix |  |
| hard50 / ChatGPT / unaug | 36/50 | 72.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Claude / aug | 38/40 | 95.0 | mcp | benchmark/source_sets/hard50/claude_augmented_grade_results.json | ecbb39aa68e1 | post-fix |  |
| hard50 / Claude / unaug | 36/50 | 72.0 | mcp | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json | 03ebde500fa9 | post-fix |  |
| hard50 / DeepSeek / aug | 46/50 | 92.0 | mcp | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json | 2ae21eb719d3 | post-fix | SI-only |
| hard50 / DeepSeek / unaug | 34/50 | 68.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Gemini-CLI / aug | 50/50 | 100.0 | agy-cli | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json | bedea6408f77 | post-fix |  |
| hard50 / Gemini-CLI / unaug | 48/50 | 96.0 | agy-cli | benchmark/source_sets/hard50/antigravity_unaugmented_grade_results.json | bde27dffc7f2 | post-fix |  |
| hard50 / Gemini / unaug | 31/47 | 66.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Qwen / unaug | 33/50 | 66.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
