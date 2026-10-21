# Per-source-set results (internal diagnostic view)

> These are the authoring groups behind the single SurfBench-200 benchmark;
> the headline table for reporting is `results/TABLE200.md`. (generated — rerun
> benchmark/build_tables.py)

## Unaugmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 34 | 50 | 68.0 |  |
| GroundZero-50 | ChatGPT | 32 | 50 | 64.0 |  |
| GroundZero-50 | DeepSeek | 34 | 50 | 68.0 |  |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 20 | 50 | 40.0 |  |
| GroundZero-50 | Qwen | 32 | 50 | 64.0 |  |
| Hard50 (new) | Claude | 36 | 50 | 72.0 |  |
| Hard50 (new) | ChatGPT | 36 | 50 | 72.0 | count-only |
| Hard50 (new) | DeepSeek | 34 | 50 | 68.0 | count-only |
| Hard50 (new) | Gemini (unaug: web; aug: CLI) | 32 | 50 | 64.0 | count-only |
| Hard50 (new) | Qwen | 33 | 50 | 66.0 | count-only |
| ExtraHard-50 | Claude | 33 | 50 | 66.0 | count-only |
| ExtraHard-50 | ChatGPT | 33 | 50 | 66.0 | count-only |
| ExtraHard-50 | DeepSeek | 28 | 50 | 56.0 | count-only |
| ExtraHard-50 | Gemini (unaug: web; aug: CLI) | 27 | 50 | 54.0 | count-only |
| ExtraHard-50 | Qwen | 21 | 50 | 42.0 | count-only |

## Augmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 44 | 50 | 88.0 |  |
| GroundZero-50 | ChatGPT | 40 | 50 | 80.0 |  |
| GroundZero-50 | DeepSeek | 42 | 50 | 84.0 | SI-only |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 48 | 50 | 96.0 |  |
| Hard50 (new) | Claude | 48 | 50 | 96.0 |  |
| Hard50 (new) | ChatGPT | 47 | 50 | 94.0 |  |
| Hard50 (new) | DeepSeek | 46 | 50 | 92.0 | SI-only |
| Hard50 (new) | Gemini (unaug: web; aug: CLI) | 50 | 50 | 100.0 |  |
| ExtraHard-50 | Claude | 50 | 50 | 100.0 | count-only |
| ExtraHard-50 | ChatGPT | 50 | 50 | 100.0 | count-only |
| ExtraHard-50 | DeepSeek | 47 | 50 | 94.0 | count-only |
| ExtraHard-50 | Gemini (unaug: web; aug: CLI) | 50 | 50 | 100.0 | Gemini CLI pro-high tier (retired); X-19/X-41/X-47 re-run on flash-high — disclosed |

## Augmentation lift

| Source set | Model | unaug % | aug % | Δ pts | comparability |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 68.0 | 88.0 | +20.0 | ok |
| GroundZero-50 | ChatGPT | 64.0 | 80.0 | +16.0 | ok |
| GroundZero-50 | DeepSeek | 68.0 | 84.0 | +16.0 | ok |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 40.0 | 96.0 | +56.0 | harness differs |
| Hard50 (new) | Claude | 72.0 | 96.0 | +24.0 | ok |
| Hard50 (new) | ChatGPT | 72.0 | 94.0 | +22.0 | harness differs |
| Hard50 (new) | DeepSeek | 68.0 | 92.0 | +24.0 | harness differs |
| Hard50 (new) | Gemini (unaug: web; aug: CLI) | 64.0 | 100.0 | +36.0 | harness differs |
| ExtraHard-50 | Claude | 66.0 | 100.0 | +34.0 | ok |
| ExtraHard-50 | ChatGPT | 66.0 | 100.0 | +34.0 | ok |
| ExtraHard-50 | DeepSeek | 56.0 | 94.0 | +38.0 | ok |
| ExtraHard-50 | Gemini (unaug: web; aug: CLI) | 54.0 | 100.0 | +46.0 | Gemini CLI pro-high tier; harness differs |

## Declared exclusions

| Source set | n | n cells | excluded | reason |
|---|---|---|---|---|
| GroundZero-50 | 50 | 50 | — | — |
| Hard50 (new) | 50 | 50 | — | — |
| ExtraHard-50 | 50 | 50 | — | — |

## Provenance ledger

| cell | passed/n | % | harness | source | sha | gold | notes |
|---|---|---|---|---|---|---|---|
| extrahard / ChatGPT / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json | bb19fd79144a | post-fix | count-only |
| extrahard / ChatGPT / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Claude / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 03c5689f29cd | post-fix | count-only |
| extrahard / Claude / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / DeepSeek / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 03c5689f29cd | post-fix | count-only |
| extrahard / DeepSeek / unaug | 28/50 | 56.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Gemini / aug | 50/50 | 100.0 | agy-cli | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json | 6043b0793db0 | post-fix | Gemini CLI pro-high tier (retired); X-19/X-41/X-47 re-run on flash-high — disclosed |
| extrahard / Gemini / unaug | 27/50 | 54.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| extrahard / Qwen / unaug | 21/50 | 42.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | 5455c165d6ef | post-fix | count-only |
| groundzero / ChatGPT / aug | 40/50 | 80.0 | mcp | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json | 26bf68e59aa6 | post-fix |  |
| groundzero / ChatGPT / unaug | 32/50 | 64.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix |  |
| groundzero / Claude / aug | 44/50 | 88.0 | mcp | benchmark/source_sets/groundzero/claude_augmented_grade_results.json | a0bc7ccb7625 | post-fix |  |
| groundzero / Claude / unaug | 34/50 | 68.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix |  |
| groundzero / DeepSeek / aug | 42/50 | 84.0 | mcp | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json | 42ac653966c4 | post-fix | SI-only |
| groundzero / DeepSeek / unaug | 34/50 | 68.0 | mcp | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json | 2e8f8288e360 | post-fix |  |
| groundzero / Gemini / aug | 48/50 | 96.0 | agy-cli | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json | 5a753572989a | post-fix |  |
| groundzero / Gemini / unaug | 20/50 | 40.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | af7478178bb7 | post-fix |  |
| groundzero / Qwen / unaug | 32/50 | 64.0 | mcp | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json | f8b8615eb935 | post-fix |  |
| hard50 / ChatGPT / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json | 33cd934a4128 | post-fix |  |
| hard50 / ChatGPT / unaug | 36/50 | 72.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Claude / aug | 48/50 | 96.0 | mcp | benchmark/source_sets/hard50/claude_augmented_grade_results.json | ecbb39aa68e1 | post-fix |  |
| hard50 / Claude / unaug | 36/50 | 72.0 | mcp | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json | 03ebde500fa9 | post-fix |  |
| hard50 / DeepSeek / aug | 46/50 | 92.0 | mcp | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json | 2ae21eb719d3 | post-fix | SI-only |
| hard50 / DeepSeek / unaug | 34/50 | 68.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Gemini / aug | 50/50 | 100.0 | agy-cli | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json | bedea6408f77 | post-fix |  |
| hard50 / Gemini / unaug | 32/50 | 64.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
| hard50 / Qwen / unaug | 33/50 | 66.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | c39dc851fee0 | not field-verifiable | count-only |
