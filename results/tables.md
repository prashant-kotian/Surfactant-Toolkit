# Per-source-set results (internal diagnostic view)

> These are the authoring groups behind the single SurfBench-200 benchmark;
> the headline table for reporting is `results/TABLE200.md`. (generated — rerun
> benchmark/build_tables.py)

## Unaugmented

| Source set | Model | passed | n | % | flags |
|---|---|---|---|---|---|
| GroundZero-50 | Claude | 35 | 50 | 70.0 |  |
| GroundZero-50 | ChatGPT | 32 | 50 | 64.0 |  |
| GroundZero-50 | DeepSeek | 34 | 50 | 68.0 |  |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 21 | 50 | 42.0 |  |
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
| GroundZero-50 | Claude | 45 | 50 | 90.0 |  |
| GroundZero-50 | ChatGPT | 40 | 50 | 80.0 |  |
| GroundZero-50 | DeepSeek | 42 | 50 | 84.0 | SI-only |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 47 | 50 | 94.0 |  |
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
| GroundZero-50 | Claude | 70.0 | 90.0 | +20.0 | ok |
| GroundZero-50 | ChatGPT | 64.0 | 80.0 | +16.0 | ok |
| GroundZero-50 | DeepSeek | 68.0 | 84.0 | +16.0 | ok |
| GroundZero-50 | Gemini (unaug: web; aug: CLI) | 42.0 | 94.0 | +52.0 | harness differs |
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
| extrahard / ChatGPT / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json | 96d767dd403d | post-fix | count-only |
| extrahard / ChatGPT / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | ab4f2dc6cf86 | post-fix | count-only |
| extrahard / Claude / aug | 50/50 | 100.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 95b75a6c6770 | post-fix | count-only |
| extrahard / Claude / unaug | 33/50 | 66.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | ab4f2dc6cf86 | post-fix | count-only |
| extrahard / DeepSeek / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json | 95b75a6c6770 | post-fix | count-only |
| extrahard / DeepSeek / unaug | 28/50 | 56.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | ab4f2dc6cf86 | post-fix | count-only |
| extrahard / Gemini / aug | 50/50 | 100.0 | agy-cli | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json | cd584b05a3bf | post-fix | Gemini CLI pro-high tier (retired); X-19/X-41/X-47 re-run on flash-high — disclosed |
| extrahard / Gemini / unaug | 27/50 | 54.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | ab4f2dc6cf86 | post-fix | count-only |
| extrahard / Qwen / unaug | 21/50 | 42.0 | mcp | benchmark/transcripts/extrahard50_unaugmented/grade_results.json | ab4f2dc6cf86 | post-fix | count-only |
| groundzero / ChatGPT / aug | 40/50 | 80.0 | mcp | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json | a30e2044075e | post-fix |  |
| groundzero / ChatGPT / unaug | 32/50 | 64.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | fac4e1fa97cc | post-fix |  |
| groundzero / Claude / aug | 45/50 | 90.0 | mcp | benchmark/source_sets/groundzero/claude_augmented_grade_results.json | 60d3e6534a59 | post-fix |  |
| groundzero / Claude / unaug | 35/50 | 70.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | fac4e1fa97cc | post-fix |  |
| groundzero / DeepSeek / aug | 42/50 | 84.0 | mcp | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json | ffa06b9383e3 | post-fix | SI-only |
| groundzero / DeepSeek / unaug | 34/50 | 68.0 | mcp | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json | 09ca70854138 | post-fix |  |
| groundzero / Gemini / aug | 47/50 | 94.0 | agy-cli | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json | 6b0bb113656d | post-fix |  |
| groundzero / Gemini / unaug | 21/50 | 42.0 | mcp | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json | fac4e1fa97cc | post-fix |  |
| groundzero / Qwen / unaug | 32/50 | 64.0 | mcp | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json | d9dc2a480489 | post-fix |  |
| hard50 / ChatGPT / aug | 47/50 | 94.0 | mcp | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json | 9a672cc3b8a5 | post-fix |  |
| hard50 / ChatGPT / unaug | 36/50 | 72.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 94eea68cede4 | not field-verifiable | count-only |
| hard50 / Claude / aug | 48/50 | 96.0 | mcp | benchmark/source_sets/hard50/claude_augmented_grade_results.json | 7c1bc1c1d001 | post-fix |  |
| hard50 / Claude / unaug | 36/50 | 72.0 | mcp | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json | 2c4aba3d3806 | post-fix |  |
| hard50 / DeepSeek / aug | 46/50 | 92.0 | mcp | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json | ced165808045 | post-fix | SI-only |
| hard50 / DeepSeek / unaug | 34/50 | 68.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 94eea68cede4 | not field-verifiable | count-only |
| hard50 / Gemini / aug | 50/50 | 100.0 | agy-cli | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json | 27c80d49afea | post-fix |  |
| hard50 / Gemini / unaug | 32/50 | 64.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 94eea68cede4 | not field-verifiable | count-only |
| hard50 / Qwen / unaug | 33/50 | 66.0 | cli | benchmark/source_sets/hard50/unaugmented_grade_stats.json | 94eea68cede4 | not field-verifiable | count-only |
