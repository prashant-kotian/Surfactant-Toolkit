# Grading ledger, per source set (generated - do not hand-edit; rerun benchmark/build_results.py). Headline table: results/TABLE200.md

## unaug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| extrahard | ChatGPT | 33/50 | 66.0 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Claude | 33/50 | 66.0 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | DeepSeek | 28/50 | 56.0 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Gemini | 27/50 | 54.0 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Qwen | 21/50 | 42.0 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| groundzero | ChatGPT | 32/50 | 64.0 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | Claude | 35/50 | 70.0 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | DeepSeek | 34/50 | 68.0 | post-fix |  | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json |
| groundzero | Gemini | 21/50 | 42.0 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | Qwen | 32/50 | 64.0 | post-fix |  | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json |
| hard50 | ChatGPT | 36/50 | 72.0 | not field-verifiable | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Claude | 36/50 | 72.0 | post-fix |  | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json |
| hard50 | DeepSeek | 34/50 | 68.0 | not field-verifiable | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Gemini | 32/50 | 64.0 | not field-verifiable | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Qwen | 33/50 | 66.0 | not field-verifiable | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |

## aug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| extrahard | ChatGPT | 50/50 | 100.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json |
| extrahard | Claude | 50/50 | 100.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json |
| extrahard | DeepSeek | 47/50 | 94.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json |
| extrahard | Gemini | 50/50 | 100.0 | post-fix |  | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json |
| groundzero | ChatGPT | 40/50 | 80.0 | post-fix |  | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json |
| groundzero | Claude | 45/50 | 90.0 | post-fix |  | benchmark/source_sets/groundzero/claude_augmented_grade_results.json |
| groundzero | DeepSeek | 42/50 | 84.0 | post-fix |  | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json |
| groundzero | Gemini | 47/50 | 94.0 | post-fix |  | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json |
| hard50 | ChatGPT | 47/50 | 94.0 | post-fix |  | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json |
| hard50 | Claude | 48/50 | 96.0 | post-fix |  | benchmark/source_sets/hard50/claude_augmented_grade_results.json |
| hard50 | DeepSeek | 46/50 | 92.0 | post-fix |  | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json |
| hard50 | Gemini | 50/50 | 100.0 | post-fix |  | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json |

## gaps

- groundzero Gemini unaug: repo transcripts for the Sept-21 Gemini direct run: 11 no-content records (ungradeable, excluded); long-prose fields resolved via the grader's JUDGMENT_RESOLUTIONS table (none swung the total, n=33). The 200Q report's 22/50 came from a fuller backend grading pass; transcript confirmed absent on benchmark machine (RUNS_REPORT 2026-10-18 search).
