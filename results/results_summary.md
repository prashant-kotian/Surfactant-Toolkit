# Grading ledger, per source set (generated - do not hand-edit; rerun benchmark/build_results.py). Headline table: results/TABLE200.md

## unaug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| extrahard | ChatGPT | 32/49 | 65.3 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Claude | 33/49 | 67.3 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | DeepSeek | 28/49 | 57.1 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Gemini-CLI | 25/49 | 51.0 | post-fix |  | benchmark/source_sets/extrahard/antigravity_unaugmented_grade_results.json |
| extrahard | Gemini | 27/49 | 55.1 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| extrahard | Qwen | 21/49 | 42.9 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| groundzero | ChatGPT | 22/43 | 51.2 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | Claude | 5/10 | 50.0 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | DeepSeek | 29/42 | 69.0 | post-fix |  | benchmark/source_sets/groundzero/deepseek_unaugmented_grade_results.json |
| groundzero | Gemini-CLI | 42/43 | 97.7 | post-fix |  | benchmark/source_sets/groundzero/antigravity_unaugmented_grade_results.json |
| groundzero | Gemini | 10/33 | 30.3 | post-fix |  | benchmark/source_sets/groundzero/unaugmented_2026-09-21/grade_results.json |
| groundzero | Qwen | 27/43 | 62.8 | post-fix |  | benchmark/source_sets/groundzero/qwen_unaugmented_grade_results.json |
| hard50 | ChatGPT | 34/47 | 72.3 | post-fix | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Claude | 33/47 | 70.2 | post-fix |  | benchmark/source_sets/hard50/claude_unaugmented_grade_results.json |
| hard50 | DeepSeek | 31/47 | 66.0 | post-fix | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Gemini-CLI | 46/47 | 97.9 | post-fix |  | benchmark/source_sets/hard50/antigravity_unaugmented_grade_results.json |
| hard50 | Gemini | 29/47 | 61.7 | post-fix | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |
| hard50 | Qwen | 30/47 | 63.8 | post-fix | yes | benchmark/source_sets/hard50/unaugmented_grade_stats.json |

## aug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| extrahard | ChatGPT | 50/50 | 100.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json |
| extrahard | Claude | 50/50 | 100.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json |
| extrahard | DeepSeek | 48/50 | 96.0 | post-fix | yes | benchmark/source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json |
| extrahard | Gemini-CLI | 49/49 | 100.0 | post-fix |  | benchmark/source_sets/extrahard/antigravity_augmented_grade_results.json |
| groundzero | ChatGPT | 31/43 | 72.1 | post-fix |  | benchmark/source_sets/groundzero/chatgpt_augmented_grade_results.json |
| groundzero | Claude | 36/43 | 83.7 | post-fix |  | benchmark/source_sets/groundzero/claude_augmented_grade_results.json |
| groundzero | DeepSeek | 32/40 | 80.0 | post-fix |  | benchmark/source_sets/groundzero/deepseek_augmented_grade_results.json |
| groundzero | Gemini-CLI | 41/43 | 95.3 | post-fix |  | benchmark/source_sets/groundzero/antigravity_augmented_grade_results.json |
| hard50 | ChatGPT | 44/47 | 93.6 | post-fix |  | benchmark/source_sets/hard50/chatgpt_augmented_grade_results.json |
| hard50 | Claude | 45/47 | 95.7 | post-fix |  | benchmark/source_sets/hard50/claude_augmented_grade_results.json |
| hard50 | DeepSeek | 44/47 | 93.6 | post-fix |  | benchmark/source_sets/hard50/deepseek_augmented_grade_results.json |
| hard50 | Gemini-CLI | 47/47 | 100.0 | post-fix |  | benchmark/source_sets/hard50/antigravity_augmented_grade_results.json |

## gaps

- groundzero Gemini unaug: repo transcripts for the Sept-21 Gemini direct run: 11 no-content records (ungradeable, excluded) + 8 NEEDS_JUDGMENT scored as conservative FAIL inside n=33. The 200Q report's 22/50 came from a fuller backend grading pass; transcript confirmed absent on benchmark machine (RUNS_REPORT 2026-10-18 search).
