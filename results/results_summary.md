# RESULTS_MASTER (generated - do not hand-edit; rerun benchmark/build_results.py)

## unaug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| GZ50 | ChatGPT | 22/43 | 51.2 | post-fix |  | benchmark/groundzero_50/unaugmented_2026-09-21/grade_results.json |
| GZ50 | Claude | 5/10 | 50.0 | post-fix |  | benchmark/groundzero_50/unaugmented_2026-09-21/grade_results.json |
| GZ50 | DeepSeek | 29/42 | 69.0 | post-fix |  | benchmark/groundzero_50/deepseek_unaugmented_grade_results.json |
| GZ50 | Gemini-CLI | 42/43 | 97.7 | post-fix |  | benchmark/groundzero_50/antigravity_unaugmented_grade_results.json |
| GZ50 | Gemini | 10/33 | 30.3 | post-fix |  | benchmark/groundzero_50/unaugmented_2026-09-21/grade_results.json |
| GZ50 | Qwen | 27/43 | 62.8 | post-fix |  | benchmark/groundzero_50/qwen_unaugmented_grade_results.json |
| H50new | ChatGPT | 34/47 | 72.3 | post-fix | yes | benchmark/hard50/unaugmented_grade_stats.json |
| H50new | Claude | 33/47 | 70.2 | post-fix |  | benchmark/hard50/claude_unaugmented_grade_results.json |
| H50new | DeepSeek | 31/47 | 66.0 | post-fix | yes | benchmark/hard50/unaugmented_grade_stats.json |
| H50new | Gemini-CLI | 46/47 | 97.9 | post-fix |  | benchmark/hard50/antigravity_unaugmented_grade_results.json |
| H50new | Gemini | 29/47 | 61.7 | post-fix | yes | benchmark/hard50/unaugmented_grade_stats.json |
| H50new | Qwen | 30/47 | 63.8 | post-fix | yes | benchmark/hard50/unaugmented_grade_stats.json |
| X50 | ChatGPT | 32/49 | 65.3 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| X50 | Claude | 33/49 | 67.3 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| X50 | DeepSeek | 28/49 | 57.1 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| X50 | Gemini-CLI | 25/49 | 51.0 | post-fix |  | benchmark/extrahard50/antigravity_unaugmented_grade_results.json |
| X50 | Gemini | 27/49 | 55.1 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |
| X50 | Qwen | 21/49 | 42.9 | post-fix | yes | benchmark/transcripts/extrahard50_unaugmented/grade_results.json |

## aug

| bank | model | passed/n | % | gold | count-only | source |
|---|---|---|---|---|---|---|
| GZ50 | ChatGPT | 31/43 | 72.1 | post-fix |  | benchmark/groundzero_50/chatgpt_augmented_grade_results.json |
| GZ50 | Claude | 36/43 | 83.7 | post-fix |  | benchmark/groundzero_50/claude_augmented_grade_results.json |
| GZ50 | DeepSeek | 32/40 | 80.0 | post-fix |  | benchmark/groundzero_50/deepseek_augmented_grade_results.json |
| GZ50 | Gemini-CLI | 41/43 | 95.3 | post-fix |  | benchmark/groundzero_50/antigravity_augmented_grade_results.json |
| H50new | ChatGPT | 44/47 | 93.6 | post-fix |  | benchmark/hard50/chatgpt_augmented_grade_results.json |
| H50new | Claude | 45/47 | 95.7 | post-fix |  | benchmark/hard50/claude_augmented_grade_results.json |
| H50new | DeepSeek | 44/47 | 93.6 | post-fix |  | benchmark/hard50/deepseek_augmented_grade_results.json |
| H50new | Gemini-CLI | 47/47 | 100.0 | post-fix |  | benchmark/hard50/antigravity_augmented_grade_results.json |
| X50 | ChatGPT | 50/50 | 100.0 | post-fix | yes | benchmark/extrahard50/augmented_runs_2026-10-12/grade_results_chatgpt.json |
| X50 | Claude | 50/50 | 100.0 | post-fix | yes | benchmark/extrahard50/augmented_runs_2026-10-12/grade_results_augmented.json |
| X50 | DeepSeek | 48/50 | 96.0 | post-fix | yes | benchmark/extrahard50/augmented_runs_2026-10-12/grade_results_augmented.json |
| X50 | Gemini-CLI | 49/49 | 100.0 | post-fix |  | benchmark/extrahard50/antigravity_augmented_grade_results.json |

## gaps

- GZ50 Gemini unaug: repo transcripts for the Sept-21 Gemini direct run: 11 no-content records (ungradeable, excluded) + 8 NEEDS_JUDGMENT scored as conservative FAIL inside n=33. The 200Q report's 22/50 came from a fuller backend grading pass; transcript confirmed absent on benchmark machine (RUNS_REPORT 2026-10-18 search).
