# Pending repeat-run prompts (batch stopped 2026-10-23 18:28 IST; 101 left)

Prompt text: the SB-xxx_*.md file in this folder, pasted exactly as is. SB-038 also needs SB-038_chart.png attached.

Rules for every prompt:
- New chat for each prompt (fresh session, no memory or history).
- No web search, no tools, no MCP or extensions.

## Settings per model
| Model | Setting |
|---|---|
| ChatGPT | GPT-5.6 Sol, thinking effort High |
| Claude | Opus 4.8, effort High, incognito chat |
| DeepSeek | DeepThink on, Search off |
| Gemini | Gemini 3.1 Pro, thinking High |
| Qwen | Qwen3.8-Max |

## Pending prompts
| Model | Round | Count | Questions |
|---|---|---|---|
| ChatGPT | r5 | 11 | all 11 |
| DeepSeek | r5 | 11 | all 11 |
| Claude | r3 | 1 | SB-061 |
| Claude | r4 | 10 | all except SB-017 |
| Claude | r5 | 11 | all 11 |
| Gemini | r1 | 3 | SB-017, SB-011, SB-012 |
| Gemini | r2 | 1 | SB-038 |
| Gemini | r3, r4, r5 | 33 | all 11, each round |
| Qwen | r4 | 9 | all except SB-011, SB-042 |
| Qwen | r5 | 11 | all 11 |

Total: ChatGPT 11, DeepSeek 11, Claude 22, Gemini 37, Qwen 20 = 101.

"All 11" = SB-017, SB-011, SB-012, SB-042, SB-061, SB-065, SB-069, SB-070, SB-073, SB-169, SB-038.

## Saving answers
Save the full reply text, including the FINAL_JSON block, to `manual_answers/r<N>/<model>_unaug/<QID>.txt`. QID is the file suffix: GZ17, GZ11, GZ12, GZ42, H11, H15, H19, H20, H23, Q19, GZ38.

Example: `manual_answers/r5/chatgpt_unaug/GZ17.txt`.

Also note the date, and any error or refusal the model gave.
