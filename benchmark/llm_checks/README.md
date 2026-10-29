# Language models on the five literature checks

Five published cases (LC1-LC5 = cases 1-5 of `../literature_checks`) given to two language models in
three stages, 30 chats in all (5 cases x 3 stages x 2 models, one run each). Paper: Section 6.6;
Supporting Information: Section S10.

| Stage | Prompt contains | Toolkit |
|---|---|---|
| S1 | the paper's printed data, its printed equations and the question the paper answers | off |
| S2 | S1 plus the paper's stated result and conclusion, quoted, and the question whether it follows | off |
| S3 | S2 word for word plus the line "The SurfactantKit tools are available to you." | on (MCP) |

Prompts contain no expected values, verdicts, method hints or tool names beyond what the quoted paper
prints; each ends by asking for a `FINAL_JSON` object. Models: ChatGPT (GPT-5.6 Sol, thinking effort
High; fresh chat, saved-memory and chat-history reference off) and Claude (Opus 4.8, effort High;
incognito chat). Web search off for both. The SurfactantKit MCP server was connected in S3 only.

## Files

- `prompts/LC<n>_S<k>.txt`: the 15 prompts, as sent.
- `run_manifest.tsv`: the 30 runs (run id, case, stage, model, prompt, mcp, web search).
- `truth/LC<n>.json`: reference values and scoring rules per key and stage, generated from
  `../literature_checks/case*_results.json` and fixed before any run.
- `answers/<run_id>/answer.txt`: the model's reply; `meta.json`: model and setting as shown, chat
  mode, toolkit calls (Claude: shown in the page; ChatGPT: the MCP server request log), number of
  ChatGPT built-in Python steps, timestamps.
- `score.py`: scores the last `FINAL_JSON` object of each answer against `truth/` and writes
  `summary.tsv`. Key kinds: `numeric` (value +/- tolerance), `interval` (inside the 95% Monte Carlo
  interval of the case result), `sign`, and string keys (yes/no/cannot_determine/no_solution).
  Parsing rules: curly double quotes in the JSON object read as plain quotes; JSON true/false on a
  yes/no key read as yes/no; a string verdict that starts with yes, no or cannot_determine is scored
  on that word.
- `test_score.py`: self-test with correct, paper-repeating and malformed fixtures for every case and
  stage and one fixture per parsing rule.

    python score.py
    python -X utf8 test_score.py
