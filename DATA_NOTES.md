# Data notes

## Benchmark composition

SurfBench-200 is one 200-question benchmark (ids `SB-001..SB-200`). It was
authored as four 50-question working groups — groundzero (definitional /
standard practice), hard50 (multi-step computation), extrahard (adversarial),
and orig_hard50 (toolkit-execution items carried in from the project's first
benchmark release) — and merged by `benchmark/build_question_set.py`. The
group names survive only as `source_sets/` directories and `authoring_alias`
columns for provenance; all reporting uses the merged ids and the question
classification (`skill type` and `subconcept` per question) it carries.
Question texts are frozen; versions are tracked through the gold files, not
silent edits. The orig_hard50 grading JSONs are committed frozen inputs from
the first-release harness; their raw transcripts are not re-committed here.

## Scoring

Field-level comparison against `gold_answers.json` (numeric tolerance by
magnitude class; categorical labels by token-set equality after stripping
justifications; NO-JSON detection per transcript). The rule set lives in the
per-source-set graders and is exercised by their self-tests.

## Exclusions and 2026-10-20/21 corrections

**Every model arm is scored on all 200 questions (2026-10-21).** Until 2026-10-20 the chart/image
questions GZ-15, GZ-36..40 and H-36..38, the retired GZ-26 and X-19 (authoring typo until
2026-10-16) were excluded for every arm (`EXCLUSIONS` in `benchmark/build_results.py`;
`formerly_excluded` in `benchmark/question_set.json`). They were re-run and are now scored. A
corrections layer (`benchmark/repair_2026-10-20/corrections.py`) applies the re-runs and re-grades
to each cell, and every change it makes is logged in that cell's `corrections` list in `results/results.json`:

- Chart/typo repair. The nine chart questions and X-19 were re-run on 2026-10-20/21 with the chart
  attached and the corrected X-19 text, under the original model settings, with the condition checked
  on the live page before each send (`answers/*/*.meta.json`). DeepSeek augmented and Gemini
  augmented (Gemini CLI) originally received the charts and keep their original chart verdicts. A question stays excluded
  for an arm only when that arm has no gradable re-run answer (reason logged).
- Counterion-factor correction. The pseudo-phase ionic factor is (1 + beta) = (2 - alpha); earlier
  toolkit text and 13 question texts (GZ-16, GZ-20, H-02, H-03, H-08, H-16, H-19, H-31, H-38, H-47,
  H-49, X-41, X-47) stated the inverse. Texts and golds were corrected (H-46 gold, H-04 label as
  well) and every arm was re-run on the corrected texts; those verdicts replace the originals. An arm
  without a gradable re-run answer is excluded on that question, never scored on the stale text.
  H-46 (gold-only error) is regraded from stored answers where a record holds them (Claude
  unaugmented) and re-run otherwise.
- GZ-16/GZ-20 electrolyte condition. These texts did not state the electrolyte condition. Arms that
  assumed no added salt keep their answers; ChatGPT augmented, which withheld Gamma_max/A_min for
  that reason, was re-run with the sentence "The solution contains no added salt; the surfactant is
  the only electrolyte present." (author decision; recorded in each answer's `question_variant`).
- Choice-field regrade (groundzero). Golds of the form "<choice> (<justification>)" (GZ-01/GZ-20
  isotherm_model, GZ-29 selected_model, GZ-23/GZ-35 valid_hlb_method) were scored by token overlap,
  which failed correct one-word answers and passed opposite choices. Rule for every arm: the field
  passes iff the answer's first word is the gold's choice ("Szyszkowski", the Langmuir equation of
  state, counts as Langmuir). A record is changed only if its stored verdict is reproducible from its
  own field flags.
- GZ-50 Delta G fields: the question names no unit; a value given in J/mol is converted to kJ/mol before the
  standard comparison (author decision 2026-10-21; affects Gemini web, now PASS).
- GZ-15 estimated CMC is graded at +/-10% (author decision; the gold note asked for a wider, unquantified
  tolerance). Free-text GZ-39 fields carry one written judgment per arm in `grade_repair.py`.

Other standing rules: orig_hard50 NOT_ANSWERED records (Claude unaugmented 1, ChatGPT unaugmented 4)
are excluded from n as ungradable; the 2026-09-21 Gemini unaugmented groundzero run has 11 no-content
records removed from n, and an earlier alternative grading pass (22/50) is not reproducible from this
repository and is not reported.

## Web search in the unaugmented condition (must be stated in the paper)

In the unaugmented condition Claude was run with web search disabled and answered every question from its own reasoning, whereas the other models (ChatGPT, DeepSeek, Gemini web, Qwen) used web search on some questions; the unaugmented comparison is therefore, if anything, biased against Claude. Example: the 2026-10-21 orig_hard50 re-runs (ChatGPT Q24/Q25/Q29 cite web sources; Claude Q07 used none; recorded per answer as `web_search_used` in `answers_fill/*/*.meta.json`).

## Gemini arms

Gemini is reported as one model with two clients: the unaugmented arm is Gemini web (gemini.google.com,
no tools), the augmented arm is Gemini CLI (Antigravity) with the surfactant-computation MCP server.
Gemini CLI unaugmented runs were made but are not part of the results: the CLI keeps its own coding
and web-search tools, so they were not a true unaugmented baseline. Those runs are kept in the repository for
reference (`source_sets/*/antigravity_unaug*`, `groundzero/gemini_cli_2026-09-19/`,
`repair_2026-10-20/answers_dg/agy_unaug/`) but are not used in any comparison (2026-10-21). The Gemini CLI extrahard tier used originally (pro-high) was retired; X-19, X-41 and X-47
were re-run on gemini-3.8-flash-high, the model used for all Gemini CLI re-runs.

## Net-200 fill (2026-10-21)

GZ-26 (retired 2026-10-18) was re-run on its current text for every arm and graded against the
corrected gold; answers missing from the original runs (Gemini web GZ-01/05/11/13/17/27/29/33/50,
DeepSeek GZ-04/17/49, orig_hard50 Claude Q07 and ChatGPT Q24/25/29/44) were re-run the same way
(`benchmark/repair_2026-10-20/answers_fill/`, graded by `grade_repair.py --fill`; orig_hard50 items by
the first-release rule in `orig_hard50_grade.py`, which reproduces all 250 original verdicts). The
orig_hard50 re-runs use the tool-free question text with a one-line FINAL_JSON answer format.

## Classification

Labels are rule-derived from gold field names, declared failure axes, and
question text (see `build_question_set.py`), with 28 manual corrections in
`question_set_overrides.py` after full review of the generated assignment.
Small categories (n ≤ 3: rheology, solubilization, data_fitting,
analysis_synthesis, foaming_films, micelle_thermodynamics, surface_wetting)
are reported but should be read as descriptive only.

## Corrections

H-04 (reference-temperature mismatch), H-26 (surfactant identity), GZ-48
(factor-of-2 row mix-up): golds corrected on 2026-10-18 with per-field
`remediation_note`. The GZ-48 correction revealed a masked model error
(in the since set-aside Gemini CLI unaugmented runs, PASS->FAIL): reported rather than smoothed. See
`results/REPRODUCIBILITY.md`.

Re-grade completeness (2026-10-19): the DeepSeek augmented groundzero and
hard50 cells had NOT been re-graded after the gold corrections (their graders
sat outside the documented chain; the hard50 one also had a stale import and
answer path). Re-graded: GZ-48 FAIL->PASS, H-04 FAIL->PASS, H-26 PASS->FAIL
(DeepSeek's 42-carbon answer matched the superseded gold). `gold_state` is now
verified from the gold values stored in each cell's records rather than
inferred from file dates; every per-field cell verifies as post-fix, and the
count-only hard50 unaugmented cells are labelled "not field-verifiable". All
graders that produce a reported cell now run in REPRODUCE.md section 1, and
Claude's hard50 unaugmented transcript is committed in-repo
(`benchmark/source_sets/hard50/claude_unaugmented/`).

Claude's groundzero-unaugmented row (2026-10-19): the 10-question partial used
through 2026-10-18 was a different, genuinely truncated Sept-2026 capture; a
complete GZ-01..50 transcript was subsequently located and verified (all 50
answer blocks structurally intact, independently-derived rather than copied
from any gold version). Grading it surfaced that `grade_unaugmented.py` had no
mechanism to resolve long-prose judgment fields -- `results/results.json`'s
verdict conversion silently scored every unresolved field as a FAIL, which
also understated ChatGPT's and Gemini's groundzero-unaugmented rows (Gemini's
score did not move; no judgment field happened to swing its total). A
`JUDGMENT_RESOLUTIONS` table (31 entries, one rationale per field, same
discipline as the Antigravity graders) now resolves every pending field.
