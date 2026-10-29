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
the first-release harness; their raw transcripts are in
`benchmark/transcripts/orig_hard50_first_release/`.

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

Four Gemini web fill answers were given two to a chat (GZ-13 with GZ-17, GZ-33 with GZ-50; the paired
questions are on unrelated topics and share no data). GZ-17 is kept: asked five more times, each in a
new chat (`benchmark/reproducibility/sb_extension/`), Gemini answered it correctly four times. GZ-13,
GZ-33 and GZ-50 were asked again on 2026-10-26, each in a new chat with the same settings. GZ-50 answered and its new
answer replaces the shared-chat one (PASS -> PASS; the old answer is kept in
`answers_fill/gemini_unaug/superseded_shared_chat/`). GZ-13 (three attempts) and GZ-33 (three
attempts) returned no answer, only a Gemini server error, so their shared-chat answers are kept;
each `.meta.json` records this. No verdict changes.

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

GZ-17 / SB-017 (2026-10-23): the gold CMC was 8.185 mM, the raw piecewise
break-point (`OrchestratorResult.cmc.cmc_mM`), which sits at the geometric mean
of the 4.09/16.37 mM gap because the three-segment fit accepted a two-point,
rising "plateau". The toolkit's authoritative value (`OrchestratorResult.cmc_mM`,
the value the MCP tool returns) flags that plateau and refines it to 16.37 mM;
the development repository had made this correction on 2026-09-23, and it was
lost when the release gold was assembled. Gold restored to 16.37 mM (Gamma_max
and A_min unchanged); `tests/test_benchmark_golds.py` now ties the orchestrated
GroundZero CMC golds to the authoritative value. Re-graded: Claude unaugmented
and augmented FAIL->PASS, Gemini web unaugmented FAIL->PASS, Gemini CLI
augmented (reported 8.19) PASS->FAIL; every other arm still fails GZ-17 on
another field. Headline: Claude 65.0% / 96.5%, Gemini 51.5% / 98.5%. The
groundzero cells for Claude and ChatGPT augmented and Qwen and DeepSeek
unaugmented are read from frozen grading files that no grader in this
repository writes (`grade_chat_runs.py` grades the same transcripts into
`grade_results_session_*.json`); their GZ-17 records were replaced by the
re-graded ones after checking that the model answers are identical.

Toolkit corrections found by the Supporting Information equation checks (2026-10-25):
- `dynamic_tension.dynamic_surface_tension_forward` / `_ode` converted the surface excess with
  gamma0 - nRT*Gamma although adsorption follows the Langmuir isotherm, whose equation of state is
  the Szyszkowski form gamma0 + nRT*Gamma_max*ln(1 - Gamma/Gamma_max); the surface step also
  overshot under stiff kinetics at low concentration. Both corrected (exact Riccati surface step,
  Crank-Nicolson substep capped at r = 0.5); tests now tie the long-time limit to
  `szyszkowski_surface_tension` and the surface excess to the independent Ward-Tordai solver.
- SB-197 (orig Q47) asks for this tool's gamma at 0.1 s; its reference value moves from 64.441 to
  57.982 mN/m (`repair_2026-10-20/q47_regrade.py`, recomputed from the toolkit). Unaugmented arms
  are graded from their stored answers (all still FAIL). Augmented arms were re-run against the
  corrected tool (2026-10-25/26; answers and meta in `repair_2026-10-20/answers_q47/`): Claude,
  ChatGPT, DeepSeek and Gemini CLI all answered 57.982 (PASS). `q47_regrade_results.json` is
  written and applied; every SB-197 verdict is unchanged, so no reported cell moves.
- `kinetics.aniansson_wall_two_regime` placed the shared factor (1 + sigma^2 a/n) on the wrong
  side when solving for R; affects only c_total != CMC; no benchmark question or application uses it.
- `emulsion.simulate_emulsion_stability` counted a size bin as separated only after its full travel
  time H/v (all drops starting at the bottom); drops now start uniformly distributed, so a bin's
  creamed share grows as t/t_c and a monodisperse t50 is H/(2v). No benchmark question uses this tool;
  the application P4 reference value is recomputed (benchmark/applications/compute_p4.py).
- orig_hard50 reference values checked against the current toolkit (every augmented wording names
  one tool with its inputs; paper3 si/check_orig_golds.py executes them): 48 of 50 reproduce exactly.
  SB-197 is the dynamic-tension item above. SB-167 (Q17, hlb_davies_guo_ecl(12, 6)) now returns 9.62
  (Guo EO-unit weight 0.33 since 2026-10-20) against the reference 33.29. The two wordings define
  their answers differently: the augmented wording asks for the tool's output, the unaugmented
  wording states every weight of the calculation itself (ether oxygen 1.3), for which 33.29 is the
  answer. Reference values are therefore per wording: augmented 9.62 (release toolkit), unaugmented
  33.29 (recomputed from the stated calculation). Neither wording is changed after the runs.
  Augmented arms re-run against the current tool (`repair_2026-10-20/answers_q17/<arm>/`) and
  regraded by `repair_2026-10-20/q17_regrade.py` (-> `q17_regrade_results.json`, applied through
  `corrections.orig_regrade`) on 2026-10-26: Claude, ChatGPT, DeepSeek and Gemini CLI all answered
  9.62 (PASS); every SB-167 verdict is unchanged, so no reported cell moves.
- `orig_hard50/question_bank_unaugmented.txt` added: the unaugmented wording of SB-151..SB-200 used in
  the first-release runs (from the development repository, HARD60_UNAUGMENTED_QUESTIONS.txt).
- SDS alpha = 0.272 (J. Phys. Chem. B 105 (2001) 6798) is credited to B. L. Bales alone in
  `thermodynamics.py` (the returned `source` string) and in the groundzero gold source notes; it was
  given as a four-author paper. No value changes. The question banks keep their as-run reference notes.
- `hlb.DAVIES_HYDROPHILIC_GROUPS["SO4Na"]` 38.6 -> 38.7 (Davies 1957). The benchmark questions
  that use this group state 38.6 in their own text, so no reference value changes.

SB-166, SB-167 and SB-197: one target per question in both wordings (2026-10-28, external audit):
- `dynamic_tension.dynamic_surface_tension_forward` was not converged at its default settings: the
  uniform grid was sized from the latest output time only, so early times were under-resolved
  (SB-197: 57.982 at t = 0.1 s against a converged 57.933 mN/m; up to 0.6 mN/m in other cases).
  The solver now uses a stretched grid resolving the earliest output time, backward-Euler steps
  proportional to t with the Langmuir surface kinetics solved in the same implicit step, and
  Richardson extrapolation; defaults agree with a 4x finer grid and step to about 1e-3 mN/m.
  `n_cells` keeps its meaning (number of cells); `max_substep_s` is now optional and `rel_step`
  (default 0.01; 0.04 in the fit) sets the step. The fit wrapper passes both through.
- `adsorption.ward_tordai_adsorption` (Langmuir mode) held the last solved subsurface concentration
  over the interval being solved and assumed a clean surface up to the first output time; Gamma was
  biased by about 1 % at all times and up to 7 % early, independent of the grid. It now uses product
  integration (c_s piecewise linear on geometric internal nodes, exact kernel weights, closed-form
  Langmuir closure). It reproduces the exact linear-isotherm solution to 2e-4 and agrees with the
  finite-difference solver to 2e-4 under diffusion control. `n_sqrt` is now the number of internal
  nodes; `bisection_tol` is accepted and unused.
- SB-166 (Q16, Ward-Tordai Gamma at 1 s, c0 = 0.02 M): reference 2.8566 -> 2.8568 umol/m^2 (both
  solvers; the long-time asymptote Gamma_eq - (dGamma/dc_s) Gamma/sqrt(pi D t) gives 2.85679).
  Wording unchanged. Unaugmented answers regraded from the stored records (Claude and Gemini web
  answered 2.8568: FAIL -> PASS); augmented arms re-asked against the corrected tool, all 2.8568.
- SB-167 (Q17): the unaugmented wording now states the Davies group number of an EO unit (0.33)
  instead of the ether-oxygen number (1.3), so both wordings have the target 9.62. Unaugmented arms
  answered the new wording; the augmented answers (9.62, wording unchanged) stand.
- SB-197 (Q47): both wordings now state every setting (gamma0 72 mN/m, n = 1, T = 298.15 K,
  k_off = 0.5 1/s; K = 1500 L/mol named as the Langmuir equilibrium constant with k_a = K k_off)
  and ask for 2 decimals; target 57.93. All arms answered the new wording.
- Answers and metadata: `repair_2026-10-20/answers_2026-10-28/<arm>/` (prompts of the unaugmented
  wording in `prompts_unaug/`); graded by `repair_2026-10-20/orig_q16_q17_q47_grade.py` ->
  `orig_q16_q17_q47_results.json`, applied through `corrections.orig_regrade` ahead of the earlier
  q17/q47 regrade files. Verdict changes (all unaugmented): Claude SB-166 and SB-197, Gemini web
  SB-166 and SB-197, ChatGPT and DeepSeek SB-197 FAIL -> PASS. Unchanged: Gemini web SB-167 (9.64),
  Qwen SB-197 (57.75), every augmented verdict. Gemini web returned no answer to SB-197 in three
  fresh chats and answered in a later one; one capture that showed the previous chat was discarded.
- Application P3 (two solvers compared) recomputed: the solvers coincide at 5 s; the remaining early
  difference (up to 7 % in Gamma at 0.1 mM) is the finite adsorption rate (k_off = 200 1/s) of the
  finite-difference solver.
- `build_figure_data.py` now applies the orig_hard50 regrade layer to the grand totals, as
  `build_bench200.py` does (the earlier regrades changed no verdict, so the omission was silent).

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
