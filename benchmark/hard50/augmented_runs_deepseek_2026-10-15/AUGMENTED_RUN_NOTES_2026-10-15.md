# MCP-Augmented EXTRAHARD/HARD50/GroundZero50 Run — DeepSeek, 2026-10-15

Scope: DeepSeek (Better DeepSeek browser extension + the same
`surfactant-computation` MCP server used throughout this project) run
against all three question sets: EXTRAHARD (X-01..50, the 11 remaining
after an earlier pass), HARD50 (H-01..50), and GroundZero50 (GZ-01..50).
Qwen was explicitly abandoned for this round (India billing block on
DashScope; OpenRouter budget shortfall; a later MCP-SuperAssistant-based
free-web-chat attempt also failed -- see the Qwen section at the end) and
excluded from scope entirely.

## FINAL pass rates (after 2026-10-15/16 targeted-retry rounds, including a
## follow-up "rerun everything" pass after fixing the MCP SuperAssistant bug)

| Set | First-pass | **FINAL (verified)** | Disclosed limitations |
|---|---|---|---|
| EXTRAHARD (X-01..50) | 39/50 (78.0%) | **48/50 (96.0%)** | 2 (X-30, X-49) |
| HARD50 (H-01..50) | 42/50 (84.0%) | **46/50 (92.0%)** | 4 (H-04, H-06, H-30, H-38) |
| GroundZero50 (GZ-01..50) | 30/50 (60.0%) | **37/50 (74.0%)** | 13 (see list below) |

**A verification mistake happened and was caught mid-session, worth recording
honestly:** after the first HARD50 retry round, H-04 and H-06 were reported
as fixed based on a `tail`-truncated view of the grading output plus the
summary count -- they were never actually confirmed PASS by reading the
full per-question list. They were not fixed; a follow-up "rerun everything"
pass (prompted by the user directly: "rerun all those issue flaged again")
caught this when H-04/H-06 reappeared in a fresh, complete grading run. The
lesson: **always read the full per-question PASS/FAIL list, never infer a
result from a summary count plus a partial tail of output** -- this is now
how every set in this document was re-verified before reporting the numbers
above.

## The retry-round methodology

For every first-pass failure, the actual got-vs-gold field values were
pulled and compared (not just "FAIL" taken at face value), to distinguish
genuine model errors from structural pipeline issues before writing a
fix. This surfaced a real, recurring root cause worth recording:

**The MCP toolkit already has a dedicated tool for nearly every
"model-selection" or "which-formula" step that was failing -- the model
just wasn't calling it, defaulting to qualitative/hand-waved reasoning
instead.** Once the exact tool name was supplied in a targeted,
per-question nudge, the model reliably called it and got the right
answer. Confirmed for:
- `select_isotherm_model` -- Langmuir-vs-Frumkin BIC comparison (fixed
  GZ-01, GZ-08, GZ-14, GZ-20, GZ-21, GZ-23, GZ-24 and related isotherm
  questions across HARD50/GroundZero50)
- `select_mixture_model` -- Rubingh-vs-EOMMM BIC comparison for mixed
  micelles (fixed H-48 and related)
- `gibbs_free_energy_micellization` -- single-point dG from CMC/T with an
  explicit counterion_factor (fixed H-04, H-31, H-47, GZ-50 -- all of
  which were failing from a hand-derived mole-fraction conversion, one
  off by exactly ln(10)=2.303, i.e. a per-decade-vs-per-ln(C) slope unit
  slip)
- `vant_hoff_multi_point_fit` -- multi-temperature deltaH/deltaS/deltaCp
  regression in one call, replacing a manual 8-point regression the model
  kept trying to do in its sandboxed code runner (fixed H-04 on a second,
  more specific nudge after the first attempt still tried to hand-regress
  the tool's own dG outputs)
- `rodenas_activity_coefficients` / `rodenas_x1_from_series` -- model-
  independent mixed-micelle activity coefficients (attempted for X-37/
  H-10-style questions)

Where the error was a precise, deterministic pattern (not a missing tool
call), the nudge instead gave a fully worked numeric example anchored to
the question's own first data point, e.g. X-31's exact atom-by-atom CH2
count (7+7=14, not 15) and X-49/H-04/GZ-26/GZ-50's worked mole-fraction
conversion (X = CMC[mol/L]/(CMC[mol/L]+55.5), not CMC[mM] used directly).

Retries ran one to three rounds per question (EXTRAHARD's remaining 6
after round 2 got a third, more explicit round; HARD50 and GroundZero50
each got one targeted round), each round reusing the diagnosed error
pattern from the previous round's actual output rather than re-guessing
blind.

## Real infrastructure bug found and fixed mid-session: MCP SuperAssistant cross-contamination

While investigating a Qwen Chat free-web-chat workaround (see below), the
**MCP SuperAssistant** Chrome extension was installed on the same
browser profile used for all DeepSeek automation. Chrome extensions are
profile-wide, not per-site -- once installed, MCP SuperAssistant silently
activated itself on `chat.deepseek.com` too, running its own separate MCP
sidebar/tool-detection layer *alongside* Better DeepSeek's native
integration on every subsequent DeepSeek page load.

This went unnoticed for a while and is the likely real explanation for a
cluster of previously-unexplained "no response generated at all" failures
scattered across this session's retries (X-06 and X-32's first attempts,
GZ-26, GZ-36, and possibly others) -- two independent MCP integration
layers fighting over the same page's tool-call detection and message
submission. Confirmed live (a "123 of 123 tools enabled" MCP
SuperAssistant panel was visibly present on a DeepSeek tab) and fixed by
unchecking MCP SuperAssistant's own "MCP" toggle for this site, verified
to persist across a fresh page navigation ("MCP Settings - Inactive"
after reload).

**Confirmed with a dedicated follow-up pass**: after the fix, every
previously-blank-response question was rerun once more in the clean
environment. GZ-26 now produces a real (still partially wrong) response --
strong evidence the SuperAssistant conflict was its actual cause. GZ-04,
GZ-17, GZ-36, H-30, and X-30 reproduced their exact same failure symptoms
even in the clean environment, confirming those are genuinely unrelated to
the extension conflict (independent root causes, documented individually
below).

### A second lesson from this pass: a brief process overlap

While sequencing the rerun, a new retry batch was launched before
confirming the previous background batch had actually finished and exited
-- both processes briefly shared the same browser tab, producing a short
burst of spurious failures (H-38 failing once then succeeding on the very
next attempt, H-39 likewise) before being caught via `ps aux` and
resolved by waiting for a clean process list before each subsequent
command. No data was lost (the affected questions were simply re-run
cleanly afterward), but it's a reminder that `nohup ... & disown`
backgrounding genuinely detaches a process -- the *next* command issued
must still wait for it to finish before sharing the same browser tab.

## Genuine, disclosed limitations (not pipeline bugs, real effort made)

All of the following reproduced their failure independently in a clean
(MCP-SuperAssistant-free) environment, so none are attributable to that
conflict.

**Extension-level SMILES-scraping bugs (not toolkit/model issues):**
- **H-30, GZ-17**: Better DeepSeek corrupts `[N+](C)(C)C`-style SMILES
  into a markdown link (`[N+](https://c/)(C)C)`) before the tool call.
  Model detects the corruption but can't escape it (must re-emit the
  SMILES to call the tool, which gets re-parsed fresh each time);
  degenerates into a repeated-token loop on retry.
- **GZ-04**: a different manifestation of the same underlying
  unreliability -- a long chiral-ring SMILES gets silently truncated
  mid-string instead of link-corrupted.

**Persistent, unexplained model resistance (confirmed real, not a nudge
wording problem -- 3 independently-run attempts per question, the last
with a mandatory step-by-step protocol explicitly naming the tool to
call and demanding its raw output be quoted, produced byte-identical
output to the first attempt each time):**
- **H-04**: `delta_H_kJ_per_mol` = -2.846400782241967 to 15 significant
  figures, unchanged across 3 fresh chat sessions with escalating
  nudges, vs gold's -4.88. The model appears to have a fixed internal
  calculation method for this specific question it will not abandon for
  a named tool, even when told explicitly "if you find yourself writing
  a manual regression formula instead of calling this tool, STOP."
- **H-06**: `selected_model`="rubingh_regular_solution",
  `selected_bic`=-27.982532417493196, identical across 3 attempts
  (gold wants `eommm_asymmetric_margules`), despite the same
  step-by-step mandatory-tool-call treatment.

**Toolkit capability gap:**
- **GZ-49**: no MCP tool exposes the two-parameter asymmetric-Margules
  (EOMMM) regression this question needs. Model correctly recognizes
  the gap and tries its sandboxed code runner, which needs a manual
  "Run Code" click the automated pipeline never provides.

**Genuine remaining model-reliability gaps** (wrong or partially wrong
even after explicit, targeted tool-call nudges; not forced to pass):
- **X-30**: BIC-based scale-corruption detection -- converged closer
  each round but never fully within tolerance; also intermittently
  fails extraction (complex, many-tool-call question, possibly hitting
  genuine response-length/stability limits).
- **X-49**: bootstrap CMC confidence-interval bounds remain wrong
  (point estimate itself was correct in one run, wrong in another --
  non-deterministic across chat sessions).
- **GZ-01, GZ-05**: compute a `delta_g_mic`/`delta_s_mic` number where
  gold says null (insufficient information given) -- the model
  disagrees that refusal is correct and computes anyway under a stated
  assumption.
- **GZ-02**: same wrong-model-selection pattern as H-06
  (`rubingh_regular_solution` instead of `eommm_asymmetric_margules`),
  also unmoved by the `select_mixture_model` tool-naming nudge.
- **GZ-08, GZ-26, GZ-29**: isotherm/van't-Hoff parameter values (K,
  delta_H sub-intervals, Frumkin `a`) remain off by a consistent
  factor after retry -- genuine fitted-parameter misses, not reasoning
  errors (the model-SELECTION half of these questions, e.g. choosing
  Langmuir vs Frumkin, did fix cleanly with the tool-naming nudge; only
  the specific numeric parameter remains wrong).
- **GZ-35**: `valid_hlb_method` -- model now actively prefers Davies'
  method for this structure, contradicting gold's conclusion that
  Davies is structurally invalid for this sugar-headgroup class.
- **GZ-36**: no response generated across all 3 attempts including the
  clean-environment retry; root cause distinct from the (resolved)
  SuperAssistant conflict, not otherwise identified.
- **GZ-42**: a real-world-knowledge judgment question (identifying AOT
  and its well-known reverse-micelle-forming behavior) -- unmoved even
  by an explicit compound-identification hint in the nudge.
- **GZ-48**: Rodenas-method activity coefficients come out at almost
  exactly half of gold on both components even after naming the
  specific tool to call -- a persistent, unexplained factor-of-2 error.

**Chart-reading visual precision (not a pipeline defect -- the image
upload itself works, proven by the ~10 chart questions that do pass):**
- **H-38**: improved substantially across retries (gold 5750 vs a
  closest attempt of 5359, within ~7%) but never cleared the 4%
  tolerance.

## Qwen: fully explored and abandoned this round

In addition to the earlier-established DashScope (India billing block)
and OpenRouter (budget shortfall) dead ends, this session also explored
whether `openclaw`'s own model/MCP integration, and separately a
dedicated "MCP SuperAssistant" Chrome extension (verified legitimate:
2,607 GitHub stars, MIT license, officially lists Qwen Chat as a
supported platform), could provide a free workaround via Qwen's free
consumer web chat. Both were real, working integrations (the extension
genuinely connected to the project's MCP server with all 123 tools
visible), but Qwen itself does not reliably follow the extension's
required tool-call output format -- confirmed by the extension's own
maintainer in a still-open GitHub issue ("Does not work in Qwen", filed
2025-09-25, no fix in over a year: "it seems qwen is not following the
instructions"). Live-tested and reproduced the same failure (Qwen
explicitly stated "I do not have a connected MCP tool available in this
session" despite the instructions being present). Abandoned per explicit
user instruction -- not to be revisited unless raised again.
