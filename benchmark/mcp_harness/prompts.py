"""Prompt construction for the SurfBench harness.

The base system + user prompt is IDENTICAL across all conditions -- the only
experimental variable is whether tools are attached to the API call. That is
what makes the unaugmented-vs-augmented comparison clean. The one exception is
that the tool-augmented condition prepends the SurfMCP tool-advertisement note,
because in real MCP use the server advertises its tools to the model as part of
'having the tool available'; withholding that would misrepresent the augmented
condition, not make it fairer.
"""

from __future__ import annotations
import json

SYSTEM_BASE = (
    "You are an expert in surfactant and interfacial physical chemistry. Solve the "
    "problem carefully, showing your reasoning. Pay close attention to units, to which "
    "formula and convention actually applies, and to whether the problem is even "
    "physically solvable as stated -- some problems have no valid solution, and the "
    "correct response in that case is to say so explicitly rather than force a number.\n\n"
    "End your response with a single line beginning exactly with 'FINAL_JSON:' followed "
    "by a JSON object. Put your final numeric answer(s) there. If the problem has no "
    "valid solution, set the JSON to {\"answer\": \"<brief explanation of why no valid "
    "solution exists>\"} instead of inventing numbers."
)

SURFMCP_NOTE = (
    "You have access to a set of validated surfactant-science calculation tools "
    "(Clint/Rubingh/Rosen mixed-micelle theory, Gibbs adsorption, HLB, critical packing "
    "parameter and aggregation number, micellization thermodynamics, Debye length, zeta "
    "potential, hydrodynamic radius, wetting, solubilization, and more). Prefer calling "
    "these tools over computing or recalling these values yourself -- they handle the "
    "exact formulas, constants, unit conversions, and no-solution checks correctly. "
    "Several tools return a categorical classification string alongside a numeric result "
    "(e.g. flow_regime, predicted_morphology, synergy_classification) -- when you report "
    "that classification in your final answer, use the tool's own exact wording verbatim, "
    "not a paraphrase or a different term for the same idea, even if your paraphrase seems "
    "equally correct to you. Some tools require a parameter that encodes real domain "
    "knowledge the question may not give you (e.g. gibbs_surface_excess's system_type, "
    "hlb_from_groups' group names) -- if the question does not state enough to determine "
    "that parameter with confidence, do not guess a plausible-looking value just to get the "
    "tool to run. Report that the quantity cannot be determined from the given information "
    "instead; a specific number produced from an unstated assumption is a wrong answer, not "
    "a correct one."
)


def expected_keys(question: dict):
    gold = question["gold_answer"]
    if isinstance(gold, dict):
        return list(gold.keys())
    if isinstance(gold, (int, float)) and not isinstance(gold, bool):
        return ["value"]
    return None  # category/no-solution question -- no fixed numeric keys


def build_system(condition: str) -> str:
    if condition == "surfmcp":
        return SURFMCP_NOTE + "\n\n" + SYSTEM_BASE
    return SYSTEM_BASE


def build_user(question: dict) -> str:
    parts = [question["question_text"]]
    given = question.get("given_data")
    if given:
        parts.append("\nGiven data (JSON): " + json.dumps(given))
    keys = expected_keys(question)
    if keys is not None:
        parts.append(
            "\nReport your final answer in FINAL_JSON using these keys: "
            + ", ".join(keys)
            + ". Use plain numeric values (not strings) for numeric quantities."
        )
    else:
        parts.append(
            "\nReport your final answer in FINAL_JSON. If the problem is solvable, give "
            "the requested numeric quantities; if it is not, use {\"answer\": \"...\"} to "
            "explain why no valid solution exists."
        )
    return "\n".join(parts)
