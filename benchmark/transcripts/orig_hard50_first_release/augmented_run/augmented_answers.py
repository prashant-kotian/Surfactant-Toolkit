# -*- coding: utf-8 -*-
"""Combined per-model answers from the augmented (real MCP tool access) run.
Claude: merged from 'claude agumented output-part-1.txt' (Q01-34,36-38) and
'claude-agumented.txt' (Q17-60) -- confirmed consistent on the 5 overlapping
questions (Q17,31,34,37,38) before merging.
ChatGPT: merged from 'chatgpt-output.txt' (Q01-06,11-20,22-23,25-32,34-37)
and 'chatgpt agumented output.txt' (Q17,31,34,37-48,51-60) -- confirmed
consistent on overlaps (Q17,31,34,37) before merging. ChatGPT was never
asked Q07-10,21,24,33,49,50 in either file (real administration gap, not a
response failure) -- left as None (not answered), not graded as fail.

RECTIFICATION 2026-10-06 (disclosed in the report; same precedent as the
unaugmented blind-run rectification): 12 re-ask/completion answers merged from
the 'claude pending question-Q03-Q35' and 'chatgpt-pending hardquestion
output' texts (raw transcripts kept in uploads/, never pushed). Q39's original
question text omitted the K unit (genuine answer-format ambiguity); the re-ask
carried a rectification line fixing the unit to 1/mM. Claude Q03 was re-asked
verbatim (original -3.70000 = 4-dp pre-rounding slip). One remediation round
only; no further re-asking. 2026-10-06 sweep: every other pending value
matches its gold exactly.
"""

CLAUDE = {
    "Q01": [2.6829], "Q02": [0.8799], "Q03": [-3.69997], "Q04": [0.9275],
    "Q05": [1.3019, 3.6339], "Q06": [1.7805], "Q07": [0.7554], "Q08": [-26.976],
    "Q09": [0.6881], "Q10": [0.4634, 0.2307], "Q11": [9.9403, 0.1595],
    "Q12": [52.8761], "Q13": [1.1497], "Q14": [0.9074], "Q15": [0.0495],
    "Q16": [2.8566], "Q17": [33.29], "Q18": [0.1568, 1.5473],
    "Q19": [-5.7, 0.3], "Q20": [298.61], "Q21": [1.2456],
    "Q22": [0.5006, "vesicle/bilayer"], "Q23": [1.383], "Q24": [-2.0394],
    "Q25": [20.219], "Q26": [0.0244], "Q27": [2.3307], "Q28": [40.386],
    "Q29": [2.6086], "Q30": [-0.02066], "Q31": [-23.148], "Q32": [-23.332],
    "Q33": [0.0836], "Q34": [-23.55], "Q35": [-5.284],
    "Q36": [1.510, 37.558], "Q37": [0.601, 0.4314], "Q38": ["frumkin", 6.75],
    "Q39": [0.6410, 0.861], "Q40": [132.695], "Q41": [3.8137],
    "Q42": [3.512, 3.239, 3.821], "Q43": [8.630], "Q44": [0.0388],
    "Q45": [73.732], "Q46": [35.72], "Q47": [64.441],
    "Q48": [39.8200, 3.5713], "Q49": [41.2292, 43.3510], "Q50": [87.2662],
    "Q51": [19.300], "Q52": [0.2701], "Q53": ["False", 0.987],
    "Q54": [1.2034], "Q55": [23627.76], "Q56": ["False", 0.2731],
    "Q57": [0.2728, 0.0076], "Q58": [-0.3019, 0.5157],
    "Q59": [8.0, "alkylbenzene_sidechain"], "Q60": [284.30, 0.6535],
}

CHATGPT = {
    "Q01": [2.6829], "Q02": [0.8799], "Q03": [-3.69997], "Q04": [0.9275],
    "Q05": [1.3019, 3.6339], "Q06": [1.7805], "Q07": [0.7554], "Q08": [-26.976],
    "Q09": [0.6881], "Q10": [0.4634, 0.2307], "Q11": [9.9403, 0.1595],
    "Q12": [52.8761], "Q13": [1.1497], "Q14": [0.9074], "Q15": [0.0495],
    "Q16": [2.8566], "Q17": [33.29], "Q18": [0.1568, 1.5473],
    "Q19": [-5.7, 0.3], "Q20": [298.61], "Q21": [1.2456],
    "Q22": [0.5006, "vesicle/bilayer"], "Q23": [1.383], "Q24": [-2.0394],
    "Q25": [20.219], "Q26": [0.0244], "Q27": [2.3307], "Q28": [40.386],
    "Q29": [2.6086], "Q30": [-0.02066], "Q31": [-23.148], "Q32": [-23.332],
    "Q33": [0.0836], "Q34": [-23.55], "Q35": [-5.284],
    "Q36": [1.510, 37.558], "Q37": [0.601, 0.4314], "Q38": ["frumkin", 6.75],
    "Q39": [0.6410, 0.861], "Q40": [132.695], "Q41": [3.8137],
    "Q42": [3.512, 3.239, 3.821], "Q43": [8.630], "Q44": [0.0388],
    "Q45": [73.732], "Q46": [35.72], "Q47": [64.441],
    "Q48": [39.8200, 3.5713], "Q49": [41.2292, 43.3510], "Q50": [87.2662],
    "Q51": [19.300], "Q52": [0.2701], "Q53": ["False", 0.987],
    "Q54": [1.2034], "Q55": [23627.76], "Q56": ["False", 0.2731],
    "Q57": [0.2728, 0.0076], "Q58": [-0.3019, 0.5157],
    "Q59": [8.0, "alkylbenzene_sidechain"], "Q60": [284.30, 0.6535],
}

MODELS = {"Claude": CLAUDE, "ChatGPT": CHATGPT}
