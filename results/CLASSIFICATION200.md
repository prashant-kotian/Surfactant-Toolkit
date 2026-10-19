# SurfBench-200 — classification matrix

Pass % per category, unaugmented→augmented (core-3 lift = mean pp over
Claude/ChatGPT/DeepSeek, the arms run in both conditions). Rows pool every
gradable model×condition cell; `n` is benchmark questions in the category.

## Skill type

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| conceptual | 12 |40.0%→100.0% |50.0%→91.7% |50.0%→100.0% |50.0%→— |50.0%→— |100.0%→91.7% | +50.6 pp |
| computational | 114 |62.8%→96.5% |66.7%→94.7% |62.2%→92.8% |49.1%→— |49.5%→— |73.9%→99.1% | +30.8 pp |
| analytical | 74 |64.7%→92.3% |56.9%→87.7% |60.0%→90.8% |56.9%→— |50.0%→— |87.5%→100.0% | +29.7 pp |

## Sub-concept

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| adsorption | 39 |53.3%→86.5% |55.6%→86.5% |54.1%→88.9% |35.1%→— |45.5%→— |96.2%→100.0% | +33.0 pp |
| analysis_synthesis | 1 |100.0%→0.0% |100.0%→100.0% |0.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| cmc_measurement | 17 |76.9%→100.0% |83.3%→100.0% |46.2%→76.9% |30.8%→— |30.8%→— |40.0%→100.0% | +23.5 pp |
| colloid_electrostatics | 12 |54.5%→100.0% |60.0%→100.0% |45.5%→100.0% |54.5%→— |27.3%→— |66.7%→100.0% | +46.7 pp |
| data_fitting | 3 |66.7%→100.0% |50.0%→100.0% |33.3%→100.0% |33.3%→— |66.7%→— |100.0%→100.0% | +50.0 pp |
| emulsion_hld | 19 |76.9%→100.0% |84.2%→94.7% |89.5%→94.7% |84.2%→— |64.7%→— |100.0%→100.0% | +12.9 pp |
| foaming_films | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| identification_terminology | 10 |33.3%→100.0% |60.0%→90.0% |70.0%→100.0% |70.0%→— |50.0%→— |100.0%→90.0% | +42.2 pp |
| kinetics_transport | 11 |50.0%→100.0% |36.4%→100.0% |54.5%→90.9% |27.3%→— |50.0%→— |66.7%→100.0% | +50.0 pp |
| micelle_structure | 18 |58.8%→100.0% |55.6%→100.0% |50.0%→100.0% |44.4%→— |38.9%→— |53.3%→100.0% | +45.2 pp |
| micelle_thermodynamics | 5 |0.0%→66.7% |33.3%→66.7% |33.3%→66.7% |0.0%→— |0.0%→— |100.0%→100.0% | +44.5 pp |
| mixed_micelles | 43 |64.7%→97.6% |61.0%→87.8% |65.0%→92.5% |56.1%→— |51.2%→— |75.0%→97.6% | +29.1 pp |
| packing_morphology | 13 |85.7%→100.0% |69.2%→92.3% |84.6%→100.0% |92.3%→— |90.9%→— |100.0%→100.0% | +17.6 pp |
| rheology | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| solubilization | 3 |0.0%→66.7% |0.0%→66.7% |0.0%→66.7% |0.0%→— |0.0%→— |100.0%→100.0% | +66.7 pp |
| surface_wetting | 4 |100.0%→100.0% |100.0%→100.0% |75.0%→100.0% |50.0%→— |75.0%→— |100.0%→100.0% | +8.3 pp |

## Conclusions (auto-generated from the matrices above)

- Weakest sub-concepts unaugmented: kinetics_transport (46.9%); colloid_electrostatics (53.1%); adsorption (54.4%); micelle_structure (54.7%); mixed_micelles (63.5%).
- Weakest even augmented (persistent blind spots): adsorption (87.3%); cmc_measurement (92.3%); mixed_micelles (92.6%); emulsion_hld (96.5%); kinetics_transport (97.0%).
- Largest augmentation benefit: kinetics_transport (+50.0 pp); colloid_electrostatics (+46.7 pp); micelle_structure (+45.2 pp).
- Where the toolkit helps least: emulsion_hld (+12.9 pp); packing_morphology (+17.6 pp); cmc_measurement (+23.5 pp).
- Skill-type split — conceptual: 40.0%→100.0% (Claude); computational: 62.8%→96.5% (Claude); analytical: 64.7%→92.3% (Claude).
