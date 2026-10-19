# SurfBench-200 — classification matrix

Pass % per category, unaugmented→augmented (core-3 lift = mean pp over
Claude/ChatGPT/DeepSeek, the arms run in both conditions). Rows pool every
gradable model×condition cell; `n` is benchmark questions in the category.

## Skill type

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| conceptual | 12 |58.3%→100.0% |58.3%→91.7% |50.0%→100.0% |50.0%→— |50.0%→— |100.0%→91.7% | +41.7 pp |
| computational | 114 |64.0%→96.4% |67.6%→94.6% |62.2%→93.6% |49.1%→— |49.5%→— |73.9%→99.1% | +30.3 pp |
| analytical | 74 |61.5%→92.3% |61.5%→87.7% |60.0%→90.8% |56.9%→— |50.0%→— |87.5%→100.0% | +29.3 pp |

## Sub-concept

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| adsorption | 39 |51.4%→86.5% |58.3%→86.5% |54.1%→88.9% |35.1%→— |45.5%→— |96.2%→100.0% | +32.7 pp |
| analysis_synthesis | 1 |100.0%→0.0% |100.0%→100.0% |0.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| cmc_measurement | 17 |76.9%→100.0% |83.3%→100.0% |46.2%→76.9% |30.8%→— |30.8%→— |40.0%→100.0% | +23.5 pp |
| colloid_electrostatics | 12 |54.5%→100.0% |60.0%→100.0% |45.5%→100.0% |54.5%→— |27.3%→— |66.7%→100.0% | +46.7 pp |
| data_fitting | 3 |66.7%→100.0% |50.0%→100.0% |33.3%→100.0% |33.3%→— |66.7%→— |100.0%→100.0% | +50.0 pp |
| emulsion_hld | 19 |73.7%→100.0% |84.2%→94.7% |89.5%→89.5% |84.2%→— |64.7%→— |100.0%→100.0% | +12.3 pp |
| foaming_films | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| identification_terminology | 10 |60.0%→100.0% |70.0%→90.0% |70.0%→100.0% |70.0%→— |50.0%→— |100.0%→90.0% | +30.0 pp |
| kinetics_transport | 11 |54.5%→100.0% |45.5%→100.0% |54.5%→100.0% |27.3%→— |50.0%→— |66.7%→100.0% | +48.5 pp |
| micelle_structure | 18 |61.1%→100.0% |55.6%→100.0% |50.0%→100.0% |44.4%→— |38.9%→— |53.3%→100.0% | +44.4 pp |
| micelle_thermodynamics | 5 |33.3%→66.7% |33.3%→66.7% |33.3%→66.7% |0.0%→— |0.0%→— |100.0%→100.0% | +33.4 pp |
| mixed_micelles | 43 |65.0%→97.6% |63.4%→87.8% |65.0%→95.0% |56.1%→— |51.2%→— |75.0%→97.6% | +29.0 pp |
| packing_morphology | 13 |76.9%→100.0% |76.9%→92.3% |84.6%→100.0% |92.3%→— |90.9%→— |100.0%→100.0% | +18.0 pp |
| rheology | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| solubilization | 3 |0.0%→66.7% |0.0%→66.7% |0.0%→66.7% |0.0%→— |0.0%→— |100.0%→100.0% | +66.7 pp |
| surface_wetting | 4 |100.0%→100.0% |100.0%→100.0% |75.0%→100.0% |50.0%→— |75.0%→— |100.0%→100.0% | +8.3 pp |

## Conclusions (auto-generated from the matrices above)

- Weakest sub-concepts unaugmented: kinetics_transport (51.5%); colloid_electrostatics (53.1%); adsorption (54.5%); micelle_structure (55.6%); mixed_micelles (64.5%).
- Weakest even augmented (persistent blind spots): adsorption (87.3%); cmc_measurement (92.3%); mixed_micelles (93.4%); emulsion_hld (94.7%); identification_terminology (96.6%).
- Largest augmentation benefit: kinetics_transport (+48.5 pp); colloid_electrostatics (+46.7 pp); micelle_structure (+44.4 pp).
- Where the toolkit helps least: emulsion_hld (+12.3 pp); packing_morphology (+18.0 pp); cmc_measurement (+23.5 pp).
- Skill-type split — conceptual: 58.3%→100.0% (Claude); computational: 64.0%→96.4% (Claude); analytical: 61.5%→92.3% (Claude).
