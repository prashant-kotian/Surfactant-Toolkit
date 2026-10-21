# SurfBench-200 — classification matrix

Pass % per category, unaugmented→augmented (core-3 lift = mean pp over
Claude/ChatGPT/DeepSeek, the arms run in both conditions). Rows pool every
gradable model×condition cell; `n` is benchmark questions in the category.

## Skill type

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| conceptual | 12 |58.3%→100.0% |58.3%→91.7% |50.0%→100.0% |50.0%→— |50.0%→— |100.0%→91.7% | +41.7 pp |
| computational | 114 |64.3%→98.1% |67.9%→95.6% |62.5%→92.8% |49.6%→— |50.5%→— |73.5%→99.1% | +30.6 pp |
| analytical | 74 |64.9%→92.6% |62.2%→90.5% |59.5%→90.5% |58.1%→— |52.2%→— |87.5%→100.0% | +29.0 pp |

## Sub-concept

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|---|---|
| adsorption | 39 |55.3%→93.3% |56.8%→92.1% |57.9%→86.5% |39.5%→— |51.5%→— |92.6%→100.0% | +34.0 pp |
| analysis_synthesis | 1 |100.0%→0.0% |100.0%→100.0% |0.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| cmc_measurement | 17 |76.5%→100.0% |75.0%→100.0% |47.1%→82.4% |41.2%→— |35.3%→— |57.1%→100.0% | +27.9 pp |
| colloid_electrostatics | 12 |50.0%→100.0% |63.6%→100.0% |41.7%→100.0% |50.0%→— |25.0%→— |66.7%→100.0% | +48.2 pp |
| data_fitting | 3 |66.7%→100.0% |50.0%→100.0% |33.3%→100.0% |33.3%→— |66.7%→— |100.0%→100.0% | +50.0 pp |
| emulsion_hld | 19 |73.7%→94.7% |84.2%→94.7% |84.2%→89.5% |78.9%→— |64.7%→— |100.0%→100.0% | +12.3 pp |
| foaming_films | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| identification_terminology | 10 |60.0%→100.0% |70.0%→90.0% |70.0%→100.0% |70.0%→— |50.0%→— |100.0%→90.0% | +30.0 pp |
| kinetics_transport | 11 |54.5%→100.0% |45.5%→100.0% |54.5%→100.0% |27.3%→— |44.4%→— |66.7%→100.0% | +48.5 pp |
| micelle_structure | 18 |61.1%→100.0% |55.6%→100.0% |50.0%→100.0% |44.4%→— |38.9%→— |53.3%→100.0% | +44.4 pp |
| micelle_thermodynamics | 5 |60.0%→80.0% |60.0%→80.0% |60.0%→80.0% |40.0%→— |50.0%→— |100.0%→100.0% | +20.0 pp |
| mixed_micelles | 43 |66.7%→97.4% |65.1%→88.4% |64.3%→92.9% |55.8%→— |52.4%→— |75.0%→97.6% | +27.5 pp |
| packing_morphology | 13 |76.9%→100.0% |76.9%→92.3% |84.6%→100.0% |92.3%→— |90.9%→— |100.0%→100.0% | +18.0 pp |
| rheology | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— |100.0%→— |100.0%→100.0% | +0.0 pp |
| solubilization | 3 |0.0%→66.7% |0.0%→66.7% |0.0%→66.7% |0.0%→— |0.0%→— |100.0%→100.0% | +66.7 pp |
| surface_wetting | 4 |100.0%→100.0% |100.0%→100.0% |75.0%→100.0% |50.0%→— |75.0%→— |100.0%→100.0% | +8.3 pp |

## Conclusions (auto-generated from the matrices above)

- Weakest sub-concepts unaugmented: colloid_electrostatics (51.4%); kinetics_transport (51.5%); micelle_structure (55.6%); adsorption (56.6%); mixed_micelles (65.4%).
- Weakest even augmented (persistent blind spots): adsorption (90.5%); mixed_micelles (92.7%); emulsion_hld (93.0%); cmc_measurement (94.1%); identification_terminology (96.6%).
- Largest augmentation benefit: kinetics_transport (+48.5 pp); colloid_electrostatics (+48.2 pp); micelle_structure (+44.4 pp).
- Where the toolkit helps least: emulsion_hld (+12.3 pp); packing_morphology (+18.0 pp); mixed_micelles (+27.5 pp).
- Skill-type split — conceptual: 58.3%→100.0% (Claude); computational: 64.3%→98.1% (Claude); analytical: 64.9%→92.6% (Claude).
