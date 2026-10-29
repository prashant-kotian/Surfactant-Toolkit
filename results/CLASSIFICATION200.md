# SurfBench-200 — classification matrix

Pass % per category, unaugmented→augmented (core-3 lift = mean pp over
Claude/ChatGPT/DeepSeek, the arms run in both conditions). Rows pool every
gradable model×condition cell; `n` is benchmark questions in the category.

## Skill type

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Gemini U→A | Qwen U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|
| conceptual | 12 |66.7%→100.0% |58.3%→91.7% |50.0%→100.0% |58.3%→91.7% |50.0%→— | +38.9 pp |
| computational | 114 |66.7%→98.2% |66.7%→95.6% |64.0%→92.1% |51.8%→98.2% |50.0%→— | +29.5 pp |
| analytical | 74 |64.9%→93.2% |62.2%→90.5% |59.5%→90.5% |52.7%→100.0% |58.1%→— | +29.2 pp |

## Sub-concept

| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Gemini U→A | Qwen U→A | Lift (core 3) |
|---|---|---|---|---|---|---|---|
| adsorption | 39 |64.1%→94.9% |59.0%→92.3% |61.5%→84.6% |53.8%→97.4% |41.0%→— | +29.1 pp |
| analysis_synthesis | 1 |100.0%→0.0% |100.0%→100.0% |0.0%→100.0% |100.0%→100.0% |100.0%→— | +0.0 pp |
| cmc_measurement | 17 |76.5%→100.0% |70.6%→100.0% |47.1%→82.4% |35.3%→100.0% |41.2%→— | +29.4 pp |
| colloid_electrostatics | 12 |50.0%→100.0% |58.3%→100.0% |41.7%→100.0% |25.0%→100.0% |50.0%→— | +50.0 pp |
| data_fitting | 3 |66.7%→100.0% |33.3%→100.0% |33.3%→100.0% |66.7%→100.0% |33.3%→— | +55.6 pp |
| emulsion_hld | 19 |73.7%→94.7% |84.2%→94.7% |84.2%→89.5% |68.4%→100.0% |78.9%→— | +12.3 pp |
| foaming_films | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— | +0.0 pp |
| identification_terminology | 10 |60.0%→100.0% |70.0%→90.0% |70.0%→100.0% |50.0%→90.0% |70.0%→— | +30.0 pp |
| kinetics_transport | 11 |54.5%→100.0% |45.5%→100.0% |54.5%→100.0% |54.5%→100.0% |27.3%→— | +48.5 pp |
| micelle_structure | 18 |61.1%→100.0% |55.6%→100.0% |50.0%→100.0% |38.9%→100.0% |44.4%→— | +44.4 pp |
| micelle_thermodynamics | 5 |60.0%→80.0% |60.0%→80.0% |60.0%→80.0% |40.0%→100.0% |40.0%→— | +20.0 pp |
| mixed_micelles | 43 |67.4%→97.7% |65.1%→88.4% |65.1%→93.0% |53.5%→97.7% |55.8%→— | +27.2 pp |
| packing_morphology | 13 |76.9%→100.0% |76.9%→92.3% |84.6%→100.0% |84.6%→100.0% |92.3%→— | +18.0 pp |
| rheology | 1 |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→100.0% |100.0%→— | +0.0 pp |
| solubilization | 3 |0.0%→66.7% |0.0%→66.7% |0.0%→66.7% |0.0%→100.0% |0.0%→— | +66.7 pp |
| surface_wetting | 4 |100.0%→100.0% |100.0%→100.0% |75.0%→100.0% |75.0%→100.0% |50.0%→— | +8.3 pp |

## Conclusions (auto-generated from the matrices above)

- Weakest sub-concepts unaugmented: colloid_electrostatics (50.0%); kinetics_transport (51.5%); micelle_structure (55.6%); adsorption (61.5%); cmc_measurement (64.7%).
- Weakest even augmented (persistent blind spots): adsorption (90.6%); emulsion_hld (93.0%); mixed_micelles (93.0%); cmc_measurement (94.1%); identification_terminology (96.7%).
- Largest augmentation benefit: colloid_electrostatics (+50.0 pp); kinetics_transport (+48.5 pp); micelle_structure (+44.4 pp).
- Where the toolkit helps least: emulsion_hld (+12.3 pp); packing_morphology (+18.0 pp); mixed_micelles (+27.2 pp).
- Skill-type split — conceptual: 66.7%→100.0% (Claude); computational: 66.7%→98.2% (Claude); analytical: 64.9%→93.2% (Claude).
