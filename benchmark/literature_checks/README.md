# Literature checks: five published results recomputed with the toolkit

Five recent papers (2022-2025) whose central numbers were recomputed from the values the papers
themselves print, using the toolkit's functions. Each script first reproduces a quantity the paper
got right (to show the inputs were read correctly), then recomputes the claim, with an independent
numpy calculation alongside the toolkit call. Inputs are transcribed from the published tables and
checked against both the publisher PDF and the journal's full-text XML; every input carries its
table or equation number in the script docstring. None of the papers had a published correction
or retraction notice in Crossref on 2026-10-27.

| Case | Paper | Claim examined | Recomputed result | Toolkit functions |
|---|---|---|---|---|
| 1 | Naskar, Colloids Interfaces 9 (2025) 4, doi:10.3390/colloids9010004 | van't Hoff and calorimetric enthalpies of CTAB micellization differ significantly; deltaCp -2.7 vs -0.45 kJ/(mol K) | the printed deltaG(T) varies by 1.3 kJ/mol over 290-323 K, equal to its stated 3 % uncertainty; Kirchhoff fit deltaH(298) 95 % interval [-79, 22] kJ/mol and deltaCp [-7.2, 4.0] kJ/(mol K), both containing the calorimetric values; the printed fit coefficients do not reproduce Table 1 | `gibbs_free_energy_micellization`, `vant_hoff_multi_point_fit` |
| 2 | Sedlarikova et al., ACS Omega 9 (2024) 23209, doi:10.1021/acsomega.3c08917 | Rubingh beta = +3.68, +3.23, +1.72, described as a mild synergism | every mixed CMC lies below both pure CMCs, which requires beta < 0; for the 1:1 mass mixtures beta = -8.7, -6.4, -2.6 (95 % intervals below zero) | `solve_rubingh_x`, `rubingh_beta`, `clint_ideal_cmc` |
| 3 | Sergeev et al., ACS Omega 7 (2022) 3016, doi:10.1021/acsomega.1c06047 | beta = +21.33, antagonism in mixed micelle formation | at the reported mixed CMC the PE65 monomer concentration (46 mg/L) exceeds PE65's own CMC (33 mg/L), which no pseudo-phase model allows; Rubingh's equations have no solution on a mass or molar basis; the printed (x1, beta) satisfies only one of the two equations | `solve_rubingh_x`, `clint_ideal_cmc` |
| 4 | Yaqoob et al., Molecules 28 (2023) 3952, doi:10.3390/molecules28093952 | mixed micellization is exothermic (deltaH -5.0 to -7.1 kJ/mol) and more spontaneous than CPC alone at every temperature | the printed deltaG(T), driven by alpha falling from 0.51 to 0.23, implies deltaH = +51 to +160 kJ/mol at 303-323 K by Gibbs-Helmholtz; the printed deltaH follows only if alpha is held constant; deltaG at 293.15 K recomputes to -40.9 kJ/mol (printed -32.4) | `gibbs_free_energy_micellization`, `vant_hoff_multi_point_fit` |
| 5 | Banjare et al., Heliyon 10 (2024) e35598, doi:10.1016/j.heliyon.2024.e35598 | Gamma_max 12-25 umol/m2, A_min 0.065-0.134 nm2, CPP 255-328, shape change from globule to sphere | the printed areas are below the 0.21 nm2 cross-section of one alkyl chain; from the paper's own CMC, C20 and pi_CMC, Gamma_max = 1.1-2.6 umol/m2, A_min = 0.63-1.46 nm2 and CPP = 0.14-0.33 in every solvent (spherical micelles throughout) | `szyszkowski_surface_tension`, `critical_packing_parameter` |

    python benchmark/literature_checks/run_all.py        # -> case1..5_results.json, literature_checks_summary.json
    python -m pytest benchmark/literature_checks

Each script also runs alone (`python case2_poloxamer_rubingh.py`) and takes a toolkit path as its
first argument; by default it uses this repository's `src/`.

How the five were found: open-access papers from 2019-2026 in Europe PMC were screened by script
for internal consistency of three quantities (free energy of micellization against CMC and
counterion binding, 341 papers; Rubingh parameters against the printed CMCs, 39 papers; area per
molecule against surface excess, 90 papers); flagged tables were read in full and kept only when
the inconsistency changes a stated conclusion and the paper's own correct values could be
reproduced. Case 1 came from a separate reading of the paper.
