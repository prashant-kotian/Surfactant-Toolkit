"""Manual corrections where the keyword rules misjudge a question.

Keyed by SB id -> (type, subconcept). Curated against the full 200-row
classification listing on 2026-10-19.
"""
OVERRIDES = {
    "SB-003": ("conceptual", "identification_terminology"),   # gemini family/charge ID; hlb field is null
    "SB-004": ("conceptual", "identification_terminology"),   # charge-classification at pH near pKa
    "SB-040": ("analytical", "cmc_measurement"),               # chart dip-then-rise read-off (image-excluded)
    "SB-094": ("analytical", "micelle_structure"),             # SLS-vs-DLS reliability forensics
    "SB-111": ("computational", "kinetics_transport"),         # Gillespie run, not a judgement
    "SB-114": ("computational", "micelle_structure"),          # Metropolis cluster sim, single surfactant
    "SB-118": ("computational", "micelle_structure"),
    "SB-123": ("computational", "micelle_structure"),
    "SB-126": ("analytical", "data_fitting"),                  # corrupt-row forensics
    "SB-128": ("analytical", "data_fitting"),                  # duplicate-concentration exclusion rule
    "SB-135": ("computational", "identification_terminology"), # pH speciation vocabulary
    "SB-137": ("computational", "packing_morphology"),         # dual-geometry shape report
    "SB-138": ("computational", "packing_morphology"),         # CPP transition algebra
    "SB-140": ("computational", "micelle_structure"),
    "SB-144": ("analytical", "foaming_films"),                 # documented AOT Cc (bending modulus) conflict
    "SB-149": ("computational", "cmc_measurement"),             # bootstrap CMC feeding dG (cmc is the task)
    "SB-005": ("computational", "micelle_thermodynamics"),  # interval van't Hoff delta_H
    "SB-011": ("analytical", "packing_morphology"),          # CPP vs measured-area morphology check
    "SB-028": ("computational", "micelle_thermodynamics"),   # delta_S_mic regimes
    "SB-036": ("analytical", "micelle_thermodynamics"),      # ln(CMC) vs 1/T slope reading (image-excluded)
    "SB-049": ("computational", "mixed_micelles"),           # EOMMM W12/W21 activity coefficients
    "SB-147": ("analytical", "mixed_micelles"),              # dual-method CMC/beta/dG panel
    "SB-150": ("computational", "kinetics_transport"),       # Gillespie relaxation-rate pipeline
    "SB-158": ("computational", "mixed_micelles"),           # Maeda mixed-micelle free energy
    "SB-169": ("conceptual", "identification_terminology"),  # sourced reference-value lookup
    "SB-175": ("computational", "adsorption"),               # Nagarajan equilibrium area (ionic)
    "SB-190": ("computational", "micelle_structure"),          # DLS size inversion
    "SB-195": ("computational", "kinetics_transport"),         # 1-D transport breakthrough half-time
}
