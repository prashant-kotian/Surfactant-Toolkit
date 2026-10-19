"""Structure descriptor panel — CPU computations a text-only LLM cannot perform.

Every value here requires executing a real cheminformatics algorithm on the
molecular graph (Crippen logP atom-contribution summation, TPSA surface
enumeration, topological indices). Language models reliably *recognize* these
descriptors but cannot *compute* them by hand for nontrivial SMILES -- they
either decline or emit plausible-looking fabricated numbers. This module makes
the computation a one-call worker.

Methods (all as implemented by RDKit, the package's one real dependency):

- logp_crippen      Wildman & Crippen partition coefficient, J. Am. Chem. Soc.
                    123 (2001) 3804-3805 (atom-contribution logP/logS model --
                    the same algorithm that underlies the ESOL solubility
                    estimator already used in solubilization.py).
- tpsa_A2           topological polar surface area, Ertl, Froemmel & Schuff,
                    J. Chem. Inf. Comput. Sci. 40 (2000) 1575-1580.
- molar_refractivity  Le molar refractivity = (n^2-1)/(n^2+2) * M/rho, from
                    atomic increments (RDKit MolMR).
- labute_asa        Labute accessible surface area (CCC algorithm), Labute,
                    J. Mol. Graph. Model. 18 (2000) 46-60.
- bertz_ct          Bertz complexity index, Bertz, J. Am. Chem. Soc. 103
                    (1981) 1266 (branching/connectivity complexity).
- chi0v / kappa1    Kier-Hall connectivity valence chi and Kier shape index
                    (topological indices; Kier & Hall, molecular connectivity
                    theory).
- counts            H-bond donors/acceptors, rotatable bonds, rings, sp3
                    fraction -- RDKit's standard definitions (Lipinski-style
                    rules but computed from the graph, never guessed).

Values are deterministic graph computations: same SMILES -> same outputs, no
fitting, no randomness, milliseconds on any CPU.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class DescriptorPanelResult:
    smiles: str
    canonical_smiles: str
    logp_crippen: float
    tpsa_A2: float
    molar_refractivity_cm3_per_mol: float
    mol_weight_g_per_mol: float
    exact_mass_g_per_mol: float
    labute_asa_A2: float
    bertz_ct: float
    chi0v: float
    kappa1: float
    n_hbd: int
    n_hba: int
    n_rotatable_bonds: int
    n_rings: int
    n_aromatic_rings: int
    fraction_csp3: float
    method: str
    caveat: str

    def as_dict(self) -> dict:
        return asdict(self)


def descriptor_panel(smiles: str) -> DescriptorPanelResult:
    """Compute the full descriptor panel for a SMILES string.

    Raises ValueError with a clear message if the SMILES cannot be parsed --
    never returns guessed values for an invalid structure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "descriptor_panel requires rdkit (pip install rdkit) -- "
            "graph descriptors cannot be computed by string parsing"
        ) from e

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty string")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            f"SMILES could not be parsed: {smiles!r} -- no descriptors are "
            "available for an invalid structure (do not guess)"
        )

    return DescriptorPanelResult(
        smiles=smiles,
        canonical_smiles=Chem.MolToSmiles(mol),
        logp_crippen=round(Descriptors.MolLogP(mol), 4),
        tpsa_A2=round(Descriptors.TPSA(mol), 4),
        molar_refractivity_cm3_per_mol=round(Descriptors.MolMR(mol), 4),
        mol_weight_g_per_mol=round(Descriptors.MolWt(mol), 4),
        exact_mass_g_per_mol=round(Descriptors.ExactMolWt(mol), 6),
        labute_asa_A2=round(Descriptors.LabuteASA(mol), 4),
        bertz_ct=round(Descriptors.BertzCT(mol), 4),
        chi0v=round(Descriptors.Chi0v(mol), 6),
        kappa1=round(Descriptors.Kappa1(mol), 6),
        n_hbd=int(Descriptors.NumHDonors(mol)),
        n_hba=int(Descriptors.NumHAcceptors(mol)),
        n_rotatable_bonds=int(Descriptors.NumRotatableBonds(mol)),
        n_rings=int(Descriptors.RingCount(mol)),
        n_aromatic_rings=int(Descriptors.NumAromaticRings(mol)),
        fraction_csp3=round(Descriptors.FractionCSP3(mol), 6),
        method="RDKit graph algorithms (Wildman-Crippen logP, Ertl TPSA, "
               "Labute ASA, Bertz, Kier-Hall indices)",
        caveat="Computed physicochemical ESTIMATES from molecular structure; "
               "logP is an atom-contribution prediction (typically ~0.3-0.5 log "
               "units, larger for some families such as n-alkanes), not a "
               "measurement. HBD/HBA follow RDKit's definitions.",
    )
