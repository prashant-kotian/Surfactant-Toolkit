"""LAOD / tau2 research module -- forward model and numerical tooling.

RESEARCH CODE, NOT LIBRARY CODE. This package does not ship in the
surfactantkit distribution and is not covered by the package's validation
tiers. It exists to test one falsifiable claim before any of it is promoted
into `src/surfactantkit/`:

    Can the slow micellar relaxation time tau2 be recovered from the
    higher-harmonic content of a large-amplitude oscillatory dilatational
    (LAOD) surface-tension response?

Everything here is written so that it can be *disproved*: the analytic
anchors in `tests/test_anchors.py` are run first, and if they fail the
forward model is wrong and no downstream claim means anything.

Modules
-------
lvdt        : analytic Lucassen-van den Tempel limit + the two micellar
              limiting cases (fast / slow reservoir). Closed form, no PDE.
model       : reaction-diffusion forward model. 1D monomer diffusion to a
              planar adsorbing interface with a linearised Aniansson-Wall
              micellar reservoir and a dilatational surface-area term.
harmonics   : demodulation of the periodic response, with convergence and
              noise-injection helpers.
"""
