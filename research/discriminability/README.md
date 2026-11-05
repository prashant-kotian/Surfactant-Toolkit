# Discriminability probe (exploratory)

`probe_mechanisms.py` asks whether dilatational rheology can distinguish between
relaxation mechanisms in a realistic frequency window at realistic noise.

This is a feasibility probe for the candidate in `NEXT_PAPER_GAP_ANALYSIS.md` §2,
not a result. It exists to check the question has an answer before committing
effort, because the previous project's failure mode was committing first.

Run: `PYTHONPATH=. python -m research.discriminability.probe_mechanisms`

Caveat: the barrier and reorientation models here are deliberately simple
serial-relaxation stand-ins, not the full Fainerman reorientation/compression or
Danov-Kralchevsky micellar treatments. They establish that the observable carries
information; they do not establish the exact numbers a paper would report.
