"""Makes `research.laod_tau2` importable as a package for its own tests."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
