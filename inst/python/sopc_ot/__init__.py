"""Composition-sensitivity analysis for fixed-cost balanced entropic OT."""
from .core import derivatives,predict,solve
from .batch import batch_solve,readout
from .api import CompositionAudit
__version__='0.2.0.dev0'
