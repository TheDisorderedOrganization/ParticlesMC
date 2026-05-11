# msad.jl
#
# On-the-fly Mean Squared Angular Displacement (MSAD) tracker.
# Implements three methods:
#   - Integral  : accumulates incremental rotation vectors every step
#   - Threshold : relays the reference frame when rotation exceeds θ_T
#   - Euler     : compares directly to initial frame, computed on log schedule



