"""
This file calculates nozzle geometry with isentropic flow relations. It can find the mach exit number, exit pressure ratio, exit temperature ratio
and ideal thrust coefficient (C_F).

This file calls values from flow_relations.py.
"""

import numpy as np

from flow_relations import (mach_from_area_ratio, pressure_ratio_from_mach,temperature_ratio_from_mach,)

def exit_area_relation(Ae_At, gamma, branch="Supersonic", tol=1e-10, max_iter=200):
    # Calculate the exit conditions when given Ae/At and gamma.

    if Ae_At < 1:
        raise ValueError("Ae/At must be >= 1. Iterate and re-run.")
    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")

    Mach_exit = mach_area_relation(Ae_At, gamma, branch=branch, tol=tol, max_iter=max_iter)
    pe_p0 = float(pressure_mach_relation(Mach_exit, gamma))
    Te_T0 = float(temperature_mach_relation(Mach_exit, gamma))

    return {
        "Ae_At": float(Ae_At),
        "Mach_exit": float(Mach_exit),
        "pe_p0": pe_p0,
        "Te_T0": Te_T0,
    }

def thrust_coefficient_CF(gamma, pe_p0, pa_p0, Ae_At):
    # CF = sqrt((2gamma^2/(gamma-1))*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-(pe/p0)^((gamma-1)/gamma)))+(pe/p0-pa/p0)*(Ae/At)
    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")
    if not (0 <= pe_p0 < 1):
        raise ValueError("pe_p0 must be in [0, 1). Iterate and re-run.")
    if not (0 <= pa_p0 < 1):
        raise ValueError("pa_p0 must be in [0, 1). Iterate and re-run.")
    if Ae_At < 1:
        raise ValueError("Ae/At must be >= 1. Iterate and re-run.")

    coefficient = (2*gamma**2/(gamma-1))*(2/(gamma+1))**((gamma+1)/(gamma-1))
    momentum = np.sqrt(coefficient*(1-pe_p0**((gamma-1)/gamma)))
    pressure = (pe_p0-pa_p0)*Ae_At
    return float(momentum+pressure)

def C_F_from_geometry(gamma, Ae_At, pa_p0, branch="Supersonic", tol=1e-10, max_iter=200):
    # Ae/At -> Me -> pe/p0 -> CF
    exit_conditions = exit_area_relation(Ae_At, gamma, branch=branch, tol=tol, max_iter=max_iter)
    return thrust_coefficient_CF(gamma, exit_conditions["pe_p0"], pa_p0, Ae_At)