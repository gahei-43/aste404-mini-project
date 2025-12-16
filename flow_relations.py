"""
The purpose of this file is to solve the isentropic compressible flow relations equations with user-defined variables. 
These flow relations are functions of gamma, mach numbers, pressures, temperatures, and areas.
There is also a section that will tell the user if the flow is between subsonic and super sonic range.
Each calculation has error messages built in.
Follow along with the error codes to understand how the relationships are established.
Conversions from nozzle geometry into flow relations is the main purpose of this code.
The flow of code is as follows:
- Convert nozzle geometry Ae/A* into exit Mach number Me.
- Then convert Mach number -> pressure ratio pe/p0 and temperature ratio Te/T0.

The numerical method used in this code:
- Bisection method to solve Mach from the area-Mach relation.
"""

import numpy as np

def area_ratio_relation(Mach_number, gamma):
    # A/A* = (1/M)*[(2/(gamma+1))*(1+(gamma-1)/2*M^2)]^((gamma+1)/(2(gamma-1)))
    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")
    if np.any(np.asarray(Mach_number) <= 0):
        raise ValueError("Mach number M must be > 0. Iterate and re-run.")

    Mach_number = np.asarray(Mach_number, dtype=float)
    mach_term_first = 2/(gamma+1)
    mach_term_second = 1+(gamma-1)*0.5*Mach_number**2
    mach_exponent_first = (gamma+1)/(2*(gamma-1))
    return (1/Mach_number)*(mach_term_first*mach_term_second)**mach_exponent_first

def pressure_mach_relation(Mach_number, gamma):
    # p/p0 = [1+(gamma-1)/2*Mach_number^2]^(-gamma/(gamma-1))

    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")
    if np.any(np.asarray(Mach_number) < 0):
        raise ValueError("Mach number M must be >= 0. Iterate and re-run.")

    Mach_number = np.asarray(Mach_number, dtype=float)
    return (1+(gamma-1)*0.5*M**2)**(-gamma/(gamma-1))

def temperature_mach_relation(Mach_number, gamma):
    # T/T0 = [1+(gamma-1)/2*M^2]^(-1)

    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")
    if np.any(np.asarray(Mach_number) < 0):
        raise ValueError("Mach number M must be >= 0. Iterate and re-run.")

    Mach_number = np.asarray(Mach_number, dtype=float)
    return (1+(gamma-1)*0.5*Mach_number**2)**(-1)

def mach_area_relation(A_Astar, gamma, branch="Supersonic", tol=1e-10, max_iter=200):
    """
    Solve for Mach number M given A/A* using bisection.
    To get the Mach number, inverting the area/mach relation is necessary.
    tol : float - Bisection stops when bracket size is smaller than tol.
    max_iter : int - Maximum bisection iterations.
    """
    if gamma <= 1:
        raise ValueError("gamma must be > 1. Iterate and re-run.")
    if A_Astar < 1:
        raise ValueError("A/A* must be >= 1. Iterate and re-run.")
    if branch not in ["Subsonic", "Supersonic"]:
        raise ValueError("branch must be 'Subsonic' or 'Supersonic'. Iterate and re-run.")
    if tol <= 0:
        raise ValueError("tol must be > 0. Iterate and re-run.")
    if max_iter < 1:
        raise ValueError("max_iter must be >= 1. Iterate and re-run.")

    # This is only for if at the throat, A/A* = 1 exactly and M = 1.
    if abs(A_Astar - 1) < 1e-14:
        return 1

    # Define f(M) = area_ratio_relation(M)-goal so that f(M) = 0.
    def f(Mach_number):
        return float(area_ratio_relation(Mach_number, gamma)-A_Astar)

    # Ensure that f(lower_mach_boudary) and f(higher_mach_boundary) have opposite signs.
    # That ensures that a root exists in the interval.
    if branch == "Subsonic":
        lower_mach_boudary = 1e-12
        higher_mach_boundary = 1 - 1e-12
    else:
        lower_mach_boudary = 1+1e-12
        higher_mach_boundary = 50

    f_lower_mach_boudary = f(lower_mach_boudary)
    f_higher_mach_boundary = f(higher_mach_boundary)

    # If Supersonic bracket fails, expand the higher_mach_boundary so that it will bracket correctly.
    if branch == "Supersonic" and f_lower_mach_boudary*f_higher_mach_boundary > 0:
        for _ in range(20):
            higher_mach_boundary *= 2
            f_higher_mach_boundary = f(higher_mach_boundary)
            if f_lower_mach_boudary * f_higher_mach_boundary <= 0:
                break
        else:
            raise RuntimeError("Could not bracket root. Try larger higher_mach_boundary value. Iterate and re-run.")

    if branch == "Subsonic" and f_lower_mach_boudary*f_higher_mach_boundary > 0:
        raise RuntimeError("Could not bracket Subsonic root. Iterate and re-run.")

    # This is the bisection solver loop.
    for _ in range(max_iter):
        middle_mach = 0.5*(lower_mach_boudary+higher_mach_boundary)
        f_middle_mach = f(middle_mach)

        # Once the interval is small enough, the solver will halt.
        if abs(higher_mach_boundary-lower_mach_boudary) < tol:
            return middle_mach

        # Half-root solver
        if f_lower_mach_boudary*f_middle_mach <= 0:
            higher_mach_boundary = middle_mach
            f_higher_mach_boundary = f_middle_mach
        else:
            lower_mach_boudary = middle_mach
            f_lower_mach_boudary = f_middle_mach

    raise RuntimeError("Bisection did not converge within max_iter. Iterate and re-run.")