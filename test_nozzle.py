"""
test_nozzle.py

Pytest tests for your ASTE 404 mini_project nozzle scripts.

Covers:
1) Inverting area-Mach relation (subsonic and supersonic)
2) Isentropic pressure ratio decreases as Mach increases
3) Exit state from nozzle geometry is physically reasonable
4) Thrust coefficient C_F is finite and positive for a typical case
5) Expansion regime classification behaves as expected
"""

import math
import pytest

from flow_relations import (
    area_ratio_relation,
    mach_area_relation,
    pressure_mach_relation,
    temperature_mach_relation,
)

from nozzle_geometry import (
    exit_area_relation,
    C_F_from_geometry,
)

from command_control import expansion_ratio


# -----------------------------
# Area–Mach inversion tests
# -----------------------------

def test_area_mach_inversion_supersonic():
    gamma = 1.4
    M_true = 2.0
    A = float(area_ratio_relation(M_true, gamma))

    M_solved = mach_area_relation(A, gamma, branch="Supersonic", tol=1e-12)
    assert abs(M_solved - M_true) < 1e-8


def test_area_mach_inversion_subsonic():
    gamma = 1.4
    M_true = 0.30
    A = float(area_ratio_relation(M_true, gamma))

    M_solved = mach_area_relation(A, gamma, branch="Subsonic", tol=1e-12)
    assert abs(M_solved - M_true) < 1e-8


# -----------------------------
# Isentropic ratios behavior
# -----------------------------

def test_pressure_ratio_monotonic():
    gamma = 1.4
    p1 = float(pressure_mach_relation(1.2, gamma))
    p2 = float(pressure_mach_relation(2.0, gamma))

    assert p2 < p1, "Expected p/p0 to decrease as Mach increases"


def test_temperature_ratio_in_range():
    gamma = 1.4
    T = float(temperature_mach_relation(2.0, gamma))
    assert 0.0 < T < 1.0


# -----------------------------
# Exit state + thrust coefficient checks
# -----------------------------

def test_exit_state_physical_supersonic():
    gamma = 1.4
    Ae_At = 10.0

    st = exit_area_relation(Ae_At, gamma, branch="Supersonic")

    assert st["Ae_At"] == pytest.approx(Ae_At)
    assert st["Mach_exit"] > 1.0
    assert 0.0 < st["pe_p0"] < 1.0
    assert 0.0 < st["Te_T0"] < 1.0


def test_thrust_coefficient_reasonable():
    gamma = 1.4
    Ae_At = 10.0
    pa_p0 = 0.02

    C_F = C_F_from_geometry(gamma, Ae_At, pa_p0, branch="Supersonic")

    assert math.isfinite(C_F)
    assert C_F > 0.0


# -----------------------------
# Expansion regime classification
# -----------------------------

def test_underexpanded_classification():
    regime, _ = expansion_ratio(pe_p0=0.10, pa_p0=0.05, rtol=0.02)
    assert regime == "under-expanded"


def test_overexpanded_classification():
    regime, _ = expansion_ratio(pe_p0=0.03, pa_p0=0.05, rtol=0.02)
    assert regime == "over-expanded"


def test_ideally_expanded_classification():
    regime, _ = expansion_ratio(pe_p0=0.051, pa_p0=0.05, rtol=0.05)
    assert regime == "ideally-expanded"
