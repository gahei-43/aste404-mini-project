"""
This file determines if a nozzle exit gas plume is under-expanded, ideally-expanded, or over-expanded based on pressures.
It will also detect unusual C_F values to guide the user.
This can be ran in the terminal, for example using:
  python command_control.py --gamma 1.4 --Ae_At 10 --pa_p0 0.02
"""

import argparse

from nozzle_geometry import exit_area_relation, C_F_from_geometry

def expansion_ratio(pe_p0, pa_p0, rtol=0.02):
    # Compare exit pressure to ambient pressure.
    if rtol < 0:
        raise ValueError("rtol must be >= 0.02")
    
    denominator = max(abs(pa_p0), 1e-15)
    relative_error = abs(pe_p0 - pa_p0) / denominator

    if relative_error <= rtol:
        return "ideally-expanded", f"pe ~= pa within {rtol*100:.1f}% relative tolerance."
    if pe_p0 > pa_p0:
        return "under-expanded", "pe > pa: nozzle could expand more for this ambient pressure. Iterate and re-run."
    return "over-expanded", "pe < pa: over-expanded and flow may shock. Iterate and re-run."

def C_F_warning(C_F):
    if C_F < 0:
        return "C_F is negative. Iterate and re-run."
    if C_F > 5:
        return "C_F is very large. Iterate and re-run."
    return None

def main():
    parser = argparse.ArgumentParser(description="Ideal nozzle analysis for expansion ratios.")
    parser.add_argument("--gamma", type=float, required=True, help="Specific heat ratio gamma (>1).")
    parser.add_argument("--Ae_At", type=float, required=True, help="Area ratio Ae/At (>=1).")
    parser.add_argument("--pa_p0", type=float, required=True, help="Ambient ratio pa/p0 in [0,1).")
    parser.add_argument(
        "--branch",
        type=str,
        choices=["Subsonic", "Supersonic"],
        default="Supersonic",
        help="Mach inversion to find area ratios.",
    )
    parser.add_argument("--rtol", type=float, default=02, help="Relative tolerance for ideal expansion check.")
    args = parser.parse_args()

    exit_conditions = exit_area_relation (args.Ae_At, args.gamma, branch=args.branch)
    C_F = C_F_from_geometry(args.gamma, args.Ae_At, args.pa_p0, branch=args.branch)
    expans, note = expansion_ratio(exit_conditions["pe_p0"], args.pa_p0, rtol=args.rtol)
    warning = C_F_warning(C_F)

    print("\nExit Condition")
    print(f"Ae/At = {exit_conditions['Ae_At']:.6g}")
    print(f"Me = {exit_conditions['Me']:.10g}")
    print(f"pe/p0 = {exit_conditions['pe_p0']:.10g}")
    print(f"Te/T0 = {exit_conditions['Te_T0']:.10g}")

    print("\Analysis")
    print(f"Regime = {expans}")
    print(f"Note = {note}")

    print("\nPerformance")
    print(f"C_F = {C_F:.10g}")
    if warn:
        print(warning)
    print()

if __name__ == "__main__":
    main()
