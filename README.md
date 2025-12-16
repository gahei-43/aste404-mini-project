# ASTE 404 Mini Project — Ideal Nozzle Analysis

## Overview
This project extends ASTE 404 Homework 5 by reorganizing the original isentropic rocket relations into a small, modular Python codebase, but giving more outputs for isentropic flow relations. All calculations assume steady, one-dimensional, isentropic flow. The files are as follows:

## Files
- flow_relations.py – Isentropic flow relations and bisection-based area–Mach inversion  
- nozzle_geometry.py – Exit conditions and ideal thrust coefficient calculations  
- command_control.py – Expansion regime classification and CLI interface  
- test_nozzle.py – Pytest verification of numerical and physical behavior  

## Numerical Method
The exit Mach number is obtained by numerically inverting the area–Mach relation using the bisection method for both subsonic and supersonic branches. The solvers will guide the user by giving warning messages if values or calculations are incorrect.

## Usage
Run a nozzle analysis:
python command_control.py --gamma 1.4 --Ae_At 10 --pa_p0 0.02

## Testing
Run a test of the entire program:
pytest -q
