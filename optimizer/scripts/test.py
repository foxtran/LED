#!/usr/bin/env python3

import glob

import LED
import tagarray as ta
import numpy as np
from scipy.optimize import minimize

def compute_error(x, EF, syss):
    EF.coefficients = x
    error = 0
    for sys in syss:
      error += EF.compute_error(sys, 0.5)
    print(EF.coefficients, ": ", error)
    return error

def main():
    system_files = glob.glob("*.tadump")
    systems = []
    for f in system_files:
      systems.append(LED.System(f))
    EF = LED.EnhancementFactor(1, False, False)
    print(EF.N_free_params)
    x0 = [ 0 ] * EF.N_free_params
    x0[0] = 1
    res = minimize(compute_error, x0, args=(EF, systems), method="Nelder-Mead")
    print(res)

if __name__ == '__main__':
    main()