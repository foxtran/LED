#!/usr/bin/env python3

import numpy as np
import tagarray as ta

fdir = "Be_ULED_cc-pVTZ"

ECC = ta.Container.load(f"{fdir}/CFFILE.tadump")
DFTORB = ta.Container.load(f"{fdir}/DFTGRID.tadump")

total = ta.Container(fdir)
print(total.description)
print(DFTORB.keys())
print(ECC.keys())

for k in ECC.keys():
  total[k] = ECC[k]
for k in DFTORB.keys():
  total[k] = DFTORB[k]

if np.all(total['scftype'].data == [82,72,70,32]): # RHF
  for k in [ x for x in total.keys() if ".A." in x ]:
    kA = k.replace(".A.", ".B.")
    total[k] = total[k].data * np.sqrt(0.5)
    total[kA] = total[k]

total.save(total.description + ".tadump_old")

for S in ["A", "B"]:
  nel_key = {"A" : "nalpha", "B" : "nbeta" }[S]
  for orb in range(1, total[nel_key].data[0] + 1):
    total[f'tau.{S}.{orb}'] = (total[f'ORB.{S}.{orb}.X'].data * total[f'ORB.{S}.{orb}.X'].data) + \
                              (total[f'ORB.{S}.{orb}.Y'].data * total[f'ORB.{S}.{orb}.Y'].data) + \
                              (total[f'ORB.{S}.{orb}.Z'].data * total[f'ORB.{S}.{orb}.Z'].data)
    total[f'rho.{S}.{orb}'] = (total[f'ORB.{S}.{orb}.val'].data * total[f'ORB.{S}.{orb}.val'].data)
    del total[f'ORB.{S}.{orb}.X']
    del total[f'ORB.{S}.{orb}.Y']
    del total[f'ORB.{S}.{orb}.Z']
    del total[f'ORB.{S}.{orb}.val']
    del total['gridX'], total['gridY'], total['gridZ']

total.save(total.description + ".tadump")
