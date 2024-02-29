#!/usr/bin/env python3

import os
import shutil
import argparse

BASIS_SETS_NAMES = {
"6-31G*" : "6-31Gs"
}

def get_basis_filename(basis_set):
    global BASIS_SETS_NAMES
    filename = basis_set
    if basis_set in BASIS_SETS_NAMES.keys():
      filename = BASIS_SETS_NAMES[basis_set]
    return filename

def parse_args():
    parser = argparse.ArgumentParser(description="Command-line tool for generating input files.")

    parser.add_argument("--WF", choices=["RHF", "UHF"], required=True, help="RHF/UHF")
    parser.add_argument("--basis_set", type=str, required=True, help="basis set")
    parser.add_argument("--file", type=str, required=True, help="file.geom")

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    filename = args.file
    geomname = filename.split("/")[-1].split(".")[0]
    WFtype = "R"
    if args.WF == "UHF":
      WFtype = "U"
    basis = get_basis_filename(args.basis_set)
    maindir = f"{geomname}_{WFtype}LED_{basis}"
    for calc in ["HF","CC","DFT_eca","DFT_ecab","DFT_ecba","DFT_ecb","LED"]:
      os.makedirs(f"{maindir}/{calc}", exist_ok=True)
      header = open(f"inputs/{calc}_MINP").read().replace("%WF%", args.WF).replace("%BASIS%",args.basis_set)
      geometry = open(filename).read()
      out = header + geometry
      open(f"{maindir}/{calc}/MINP","w").write(out)
    shutil.copy("inputs/run.sh",f"{maindir}/run.sh")

if __name__ == "__main__":
    main()
