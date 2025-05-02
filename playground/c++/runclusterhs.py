import os
import sys

FOLDER = "../../instances/official/hs/exact"
BINARY = "./build/ogdf_hsexact"
OUTDIR = "resultshs"

cwd = os.getcwd()
for file in os.listdir(FOLDER):
    if file.endswith(".hgr"):
        cmd = ["qsub", "-N", f"hs{file}", "-l", "bc4", "-e", "~/error.log", "-pe", "pthreads", "1", "-l", "h_vmem=8G", "./runscript.sh", BINARY, os.path.join(FOLDER, file), cwd, OUTDIR]
        cmd = " ".join(cmd)
        print(cmd)
        os.system(cmd)