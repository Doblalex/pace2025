import os
import sys

DEBUG = false
FOLDER = "../../instances/official/ds/exact"
BINARY = "./build-debug/ogdf_dsexact" if DEBUG else "./build-release/ogdf_dsexact"
OUTDIR = "resultsdebug" if DEBUG else "results"

cwd = os.getcwd()
for file in os.listdir(FOLDER):
    if file.endswith(".gr"):
        cmd = ["qsub", "-N", file, "-l", "bc4", "-e", "~/error.log", "-pe", "pthreads", "1", "-l", "h_vmem=8G", "./runscript.sh", BINARY, os.path.join(FOLDER, file), cwd, OUTDIR]
        cmd = " ".join(cmd)
        print(cmd)
        os.system(cmd)