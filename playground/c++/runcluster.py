import os
import sys

FOLDER = "../../instances/official/ds/exact"
BINARY = "./build/"
cwd = os.getcwd()
for file in os.listdir(FOLDER):
    if file.endswith(".gr"):
        cmd = ["qsub", "-N", file, "-l", "bc4", "-e", "~/error.log", "runscript.sh", "./build/dsexact", os.path.join(FOLDER, file), cwd]
        cmd = " ".join(cmd)
        print(cmd)
        os.system(cmd)