import os
import sys
import csv

DIR = sys.argv[1]

csvwriter = csv.writer(open(sys.argv[1]+".csv", "w"))

for file in sorted(os.listdir(DIR)):
    with open(os.path.join(DIR, file), "r") as f:
        dssize = ""
        runtime = ""
        for line in f.readlines():
            line = line.strip()
            if line.startswith("c DS solution size: "):
                dssize = int(line.removeprefix("c DS solution size: "))
            if line.startswith("c solve time: "):
                runtime = float(line.removeprefix("c solve time: ").removesuffix("ms"))
        csvwriter.writerow([file.removesuffix(".out"), dssize, runtime])