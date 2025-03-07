import os
import sys

FOLDER = "../../instances/official/ds/exact"
BINARY = "./build/"

for file in os.listdir(FOLDER):
    if file.endswith(".gr"):
