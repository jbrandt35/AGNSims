import os
import json
from Setup import run_bash

os.chdir("../runs/")

for BBH_separation_dir in [directory for directory in os.listdir() if "BBH_separation" in directory]:
    for perturber_separation_dir in [directory for directory in os.listdir(BBH_separation_dir) if
                                     "perturber_separation" in directory]:
        for run_number in os.listdir(os.path.join(BBH_separation_dir, perturber_separation_dir)):

            try:
                outcome = \
                json.load(open(os.path.join(BBH_separation_dir, perturber_separation_dir, run_number, "outcome.json")))[
                    "Result"]
            except FileNotFoundError:
                outcome = "Didn't Finish"

            if "Collision Encountered" not in outcome:

                os.system(f"rm -r {BBH_separation_dir}/{perturber_separation_dir}/{run_number}/*")

                BBH_separation = float(BBH_separation_dir.split("_")[-1])
                perturber_separation = float(perturber_separation_dir.split("_")[-1])

                job = f"{BBH_separation} {perturber_separation} {run_number}"
                run_bash("build.sh", job)
