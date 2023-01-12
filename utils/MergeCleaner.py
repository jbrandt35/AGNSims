import os
import json
from Setup import run_bash, return_bash

os.chdir("../runs/")

for BBH_separation_dir in [directory for directory in os.listdir() if "BBH_separation" in directory]:
    for perturber_separation_dir in [directory for directory in os.listdir(BBH_separation_dir) if
                                     "perturber_separation" in directory]:
        for run_number in os.listdir(os.path.join(BBH_separation_dir, perturber_separation_dir)):

            print(f"Looking at {BBH_separation_dir}, {perturber_separation_dir}, {run_number}")

            try:
                outcome = json.load(open(os.path.join(BBH_separation_dir, perturber_separation_dir, run_number, "outcome.json")))["Result"]
                print(f"Outcome file found. The outcome was {outcome}")
            except FileNotFoundError:
                print("No outcome file found, skipping")
                continue

            if return_bash("running_jobs.sh", path = "../utils") < 500:
                if "Collision Encountered" not in outcome:

                    print("Found a finished job that didn't end merged... Deleting Files and Rebuilding")

                    os.system(f"rm -r {BBH_separation_dir}/{perturber_separation_dir}/{run_number}/*")

                    BBH_separation = float(BBH_separation_dir.split("_")[-1])
                    perturber_separation = float(perturber_separation_dir.split("_")[-1])

                    job = f"{BBH_separation} {perturber_separation} {run_number}"
                    run_bash("build.sh", job, path = "../utils")
                else:
                    print("Found a job that did end merged. Ignoring.")
            else:
                print("Queue is Full. Stopping...")
                break
        else:
            continue
        break

    else:
        continue
    break
