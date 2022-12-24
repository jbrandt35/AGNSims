import os
import json

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
                os.system(f"rm -r {BBH_separation_dir}/{perturber_separation_dir}/{run_number}")

                if len([directory for directory in
                        os.listdir(os.path.join(BBH_separation_dir, perturber_separation_dir))]) == 0:
                    os.system(f"rm -r {BBH_separation_dir}/{perturber_separation_dir}")

                if len([directory for directory in os.listdir(BBH_separation_dir)]) == 0:
                    os.system(f"rm -r {BBH_separation_dir}")
