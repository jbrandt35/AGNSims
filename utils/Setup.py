import os
import numpy as np
from time import sleep
import sys

perturber_separation_range = (2.5, 3.5, 0.25)
BBH_separation_range = (0.1, 0.4, 0.02)
num_in_each = 10


perturber_separation_range = np.array(list(range(int(100 * perturber_separation_range[0]), int(100 * (perturber_separation_range[1] + perturber_separation_range[2])), int(100 * perturber_separation_range[2])))) / 100
BBH_separation_range = np.array(list(range(int(100 * BBH_separation_range[0]), int(100 * (BBH_separation_range[1] + BBH_separation_range[2])), int(100 * BBH_separation_range[2])))) / 100


def return_bash(script, path = "utils"):
    working_dir = os.getcwd()
    os.chdir(path)
    output_stream = os.popen("./" + script)
    output = output_stream.read().strip()
    output_stream.close()
    os.chdir(working_dir)
    return int(output)


def run_bash(script, options, path = "utils"):
    working_dir = os.getcwd()
    os.chdir(path)
    os.system(f"sh {script} {options}")
    os.chdir(working_dir)


def build_queue(BBH_separations, perturber_separations):
    q = []
    for BBH_separation in BBH_separations:
        for perturber_separation in perturber_separations:
            for i in range(1, num_in_each + 1):
                label = f"{np.round(BBH_separation,2)} {np.round(perturber_separation,2)} {i}"
                q.append(label)
    return q


def print_to_stdout(text):
    sys.stdout.write("\n" + text)


queue = build_queue(BBH_separation_range, perturber_separation_range)

while len(queue) > 0:
    print_to_stdout(f"There are {len(queue)} jobs in the queue")
    number_of_jobs_submittable = 500 - return_bash("running_jobs.sh")
    print_to_stdout(f"There are {number_of_jobs_submittable} slots to submit \n")
    for i in range(0, number_of_jobs_submittable):
        try:
            job = queue.pop()
            run_bash("build.sh", job)
        except IndexError:
            print_to_stdout("The queue is finished!")
            break
    sleep(60)








