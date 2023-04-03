import numpy as np

#This exists as a way to share variables between the Tools and Simulations modules


spin = [None, None, None]
mode = None

total_time_steps_completed = 0


outcome_record = {"Result": None, "Minimum Distance Between Binary COM and Perturber": np.inf,
                  "Minimum Distance Between BBHs": np.inf, "Minimum t_GW": np.inf, "Minimum relative t_GW": np.inf, "Events": []}

