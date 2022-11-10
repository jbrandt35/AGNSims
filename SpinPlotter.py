import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

rebound_data = pd.read_csv("spin.dat", header = None)
matlab_data = pd.read_csv("psi.txt", delimiter = " ", header = None)

t_to_plot = min(max(matlab_data.iloc[:, 0].array.tolist()), max(rebound_data.iloc[:, 0].array))

matlab_data = matlab_data.loc[matlab_data[0] <= t_to_plot]
rebound_data = rebound_data.loc[rebound_data[0] <= t_to_plot]

matlab_times = matlab_data.iloc[:, 0].array
matlab_longitudes = matlab_data.iloc[:, 1]
matlab_longitude_derivative = np.gradient(matlab_longitudes, matlab_times)

rebound_times = rebound_data.iloc[:, 0].array
rebound_longitudes = rebound_data.iloc[:, 2]
rebound_longitude_derivative = np.gradient(rebound_longitudes, rebound_times)


fig, axs = plt.subplots(2)

derivative_plot = axs[1].plot(rebound_times, rebound_longitude_derivative, label = "Rebound")
axs[1].plot(matlab_times, matlab_longitude_derivative, label = "MATLAB")
axs[1].set_xlim([0, min(max(matlab_times.tolist()), max(rebound_times.tolist()))])
axs[1].set_title("Time Derivative of Spin Longitude over Time")
axs[1].set_xlabel("Time [yr]")
axs[1].set_ylabel("Derivative [1/yr]")
axs[1].set_ylim([0.103-0.002, 0.103+0.002])
axs[1].text(0.25, 0.104, f"Average: {np.format_float_scientific(np.mean(matlab_longitude_derivative), precision = 4)}", color = "orange")
axs[1].text(0.25, 0.102, f"Average: {np.format_float_scientific(np.mean(rebound_longitude_derivative), precision = 4)}", color = "blue")

#obliquity_plot = rebound_data.plot(x = 0, y = 1, xlabel = "Time [yr]", ylabel = "Obliquity of Spin [deg]", title = "Obliquity of Spin Over Time", ax = axs[0], legend = True, sharex = True, label = "Rebound")
longitude_plot = rebound_data.plot(x = 0, y = 2, xlabel = "Time [yr]", ylabel = "Longitude", title = "Spin Longitude Over Time", ax = axs[0], label = "Rebound (Bulirsch-Stoer)", legend = True)
matlab_data.plot(x = 0, y = 1, ax = longitude_plot, label = "MATLAB", sharex = True, xlabel = "Time [yr]", xlim = (0, min(max(matlab_times.tolist()), max(rebound_times.tolist()))), ylim = (-1.8, -1.2))

fig.tight_layout()
plt.savefig("Spin.jpg")
plt.close()


