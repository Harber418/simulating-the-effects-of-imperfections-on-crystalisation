import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Read the data from a text file
# Make sure the file has a header line: sigma_input PD psi6 seed
df2 = pd.read_csv("results_area_0.8", sep=r"\s+")
df4 = pd.read_csv("area_corrected_results", sep=r"\s+")
df = pd.read_csv("results_psi6_plot", sep=r"\s+")
df3 = pd.read_csv("results_0.89", sep=r"\s+")
#results_psi6_plot
# area_corrected_results
# Drop any rows with missing data (some of your lines have only sigma values)
df = df.dropna(subset=["sigma_input", "PD", "psi6"])
df2 = df2.dropna(subset=["sigma_input", "PD", "psi6"])
df3 = df3.dropna(subset=["sigma_input", "PD", "psi6"])
df4 = df4.dropna(subset=["sigma_input", "PD", "psi6"])
# Group by sigma_input (each unique size ratio)

# Group by sigma_input (each unique size ratio)
grouped = df.groupby("sigma_input").agg(
    mean_PD=("PD", "mean"),
    mean_psi6=("psi6", "mean"),
    std_psi6=("psi6", "std"),
    count=("psi6", "count")
).reset_index()


# Plot ψ6 vs PD (or sigma)
plt.errorbar(
    grouped["mean_PD"],            # x-axis
    grouped["mean_psi6"],          # y-axis
    yerr=grouped["std_psi6"] / np.sqrt(grouped["count"]),      # error bars standard error on the mean 
    fmt='^-',                      # points with lines
    capsize=4,
    elinewidth=1,
    alpha=1.0,
    ecolor='black',
    label="ψ₆ vs Polydispersity - Fixed box "
)


plt.xlabel("Polydispersity (%)")
plt.ylabel("Global ψ₆")
plt.title("Order Parameter ψ₆ vs. Polydispersity")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()


plt.savefig("/home/oskar_harber/Desktop/plots/psi6_vs_PDF.png", dpi=600)

plt.show()


