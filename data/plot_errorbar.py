import pandas as pd
import matplotlib.pyplot as plt

# Read the data from a text file
# Make sure the file has a header line: sigma_input PD psi6 seed
df = pd.read_csv("results_psi6_plot", sep=r"\s+")

# Drop any rows with missing data (some of your lines have only sigma values)
df = df.dropna(subset=["sigma_input", "PD", "psi6"])

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
    yerr=grouped["std_psi6"],      # error bars
    fmt='^-',                      # points with lines
    capsize=4,
    ecolor='black',
    label="ψ₆ vs Polydispersity"
)

plt.xlabel("Polydispersity (%)")
plt.ylabel("Average ψ₆")
plt.title("Order Parameter ψ₆ vs. Polydispersity")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()

