import pandas as pd
import matplotlib.pyplot as plt

# Read the data from a text file
# Make sure the file has a header line: sigma_input PD psi6 seed
df = pd.read_csv("results_area_0.8", sep=r"\s+")
#results_psi6_plot
# area_corrected_results
# Drop any rows with missing data (some of your lines have only sigma values)
df = df.dropna(subset=["sigma_input", "PD", "psi6"])

# 
grouped = df.groupby("sigma_input").agg(
    mean_PD=("PD", "mean"),
    mean_psi6=("psi6", "mean"),
    std_psi6=("psi6", "std"),
    count=("psi6", "count")
).reset_index()


plt.plot(grouped["mean_PD"], grouped["mean_psi6"], '^-', color='dodgerblue', label="ψ₆ vs Polydispersity")
plt.fill_between(
    grouped["mean_PD"],
    grouped["mean_psi6"] - grouped["std_psi6"],
    grouped["mean_psi6"] + grouped["std_psi6"],
    color='dodgerblue',
    alpha=0.2,       # translucent shading
    linewidth=0
)

plt.xlabel("Polydispersity (%)")
plt.ylabel("Global ψ₆")
plt.title("Order Parameter ψ₆ vs. Polydispersity - Area Corrrected")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

