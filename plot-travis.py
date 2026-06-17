import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "power_spectrum_global.csv",
    sep=";",
    comment="#",
    header=None,
    names=["wavenumber_cm-1", "spectrum_K_cm", "integral_K"],
    skipinitialspace=True,
)

plt.plot(df["wavenumber_cm-1"], df["spectrum_K_cm"])
plt.xlim(0, 300)   # focus on InAs phonon region
plt.xlabel("Wavenumber / cm$^{-1}$")
plt.ylabel("Power spectrum / K cm")
plt.tight_layout()
plt.show()
