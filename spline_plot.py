import os
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-v0_8-darkgrid")  # красивый стиль

base = os.path.dirname(os.path.abspath(__file__))

files = [
    "spline_graphics_10.dat",
    "spline_graphics_50.dat",
    "spline_graphics_100.dat",
    "spline_graphics_1000.dat",
    "spline_graphics_10000.dat"
]

for fname in files:
    path = os.path.join(base, fname)

    data = np.loadtxt(path, skiprows=1)

    X = data[:, 0]
    Y = data[:, 1]
    D = data[:, 2]

    plt.figure(figsize=(12, 7))

    # --- Графики ---
    plt.plot(X, Y, label="Analytic", linewidth=2.5, color="#FF9500")
    plt.scatter(X, Y, s=12, color="#FF9500")

    plt.plot(X, D, label="DFT", linewidth=2.5, color="#0066FF")
    plt.scatter(X, D, s=12, color="#0066FF")

    # --- Оформление ---
    plt.title(f"Сравнение аналитического решения и DFT\n({fname})", fontsize=16)
    plt.xlabel("X", fontsize=14)
    plt.ylabel("Y", fontsize=14)
    plt.grid(True, alpha=0.3)

    plt.legend(framealpha=0.9, fontsize=12)
    plt.tight_layout()

    plt.show()
