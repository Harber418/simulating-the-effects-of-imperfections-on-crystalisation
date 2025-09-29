
import math

nx, ny = 30, 30
a = 1.0  # lattice constant

with open("hex30x30.txt", "w") as f:
    for i in range(nx):
        for j in range(ny):
            x = i * a + 0.5 * (j % 2) * a
            y = j * math.sqrt(3)/2 * a
            f.write(f"{x:.6f} {y:.6f}\n")