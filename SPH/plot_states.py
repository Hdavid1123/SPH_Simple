import glob
import numpy as np
import matplotlib.pyplot as plt

# Parámetro de escala para los vectores de velocidad
f = 0.01  

# Cada cuántos pasos graficar
step = 10  # <- puedes cambiar este número

# Buscar todos los archivos state_*.txt
state_files = sorted(glob.glob("state_*.txt"))

if not state_files:
    print("❌ No se encontraron archivos state_*.txt.")
    exit()

# Calcular rango usando solo partículas de frontera
xmin, xmax = np.inf, -np.inf
ymin, ymax = np.inf, -np.inf

for fname in state_files:
    data = np.loadtxt(fname)
    x, y = data[:, 1], data[:, 2]
    ptype = data[:, 11].astype(int)
    mask_boundary = (ptype == 1)
    if np.any(mask_boundary):
        xmin = min(xmin, x[mask_boundary].min())
        xmax = max(xmax, x[mask_boundary].max())
        ymin = min(ymin, y[mask_boundary].min())
        ymax = max(ymax, y[mask_boundary].max())

dx = xmax - xmin
dy = ymax - ymin
margin_factor = 3
xmin -= margin_factor * dx * 0.01
xmax += margin_factor * dx * 0.01
ymin -= margin_factor * dy * 0.01
ymax += margin_factor * dy * 0.01

print(f"📊 Rango centrado en fronteras: x[{xmin}, {xmax}], y[{ymin}, {ymax}]")

# Crear imágenes solo de cada 'step' archivos
for i, fname in enumerate(state_files):
    if i % step != 0:  # <-- esto hace que solo se graficen cada 'step'
        continue

    data = np.loadtxt(fname)
    x, y = data[:, 1], data[:, 2]
    vx, vy = data[:, 3], data[:, 4]
    ptype = data[:, 11].astype(int)

    plt.figure(figsize=(6, 6))
    plt.title(f"Paso {i}")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.gca().set_aspect('equal', adjustable='box')

    mask_fluid = (ptype == 0)
    mask_boundary = (ptype == 1)

    plt.scatter(x[mask_fluid], y[mask_fluid], s=5, c="blue", label="Fluido")
    plt.scatter(x[mask_boundary], y[mask_boundary], s=5, c="black", label="Frontera")

    plt.quiver(
        x[mask_fluid], y[mask_fluid],
        vx[mask_fluid], vy[mask_fluid],
        angles="xy", scale_units="xy", scale=1/f,
        color="red", alpha=0.6
    )

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"frame_{i:04d}.png", dpi=150)
    plt.close()

print("✅ Todas las imágenes fueron generadas en output/")
