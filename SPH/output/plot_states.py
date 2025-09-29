import glob
import numpy as np
import matplotlib.pyplot as plt

# Parámetro de escala para los vectores de velocidad
f = 0.01  

# Buscar todos los archivos state_*.txt
state_files = sorted(glob.glob("state_*.txt"))

if not state_files:
    print("❌ No se encontraron archivos state_*.txt.")
    exit()

# Leer todos los datos para calcular los rangos globales
xmin, xmax = np.inf, -np.inf
ymin, ymax = np.inf, -np.inf

for fname in state_files:
    data = np.loadtxt(fname)
    x, y = data[:, 1], data[:, 2]  # columnas 2 y 3: posiciones
    xmin, xmax = min(xmin, x.min()), max(xmax, x.max())
    ymin, ymax = min(ymin, y.min()), max(ymax, y.max())

print(f"📊 Rango detectado: x[{xmin}, {xmax}], y[{ymin}, {ymax}]")

# Crear imágenes
for i, fname in enumerate(state_files):
    data = np.loadtxt(fname)

    # Estructura: [id, x, y, vx, vy, ax, ay, rho, m, p, h, type]
    x, y = data[:, 1], data[:, 2]
    vx, vy = data[:, 3], data[:, 4]
    ptype = data[:, 11].astype(int)   # <- ahora sí columna correcta (última)

    plt.figure(figsize=(6, 6))
    plt.title(f"Paso {i}")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.gca().set_aspect('equal', adjustable='box')

    # Separar partículas fluido vs frontera
    mask_fluid = (ptype == 0)
    mask_boundary = (ptype == 1)

    plt.scatter(x[mask_fluid], y[mask_fluid], s=5, c="blue", label="Fluido")
    plt.scatter(x[mask_boundary], y[mask_boundary], s=5, c="black", label="Frontera")

    # Vectores de velocidad (solo para fluido)
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
