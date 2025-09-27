# initial_conditions/boundaries/export.py

from pathlib import Path
from initial_conditions.boundaries.builder import BoundaryBuilder
from initial_conditions.boundaries.visualizer import visualize_boundary


def export_boundary_particles(output_filename: str = "boundary_particles.txt", visualize: bool = False):
    # 1. Construir partículas de frontera
    builder = BoundaryBuilder()
    particles = builder.build()  # resolution y escala vienen del JSON y builder

    # 2. Ruta de salida: dentro de outputs/
    output_dir = Path(__file__).resolve().parents[1] / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    # 3. Guardar archivo de texto con todos los campos relevantes
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("posx\tposy\th\ttype\tmass\tdx\tdy\n")
        for p in particles:
            x, y = p["position"]
            h = p["h"]
            ptype = p["type"]
            mass = p["mass"]
            dx = p["dx"]
            dy = p["dy"]
            f.write(f"{x:.6f}\t{y:.6f}\t{h:.6f}\t{ptype}\t{mass:.6f}\t{dx:.6f}\t{dy:.6f}\n")

    print(f"[✓] Archivo de frontera exportado en: {output_path}")

    # 4. Visualizar si se solicita
    if visualize:
        visualize_boundary(particles, show=True)
