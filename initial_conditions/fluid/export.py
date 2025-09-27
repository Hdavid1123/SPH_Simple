# initial_conditions/fluid/export.py

from pathlib import Path
from initial_conditions.fluid.builder import FluidBuilder
from initial_conditions.fluid.particleizer import FluidParticleizer
from initial_conditions.fluid.visualizer import visualize_fluid
from initial_conditions.fluid.stats import save_stats

def export_fluid_particles(output_filename: str = "fluid_particles.txt", visualize: bool = False):
    # 1. Construir puntos del fluido
    builder = FluidBuilder()
    points = builder.build()  # Ya incluye chequeo de fronteras si lo configuraste

    # 2. Obtener espaciado del archivo de parámetros
    espaciado = builder.config["espaciado"]

    # 3. Convertir puntos en partículas SPH
    particleizer = FluidParticleizer()
    particles = particleizer.generate(points, espaciado=espaciado)

    # 4. Ruta de salida en outputs/
    output_dir = Path(__file__).resolve().parents[1] / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    # 5. Guardar como TXT con todos los campos relevantes
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

    print(f"[✓] Archivo de fluido exportado en: {output_path}")

    # 6. Guardar estadísticas
    save_stats(points, output_dir)

    # 7. Visualizar si se solicita
    if visualize:
        visualize_fluid(points, show=True)
