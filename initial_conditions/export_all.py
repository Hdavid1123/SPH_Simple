# initial_conditions/export_all.py
from pathlib import Path
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt

from initial_conditions.boundaries.builder import BoundaryBuilder
from initial_conditions.fluid.builder import FluidBuilder
from initial_conditions.fluid.particleizer import FluidParticleizer
from initial_conditions.fluid.stats import save_stats

OUTPUT_DIR_DEFAULT = Path(__file__).resolve().parent / "outputs"

# ----------------------------
# Construcción / Particleización
# ----------------------------
def build_boundary_particles(resolution: int = None,
                             particle_type: int = 1,
                             h: float = None,
                             dx: float = None,
                             dy: float = None) -> List[dict]:
    """Construye y devuelve lista de partículas de frontera (list[dict])."""
    builder = BoundaryBuilder()
    particles = builder.build(resolution=resolution,
                              particle_type=particle_type,
                              h=h, dx=dx, dy=dy)
    return particles


def build_fluid_points(border_points: np.ndarray | None = None) -> Tuple[np.ndarray, float]:
    """
    Construye puntos del fluido (np.ndarray) y devuelve también el espaciado leído del JSON.
    Retorna: (points, espaciado)
    """
    builder = FluidBuilder()
    points = builder.build(border_points=border_points)
    espaciado = builder.config.get("espaciado")
    if espaciado is None:
        raise ValueError("El parámetro 'espaciado' debe estar definido en fluid_region.json")
    return points, float(espaciado)


def particleize_fluid(points: np.ndarray, espaciado: float,
                      ptype: int = 0, h: float | None = None) -> List[dict]:
    """Convierte puntos del fluido en partículas (list[dict]) usando FluidParticleizer."""
    particleizer = FluidParticleizer()
    particles = particleizer.generate(points, espaciado=espaciado, ptype=ptype, h=h)
    return particles


# ----------------------------
# Export / utilidades
# ----------------------------
def assign_ids(all_particles: List[dict]) -> None:
    """Asigna ids consecutivos en el campo 'id' de cada dict (in-place)."""
    for i, p in enumerate(all_particles):
        p["id"] = i


def save_particles_txt(all_particles: List[dict], output_path: Path) -> None:
    """
    Guarda en texto con columnas: id posx posy h type mass dx dy
    Sobrescribe output_path si existe.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("id\tposx\tposy\th\ttype\tmass\tdx\tdy\n")
        for i, p in enumerate(all_particles):
            # Asegurarse de que cada partícula tenga los campos esperados
            pid = p.get("id", i)
            x, y = p["position"]
            h = p["h"]
            ptype = p.get("type", 0)
            mass = p.get("mass", 0.0)
            dx = p.get("dx", 0.0)
            dy = p.get("dy", 0.0)
            f.write(f"{pid}\t{x:.6f}\t{y:.6f}\t{h:.6f}\t{ptype}\t{mass:.6f}\t{dx:.6f}\t{dy:.6f}\n")


def make_numpy_for_cpp(all_particles: List[dict]) -> np.ndarray:
    """
    Construye un numpy.ndarray shape (N,8) con columnas:
    id, posx, posy, h, type, mass, dx, dy
    """
    n = len(all_particles)
    data = np.zeros((n, 8), dtype=np.float64)
    for i, p in enumerate(all_particles):
        pid = p.get("id", i)
        x, y = p["position"]
        h = p["h"]
        ptype = p.get("type", 0)
        mass = p.get("mass", 0.0)
        dx = p.get("dx", 0.0)
        dy = p.get("dy", 0.0)
        data[i] = [pid, x, y, h, ptype, mass, dx, dy]
    return data


def visualize_particles(boundary_particles: List[dict],
                        fluid_particles: List[dict]) -> None:
    """Dibuja frontera y fluido en el mismo plot."""
    fig, ax = plt.subplots()

    if boundary_particles:
        bx, by = zip(*[p["position"] for p in boundary_particles])
        ax.scatter(bx, by, s=6, color="black", label="Frontera")
    if fluid_particles:
        fx, fy = zip(*[p["position"] for p in fluid_particles])
        ax.scatter(fx, fy, s=6, color="blue", label="Fluido")

    ax.set_title("Partículas de frontera y fluido")
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale()
    ax.legend()
    plt.show()


# ----------------------------
# Función orquestadora (SRP: coordina, no hace todo internamente)
# ----------------------------
def export_all_particles(output_filename: str = "all_particles.txt",
                         output_dir: Path | None = None,
                         visualize: bool = False) -> Tuple[List[dict], np.ndarray]:
    """
    Orquesta la generación, particleización y exportación:
      - construye fronteras
      - construye puntos de fluido y los particleiza evitando solapamientos
      - guarda archivo .txt (id,posx,posy,h,type,mass,dx,dy)
      - guarda estadísticas de fluido
      - devuelve (all_particles_list, numpy_array_for_cpp)
    """
    output_dir = (Path(output_dir) if output_dir is not None else OUTPUT_DIR_DEFAULT)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    # 1) Fronteras
    boundary_particles = build_boundary_particles()

    # 2) Puntos del fluido (filtrados contra frontera)
    border_points = np.array([p["position"] for p in boundary_particles]) if boundary_particles else None
    points_fluid, espaciado = build_fluid_points(border_points=border_points)

    # 3) Particleizar fluido
    fluid_particles = particleize_fluid(points_fluid, espaciado=espaciado, ptype=0)

    # 4) Combinar
    all_particles = boundary_particles + fluid_particles

    # 5) Asignar ids consecutivos y exportar .txt
    assign_ids(all_particles)
    save_particles_txt(all_particles, output_path)

    # 6) Guardar estadísticas del fluido (usa puntos antes de particleizar)
    save_stats(points_fluid, output_dir)

    # 7) Construir numpy array listo para C++
    data_for_cpp = make_numpy_for_cpp(all_particles)

    # 8) Visualizar (opcional)
    if visualize:
        visualize_particles(boundary_particles, fluid_particles)

    return all_particles, data_for_cpp