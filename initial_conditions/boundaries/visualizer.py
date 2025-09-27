# boundaries/visualizer.py
import matplotlib.pyplot as plt
from typing import List, Tuple, Any

Point = Tuple[float, float]
Segment = List[Point]


def visualize_boundary(particles: List[dict],
                       domain: Any = None,
                       show: bool = True,
                       save_path: str = None,
                       ax: plt.Axes = None) -> plt.Axes:
    """
    Dibuja la nube de partículas SPH como representación de la frontera.

    Args:
        particles: lista de partículas, cada una con campos incluyendo 'position'.
        domain: objeto de dominio (Quadrilateral, CompositeDomain, etc.)
                que implemente .endpoints().
        show: si True muestra el plot al final (plt.show()).
        save_path: si se especifica, guarda la imagen en esa ruta.
        ax: eje existente donde dibujar (opcional).

    Returns:
        plt.Axes: el eje donde se dibujó.
    """
    own_ax = ax is None
    if own_ax:
        fig, ax = plt.subplots()

    # Partículas en negro
    xs = [p["position"][0] for p in particles]
    ys = [p["position"][1] for p in particles]
    ax.plot(xs, ys, "ko", markersize=2)

    # Extremos de las líneas en rojo con índices
    if domain is not None:
        extremos = domain.endpoints()
        for idx, (x, y) in enumerate(extremos):
            ax.plot(x, y, "ro")  # marcar extremo
            ax.text(x, y, str(idx),
                    color="red", fontsize=8,
                    ha="center", va="center",
                    bbox=dict(facecolor="white", alpha=0.7, boxstyle="round,pad=0.2"))

    ax.set_aspect("equal", "box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Frontera SPH")

    if save_path:
        plt.savefig(save_path, dpi=150)
    elif show and own_ax:
        plt.show()

    return ax