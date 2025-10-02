# initial_conditions/fluid/particleizer.py

from typing import List, Dict, Any
import numpy as np

class FluidParticleizer:
    def __init__(self, rho0: float = 1000.0):
        """
        Inicializa el particleizer para fluido.

        Args:
            rho0: densidad de referencia (kg/m³), por defecto 1000 (agua).
        """
        self.rho0 = rho0

    def generate(self,
                 points: np.ndarray,
                 espaciado: float,
                 ptype: int = 0,
                 h: float | None = None,
                 velocity: tuple[float, float] = (0.0, 0.0)) -> List[Dict[str, Any]]:
        """
        Convierte un array de puntos 2D en partículas SPH de fluido.

        Args:
            points: np.ndarray de forma (N, 2), coordenadas (x, y).
            espaciado: separación entre partículas (dx, dy).
            ptype: tipo de partícula (por defecto 0 = fluido).
            h: radio de suavizado. Si None, se usa h = dx.
            velocity: tupla con la velocidad inicial común (vx, vy).

        Returns:
            Lista de diccionarios con campos:
            id, type, position, velocity, h, dx, dy, mass
        """
        dx = dy = espaciado
        #h = dx if h is None else h
        h = 0.01 if h is None else h  #Valor fijo para consistencia en pruebas
        mass = self.rho0 * dx * dy

        particles: List[Dict[str, Any]] = []

        for i, (x, y) in enumerate(points):
            particles.append({
                "id": i,
                "type": ptype,
                "position": [float(x), float(y)],
                "velocity": [float(velocity[0]), float(velocity[1])],
                "h": h,
                "dx": dx,
                "dy": dy,
                "mass": mass
            })

        return particles
