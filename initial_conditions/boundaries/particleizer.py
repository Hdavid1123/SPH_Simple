# boundaries/particleizer.py

from typing import List, Dict, Any, Tuple

Point = Tuple[float, float]
Segment = List[Point]

class BoundaryParticleizer:
    def __init__(self):
        # Podrías parametrizar aquí settings globales si los necesitas
        pass

    def generate(self,
                 segments: List[Segment],
                 ptype: int = 1,
                 h: float = 0.01,
                 dx: float = 0.01,
                 dy: float = 0.01,
                 mass: float = 1.0) -> List[Dict[str, Any]]:
        """
        Convierte lista de segmentos (líneas de frontera) en partículas SPH.

        Args:
            segments: Cada segmento es una lista de puntos (x, y).
            ptype: entero que identifica el tipo de partículas.
            h: radio de suavizado.
            dx, dy: espaciamientos de partículas.
            mass: masa de cada partícula.

        Returns:
            List[Dict]: cada dict contiene keys: id, type, position, velocity, h, dx, dy, mass
        """
        particles: List[Dict[str, Any]] = []
        seen: set[Tuple[float, float]] = set()
        id_counter = 0

        for seg in segments:
            for x, y in seg:
                key = (round(x, 8), round(y, 8))  # Precisión configurable
                if key in seen:
                    continue
                seen.add(key)
                particles.append({
                    "id": id_counter,
                    "type": ptype,
                    "position": [x, y],
                    "velocity": [0.0, 0.0],
                    "h": h,
                    "dx": dx,
                    "dy": dy,
                    "mass": mass
                })
                id_counter += 1

        return particles
