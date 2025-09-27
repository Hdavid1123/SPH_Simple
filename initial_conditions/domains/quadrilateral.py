import numpy as np
from typing import List
from initial_conditions.domains.base import BoundaryShape, Point, Segment
from initial_conditions.domains.utils import construir_trapecio, agregar_agujero

class Quadrilateral(BoundaryShape):
    def __init__(self,
                 d1: float, d2: float, d3: float,
                 a1: float, a2: float, a3: float,
                 resolution: int = 1,
                 holes: List[dict] = None,
                 sc_base: float = 1.0):
        """
        Construye un cuadrilátero normalizado con muestreo de lados,
        opcionalmente le aplica agujeros. Escala usando sc_base.
        """
        # Construir trapecio escalado
        lados, vertices = construir_trapecio(
            d1, d2, d3,
            a1, a2, a3,
            k=resolution,
            sc_base=sc_base
        )

        # Guardamos vértices con etiquetas A, B, C, D
        self._vertices: dict[str, np.ndarray] = vertices

        # Guardamos segmentos (cada lado tiene sus puntos muestreados)
        self._segments: dict[str, np.ndarray] = lados

        self._resolution: int = resolution

        # Agregar agujeros en lados si existen
        if holes:
            extremos = {
                "AB": ("A", "B"),
                "BC": ("B", "C"),
                "CD": ("C", "D"),
                "DA": ("D", "A")
            }
            for h in holes:
                side = h["lado"]
                if side in extremos:
                    v1, v2 = extremos[side]
                    P1, P2 = vertices[v1], vertices[v2]
                    
                    self._segments[side] = agregar_agujero(
                        P1, P2,
                        longitud=h["tam"],
                        offset=h["offset"],
                        k=resolution,
                        sc_base=sc_base
                    )

    # -------------------------
    # Métodos públicos
    # -------------------------
    def segments(self) -> List[Segment]:
        return list(self._segments.values())

    def vertices(self) -> List[Point]:
        return [tuple(self._vertices[k]) for k in ["A", "B", "C", "D"]]

    def vertex_dict(self) -> dict[str, Point]:
        return {k: tuple(v) for k, v in self._vertices.items()}
