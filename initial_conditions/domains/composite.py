import numpy as np
from typing import List, Optional
from initial_conditions.domains.base import Point, Segment, BoundaryShape


class CompositeDomain(BoundaryShape):
    """
    Representa un dominio compuesto por varias formas 2D
    (cuadriláteros, círculos, etc.) y conexiones (líneas entre extremos).
    """

    def __init__(self,
                 domains: List[BoundaryShape] = None,
                 connections: List[Segment] = None,
                 sc_base: float = 1.0):
        self.domains = domains or []
        self.connections = connections or []
        self.named_points = {}  # almacena etiquetas -> coordenadas
        self.sc_base = sc_base

    # -------------------------------
    # Registro y resolución de puntos
    # -------------------------------
    def add_domain(self, domain: BoundaryShape):
        """Agrega una forma al dominio compuesto."""
        self.domains.append(domain)

    def register_point(self, name: str, coords):
        """Asocia un nombre (ej: 'A', 'B', 'P1') a coordenadas."""
        self.named_points[name] = np.array(coords, dtype=float)

    def resolve_point(self, point):
        """Convierte un nombre o coordenadas en un np.array."""
        if isinstance(point, str):
            if point in self.named_points:
                return self.named_points[point]
            raise ValueError(f"Punto '{point}' no está registrado")
        return np.array(point, dtype=float)

    # -------------------------------
    # Construcción de conexiones
    # -------------------------------
    def add_connection(self, p1: Point, p2: Point, resolution: int = 1):
        """Añade una línea de conexión entre dos puntos existentes."""
        p1, p2 = self.resolve_point(p1), self.resolve_point(p2)

        # coherencia con Quadrilateral: separación = resolution * sc_base
        dist_punto = resolution * self.sc_base
        d = np.linalg.norm(p2 - p1)
        n_points = max(2, int(np.ceil(d / dist_punto)))
        seg = np.linspace(p1, p2, n_points)
        self.connections.append(seg)

    def add_free_line(self,
                      from_point: Point,
                      to_point: Optional[Point] = None,
                      length: Optional[float] = None,
                      angle: Optional[float] = None,
                      resolution: int = 1,
                      name: Optional[str] = None):
        """
        Añade una línea que parte desde un punto y termina en:
        - un punto arbitrario (to_point), o
        - un punto definido por (length, angle) relativo a from_point.

        El punto final se registra automáticamente en named_points.
        """
        p1 = self.resolve_point(from_point)

        if to_point is not None:
            p2 = self.resolve_point(to_point)
        elif length is not None and angle is not None:
            theta = np.radians(angle)
            dx, dy = length * np.cos(theta), length * np.sin(theta)
            p2 = p1 + np.array([dx, dy])
        else:
            raise ValueError("Debe especificar `to_point` o (`length` y `angle`).")

        # coherencia con Quadrilateral: separación = resolution * sc_base
        dist_punto = resolution * self.sc_base
        d = np.linalg.norm(p2 - p1)
        n_points = max(2, int(np.ceil(d / dist_punto)))
        seg = np.linspace(p1, p2, n_points)
        self.connections.append(seg)

        if name is None:
            name = f"P{len(self.named_points)}"
        self.register_point(name, p2)

        return name

    # -------------------------------
    # Consultas
    # -------------------------------
    def segments(self) -> List[Segment]:
        """Devuelve todos los segmentos: internos + conexiones."""
        segs = []
        for d in self.domains:
            segs.extend(d.segments())
        segs.extend(self.connections)
        return segs
    
    def print_points(self):
        """Imprime todos los puntos registrados (vértices y finales de líneas)."""
        print("=== Puntos registrados en CompositeDomain ===")
        for name, coords in self.named_points.items():
            print(f"{name}: {coords}")
        print("============================================")
