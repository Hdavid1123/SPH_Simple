from abc import ABC, abstractmethod
from typing import List, Tuple

Point = Tuple[float, float]
Segment = List[Point]


class BoundaryShape(ABC):
    """Interfaz base para cualquier geometría 2D representada por segmentos."""

    @abstractmethod
    def segments(self) -> List[Segment]:
        """Retorna lista de segmentos (lados, líneas extra, etc.)"""
        pass
