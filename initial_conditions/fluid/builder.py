# initial_conditions/fluid/builder.py

import json
from pathlib import Path
import numpy as np
from initial_conditions.fluid.geometry import sample_fluid_region
from initial_conditions.fluid.filter import eliminar_solapamientos

PARAM_PATH = Path(__file__).resolve().parent.parent / "parameters" / "fluid_region.json"


class FluidBuilder:
    def __init__(self, config_file: Path | str = PARAM_PATH):
        with open(config_file, 'r', encoding='utf-8') as f:
            self.config = json.load(f)

    def build(self, border_points: np.ndarray | None = None) -> np.ndarray:
        """
        Genera la nube de puntos de la región fluida.
        Si se proporcionan border_points, elimina solapamientos con las partículas de frontera.
        Además, imprime estadísticas para depuración.
        """
        # 1. Generar nube de puntos del fluido
        puntos = sample_fluid_region(self.config)
        print(f"[DEBUG] Puntos iniciales generados: {puntos.shape[0]}")

        # Guardar coordenadas únicas antes del filtrado
        y_vals_before = np.unique(np.round(puntos[:, 1], 6))
        x_vals_before = np.unique(np.round(puntos[:, 0], 6))

        # 2. Filtrar para evitar solapamiento con la frontera
        if border_points is not None and len(border_points) > 0:
            espaciado = self.config["espaciado"]
            puntos_filtrados = eliminar_solapamientos(puntos, border_points, espaciado / 2)
            print(f"[DEBUG] Puntos después del filtrado: {puntos_filtrados.shape[0]}")

            # Detectar filas perdidas
            y_vals_after = np.unique(np.round(puntos_filtrados[:, 1], 6))
            filas_perdidas = set(y_vals_before) - set(y_vals_after)
            if filas_perdidas:
                print(f"[DEBUG] Filas eliminadas (y): {sorted(filas_perdidas)}")

            # Detectar columnas perdidas
            x_vals_after = np.unique(np.round(puntos_filtrados[:, 0], 6))
            columnas_perdidas = set(x_vals_before) - set(x_vals_after)
            if columnas_perdidas:
                print(f"[DEBUG] Columnas eliminadas (x): {sorted(columnas_perdidas)}")

            puntos = puntos_filtrados

        return puntos