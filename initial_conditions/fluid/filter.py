# fluid/filter.py
import numpy as np
from scipy.spatial import cKDTree

def eliminar_solapamientos(p_fluid, p_border, distancia_minima):
    """Elimina partículas de fluido que estén demasiado cerca de las fronteras."""
    if isinstance(p_border, list) and isinstance(p_border[0], dict):
        # Si recibimos partículas en formato dict, extraer coordenadas
        p_border = np.array([p["position"] for p in p_border])

    if len(p_border) == 0:
        return p_fluid

    tree_border = cKDTree(p_border)
    distancias, _ = tree_border.query(p_fluid, k=1)
    mascara = distancias > distancia_minima
    return p_fluid[mascara]
