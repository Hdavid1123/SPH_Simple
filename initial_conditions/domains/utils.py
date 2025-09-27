import numpy as np

def segmentar_lado(P1, P2, sc):
    """
    Devuelve un arreglo de puntos entre P1 y P2.
    sc es la distancia entre puntos (ya normalizada globalmente).
    """
    d = np.linalg.norm(P2 - P1)
    n_puntos = max(2, int(np.ceil(d / sc)))
    return np.linspace(P1, P2, n_puntos)


def construir_trapecio(d1, d2, d3, a1, a2, a3, k=1, sc_base: float = 1.0):
    # Escalar longitudes
    d1, d2, d3 = d1 * sc_base, d2 * sc_base, d3 * sc_base
    angulos = np.deg2rad([a1, a2, a3])

    A = np.array([-0.5, 0.0])
    B = A + d1 * np.array([np.cos(angulos[0]), np.sin(angulos[0])])
    C = B + d2 * np.array([np.cos(angulos[1]), np.sin(angulos[1])])
    D = C + d3 * np.array([np.cos(angulos[2]), np.sin(angulos[2])])

    # Centrar verticalmente
    offset = - (min(A[1], B[1], C[1], D[1]) + max(A[1], B[1], C[1], D[1])) / 2
    for p in [A, B, C, D]:
        p[1] += offset

    # CORRECCIÓN: usar k*sc_base para segmentar lados
    lados = {
        'AB': segmentar_lado(A, B, k * sc_base),
        'BC': segmentar_lado(B, C, k * sc_base),
        'CD': segmentar_lado(C, D, k * sc_base),
        'DA': segmentar_lado(D, A, k * sc_base),
    }

    vertices = {'A': A, 'B': B, 'C': C, 'D': D}
    return lados, vertices


def agregar_agujero(P1, P2, longitud, offset, k=1, sc_base=1.0):
    longitud *= sc_base
    offset *= sc_base

    dir_vec = P2 - P1
    dir_unit = dir_vec / np.linalg.norm(dir_vec)
    inicio = P1 + offset * dir_unit

    puntos = segmentar_lado(P1, P2, k * sc_base)
    filtrados = np.array([p for p in puntos if not (0 <= np.dot(p - inicio, dir_unit) <= longitud)])
    return filtrados
