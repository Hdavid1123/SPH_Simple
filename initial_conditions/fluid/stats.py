import numpy as np
import pandas as pd
from pathlib import Path

def count_part_file_col(points: np.ndarray, tol: float = 1e-8):
    """
    Cuenta cuántas partículas hay por columna (x) y por fila (y)
    después del filtrado.

    Args:
        points (np.ndarray): array Nx2 con coordenadas (x, y)
        tol (float): tolerancia para agrupar valores flotantes

    Returns:
        df_cols (pd.DataFrame): conteo por columna
        df_filas (pd.DataFrame): conteo por fila
    """
    # Normalizamos para evitar problemas de flotantes
    xs = np.round(points[:, 0] / tol) * tol
    ys = np.round(points[:, 1] / tol) * tol

    # Contar por columnas (agrupando por coordenada X)
    unique_x = np.unique(xs)
    col_counts = []
    for i, x in enumerate(sorted(unique_x), start=1):
        count = np.sum(xs == x)
        col_counts.append({"Col": i, "Num": count})

    # Contar por filas (agrupando por coordenada Y)
    unique_y = np.unique(ys)
    row_counts = []
    for i, y in enumerate(sorted(unique_y, reverse=True), start=1):  # de arriba hacia abajo
        count = np.sum(ys == y)
        row_counts.append({"Fil": i, "Num": count})

    import pandas as pd
    df_cols = pd.DataFrame(col_counts)
    df_filas = pd.DataFrame(row_counts)

    return df_cols, df_filas

def save_stats(points: np.ndarray, output_dir: Path):
    """
    Guarda las estadísticas de partículas de fluido en un archivo .txt
    dentro del directorio de salida.
    """
    df_cols, df_filas = count_part_file_col(points)

    stats_path = output_dir / "stats_fluid.txt"
    with open(stats_path, "w", encoding="utf-8") as f:
        f.write("stats_fluid.txt\n")
        f.write("Número de partículas por columna:\n")
        f.write(df_cols.to_string(index=False))
        f.write("\n\nNúmero de partículas por fila:\n")
        f.write(df_filas.to_string(index=False))

    print(f"[✓] Estadísticas de fluido guardadas en: {stats_path}")
