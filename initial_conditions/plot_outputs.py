import pandas as pd
import matplotlib.pyplot as plt

def plot_particles(filename: str):
    # Leer archivo .txt, separador = tabulador
    df = pd.read_csv(filename, sep="\t")

    # Separar partículas por tipo
    boundary = df[df["type"] == 1]   # Frontera
    fluid = df[df["type"] == 0]      # Fluido

    # Crear figura
    plt.figure(figsize=(6, 6))

    # Graficar frontera (negro)
    plt.scatter(boundary["posx"], boundary["posy"],
                s=10, c="black", label="Frontera (type=1)")

    # Graficar fluido (azul)
    plt.scatter(fluid["posx"], fluid["posy"],
                s=10, c="blue", label="Fluido (type=0)")

    # Ajustes
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Visualización de partículas")
    plt.legend()
    plt.gca().set_aspect("equal", adjustable="box")  # Misma escala en x e y
    plt.grid(True)
    plt.show()

# Ejemplo de uso
if __name__ == "__main__":
    #plot_particles("outputs/problema_vasos_comunicantes.txt")
    plot_particles("outputs/salpicadura_vaciado_llenado.txt")
    #plot_particles("outputs/interaccion_frontera.txt")
