import matplotlib.pyplot as plt
import numpy as np
from initial_conditions.fluid.builder import FluidBuilder

def visualize_fluid(points: np.ndarray,
                    ax: plt.Axes = None,
                    show: bool = True) -> plt.Axes:
    own_ax = ax is None
    if own_ax:
        fig, ax = plt.subplots()

    ax.scatter(points[:, 0], points[:, 1], s=4, color='blue')
    ax.set_aspect('equal', 'box')
    ax.set_title("Partículas del fluido")

    if show and own_ax:
        plt.show()

    return ax

if __name__ == "__main__":
    from fluid.builder import FluidBuilder

    builder = FluidBuilder()
    points = builder.build()
    visualize_fluid(points)