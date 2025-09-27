# initial_conditions/main.py

import sys
from pathlib import Path
import argparse

# Asegurar que el proyecto raíz está en el PYTHONPATH
PROJECT_ROOT = Path(__file__).resolve().parents[0]
sys.path.insert(0, str(PROJECT_ROOT))

# Importar funciones disponibles
from boundaries.export import export_boundary_particles
from fluid.export import export_fluid_particles
from export_all import export_all_particles


def main():
    parser = argparse.ArgumentParser(
        description="Herramienta para exportar partículas de fronteras y fluido."
    )
    parser.add_argument(
        "command",
        choices=["export_boundaries", "export_fluid", "export_all"],
        help="Función a ejecutar"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Nombre del archivo de salida (por defecto depende del comando)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Visualizar el resultado gráficamente"
    )

    args = parser.parse_args()

    if args.command == "export_boundaries":
        filename = args.output or "boundary_particles.txt"
        export_boundary_particles(output_filename=filename, visualize=args.plot)

    elif args.command == "export_fluid":
        filename = args.output or "fluid_particles.txt"
        export_fluid_particles(output_filename=filename, visualize=args.plot)

    elif args.command == "export_all":
        filename = args.output or "all_particles.txt"
        export_all_particles(output_filename=filename, visualize=args.plot)


if __name__ == "__main__":
    main()
