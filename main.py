import random
from pathlib import Path

from fem.fem_solver import FemSolver
from mesh.mesh_builder import MeshBuilder
from mesh.mesh_parameters import MeshParameters
from mesh.point import Point
from mesh.time_parameters import TimeParameters


def generate_random_points(r_min: float, r_max: float, z_min: float, z_max: float, points_count: int):
    points_list = []
    for _ in range(points_count):
        r = r_min + random.random() * (r_max - r_min)
        z = z_min + random.random() * (z_max - z_min)
        points_list.append(Point(r, z))

    with open("output/points_random", "w") as file:
        for p in points_list:
            file.write(f"{p.r} {p.z}\n")
    return points_list


def read_points():
    points = []
    path = Path("output/points_random")

    if not path.exists():
        return points

    with open(path, "r") as file:
        lines = file.readlines()
        for line in lines:
            x, y = map(float, line.split())
            points.append(Point(x, y))
    return points


parameters = MeshParameters.read_json("input/area.json")

mesh_builder = MeshBuilder(parameters)
mesh_builder.create_points()
mesh_builder.create_elements()
mesh_builder.create_boundaries()

mesh = mesh_builder.get_mesh()
solver = FemSolver(mesh)

time_config = Path("input/time.json")
if time_config.exists():
    time_parameters = TimeParameters.read_json(str(time_config))
    solver.solve_time(time_parameters)
else:
    solver.solve()

r_min = min(p.r for p in parameters.control_points)
r_max = max(p.r for p in parameters.control_points)
z_min = min(p.z for p in parameters.control_points)
z_max = max(p.z for p in parameters.control_points)

points = read_points()
if not points:
    points = generate_random_points(r_min, r_max, z_min, z_max, 20)

print(f"Node relative error: {solver.compare_solution_with_exact_in_nodes():.2e}")
residual = solver.root_mean_square(points)

print(f"RMS error in arbitrary points: {residual:.2e}")
