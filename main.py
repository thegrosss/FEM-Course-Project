import random
from mesh.mesh_parameters import MeshParameters
from mesh.mesh_builder import MeshBuilder
from mesh.point import Point
from fem.fem_solver import FemSolver

parameters = MeshParameters.read_json("input/area.json")

mesh_builder = MeshBuilder(parameters)

mesh_builder.create_points()
mesh_builder.create_elements()
mesh_builder.create_boundaries()
mesh = mesh_builder.get_mesh()
solver = FemSolver(mesh)
solver.solve()

# Дальше будем сравнивать аналитическое решение и численное
# Будем сравнивать в наборе произвольных точек
def generate_random_points(r_min: float, r_max: float, z_min: float, z_max: float, points_count: int):
    points_list = []
    for i in range(points_count):
        r = r_min + random.random() * (r_max - r_min)
        z = z_min + random.random() * (z_max - z_min)
        points_list.append(Point(r, z))

    with open("output/points_random", "w") as file:
        for p in points_list:
            file.write(f"{p.r} {p.z}\n")
    return points_list

r_min = min(p.r for p in parameters.control_points)
r_max = max(p.r for p in parameters.control_points)
z_min = min(p.z for p in parameters.control_points)
z_max = max(p.z for p in parameters.control_points)
#points = generate_random_points(r_min, r_max, z_min, z_max, 20)

def read_points():
    points = []
    with open("output/points_random", "r") as file:
        lines = file.readlines()
        for l in lines:
            x, y = map(float, l.split())
            points.append(Point(x, y))
    return points

points = read_points()

print(f"Погрешность: {solver.compare_solution_with_exact_in_nodes():.2e}")
residual = solver.root_mean_square(points)

print(f"Невязка в произвольных точках: {residual:.2e}")