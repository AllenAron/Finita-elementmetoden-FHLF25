import matplotlib as mpl
import numpy as np

from model import model
from solver import solver
from visualizer import visualizer

mpl.use('TkAgg')            

# Parameters
N = 2001 # number of timesteps
t0 = 0  # initial time
t1 = 2001   # end time
elType = 2  # type of meshing element
elSizeFactor = 0.01 # size of mesh

# Create the geometry and mesh
geometry = model.geometry()
temperature_mesh = model.create_mesh(geometry=geometry, elType=elType, dofsPerNode=1, elSizeFactor=elSizeFactor) # mesh for heat equation
stress_mesh = model.create_mesh(geometry=geometry, elType=elType, dofsPerNode=2, elSizeFactor=elSizeFactor) # mesh for von mises problem

# Solver
temperature_solver = solver(temperature_mesh) # create solver object for the heat equation
a_stationary = temperature_solver.stationary_solver() # compute solution for stationary heat equation
a_convection = temperature_solver.convection_solver()   # compute solution for convection heat equation
A_transient, tgrid, tmax, tmin = temperature_solver.transient_solver(t0, t1, N)   # compute solution for transient heat equation

stress_solver = solver(stress_mesh) # create solver object for the von mises problem
vonMises, displacements = stress_solver.von_mises_solver(A_transient[tmax])    # compute solution for the von mises problem

# Visualizer
vis = visualizer(geometry)  # create visualization object from the geometry
vis.draw_mesh(temperature_mesh)

# Visualize the solutions for the heat equation
vis.visualize_temperature_solution(data=a_stationary, mesh=temperature_mesh, vmin=np.min(a_stationary), vmax=np.max(a_stationary), draw_mesh=False, draw_geometry=False, cmap='turbo')
vis.visualize_temperature_solution(data=a_convection, mesh=temperature_mesh, vmin=np.min(a_convection), vmax=np.max(a_convection), draw_mesh=False, draw_geometry=False, cmap='turbo')
values = A_transient[tmax] # chose which time-frame to display
vis.visualize_temperature_solution(data=values, mesh=temperature_mesh, vmin=np.min(values), vmax=np.max(values), draw_mesh=False, draw_geometry=False, cmap='turbo')

# Visualize the solution for the von mises problem
vis.visualize_stress_solution(data=vonMises, disp=displacements, mesh=stress_mesh, dofsPerNode=stress_mesh.dofs_per_node, elType=elType, draw_values=True, draw_displacements=True)
vis.show()

# plot the transient temperature plot
vis.plotminmax(A_transient, tgrid, N)