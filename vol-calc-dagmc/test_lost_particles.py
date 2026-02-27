# Taken from /tests/unit_tests/dagmc/test_lost_particles.py of openmc-dev and
# modified for Volume Calculation

import numpy as np
from pathlib import Path

import openmc
import openmc.lib


openmc.reset_auto_ids()
model = openmc.Model()

### MATERIALS ###
fuel = openmc.Material(name='no-void fuel')
fuel.set_density('g/cc', 10.29769)
fuel.add_nuclide('U233', 1.0)

cladding = openmc.Material(name='clad')
cladding.set_density('g/cc', 6.55)
cladding.add_nuclide('Zr90', 1.0)

h1 = openmc.Material(name='water')
h1.set_density('g/cc', 0.75)
h1.add_nuclide('H1', 1.0)

model.materials = openmc.Materials([fuel, cladding, h1])

### GEOMETRY ###
# create the DAGMC universe using a model that has many triangles
# removed
dagmc_file = "broken_model.h5m"
pincell_univ = openmc.DAGMCUniverse(filename=dagmc_file, auto_geom_ids=True)

# create a 2 x 2 lattice using the DAGMC pincell
pitch = np.asarray((24.0, 24.0))
lattice = openmc.RectLattice()
lattice.pitch = pitch
lattice.universes = [[pincell_univ] * 2] * 2
lattice.lower_left = -pitch

# clip the DAGMC geometry at +/- 10 cm w/ CSG planes
rpp = openmc.model.RectangularParallelepiped(
    -pitch[0], pitch[0], -pitch[1], pitch[1], -10.0, 10.0, boundary_type='reflective')
bounding_cell = openmc.Cell(fill=lattice, region=-rpp)

model.geometry = openmc.Geometry(root=[bounding_cell])

# use 2x factor for tracing outside
ll = [-2.*pitch[0], -2.*pitch[1], -2.*10.]
ur = [ 2.*pitch[0],  2.*pitch[1],  2.*10.]

mat_vol    = openmc.VolumeCalculation(model.materials, 1000000, ll, ur)
mat_vol_rt = openmc.VolumeCalculation(model.materials, 1000000, ll, ur, 'ray')

# settings
settings = openmc.Settings()
settings.volume_calculations = [mat_vol, mat_vol_rt]
settings.run_mode = 'volume'

model.export_to_xml()
settings.export_to_xml()

# ensure that particles will be lost when cell intersections can't be found
# due to the removed triangles in this model
openmc.calculate_volumes()

# run this again, but with the dagmc universe as the root unvierse
# to ensure that lost particles are still caught in this case
for univ in model.geometry.get_all_universes().values():
    if isinstance(univ, openmc.DAGMCUniverse):
        model.geometry.root_universe = univ
        break

model.export_to_xml()
settings.export_to_xml()
openmc.calculate_volumes()

