# Taken from https://docs.openmc.org/en/v0.11.0/examples/cad-geom.html#using-cad-based-geometries
# and modified for Volume Calculation

import numpy as np
from pathlib import Path

import openmc
import openmc.lib

import urllib.request
import openmc.dagmc

teapot_url = 'https://tinyurl.com/y4mcmc3u' # 29 MB

def download(url):
    """
    Helper function for retrieving dagmc models
    """
    u = urllib.request.urlopen(url)

    if u.status != 200:
        raise RuntimeError("Failed to download file.")

    # save file as dagmc.h5m
    with open("dagmc.h5m", 'wb') as f:
        f.write(u.read())


openmc.reset_auto_ids()
model = openmc.Model()

### MATERIALS ###
water = openmc.Material(name="water")
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.set_density('g/cc', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.id = 41

iron = openmc.Material(name="iron")
iron.add_nuclide("Fe54", 0.0564555822608)
iron.add_nuclide("Fe56", 0.919015287728)
iron.add_nuclide("Fe57", 0.0216036861685)
iron.add_nuclide("Fe58", 0.00292544384231)
iron.set_density("g/cm3", 7.874)

air_csg = openmc.Material(name="air")
air_csg.add_nuclide("Fe54", 0.0564555822608)
air_csg.add_nuclide("Fe56", 0.919015287728)
air_csg.set_density("g/cm3", 7.874)

model.materials = openmc.Materials([iron, water, air_csg])

### GEOMETRY ###
download(teapot_url)
teapot_univ = openmc.DAGMCUniverse(filename="dagmc.h5m", auto_geom_ids=True)

inner_sphere = openmc.Sphere(r=40.0)
outer_sphere = openmc.Sphere(r=41.0, boundary_type='vacuum')

dagmc_cell = openmc.Cell(fill=teapot_univ, region=-inner_sphere)
outer_cell  = openmc.Cell(fill=air_csg, region=+inner_sphere & -outer_sphere)

model.geometry = openmc.Geometry([dagmc_cell, outer_cell])

bb = model.geometry.bounding_box
margin = .01

import numpy as np
ll = list(np.clip(bb.lower_left  + margin, -inner_sphere.r + margin,  0))
ur = list(np.clip(bb.upper_right - margin,  0, inner_sphere.r - margin))

mat_vol    = openmc.VolumeCalculation(model.materials, 300000, ll, ur)
mat_vol_rt = openmc.VolumeCalculation(model.materials, 100000, ll, ur, 'ray')

# settings
settings = openmc.Settings()
settings.dagmc = True
settings.volume_calculations = [mat_vol, mat_vol_rt]
settings.run_mode = 'volume'
model.settings = settings

model.export_to_xml()
settings.export_to_xml()

openmc.calculate_volumes()

