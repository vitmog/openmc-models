# "Spheres within spheres" model for volume calculation testing

import openmc
from math import pi

fuel = openmc.Material(1)
fuel.add_element('U', 1, enrichment=3.0)
fuel.add_element('O', 2)
fuel.set_density('g/cm3', 10.0)

water = openmc.Material(2)
water.add_nuclide('H1', 2)
water.add_nuclide('O16', 1)
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('g/cm3', 1)

wapor = openmc.Material(3)
wapor.add_nuclide('H1', 2)
wapor.add_nuclide('O16', 1)
wapor.add_s_alpha_beta('c_H_in_H2O')
wapor.set_density('g/cm3', 1)

carbon = openmc.Material(4)
carbon.add_nuclide('C12', 1)
carbon.set_density('g/cm3', 1.8)

air = openmc.Material(5)
air.add_nuclide('N14', 0.78, 'ao')
air.add_nuclide('O16', 0.22, 'ao')
air.set_density('g/cm3', 1.22e-3)

mats = openmc.Materials([fuel, water, wapor, carbon, air])
mats.export_to_xml()

s1 = openmc.Sphere(r=10., boundary_type='vacuum')
s2 = openmc.Sphere(r=1.0001)
s3 = openmc.Sphere(r=1.)
s4 = openmc.Sphere(r=.1)            # sphere in the center
s5 = openmc.Sphere(x0=0.2, r=.01)   # smaller sphere near the center
s6 = openmc.Sphere(x0=9.8, r=.1)    # sphere on periphery

# Compute real volumes of cells
vol = (
    4 * pi / 3 * pow(s5.r,3),
    4 * pi / 3 * pow(s4.r,3),
    4 * pi / 3 * (pow(s3.r,3) - pow(s4.r,3) - pow(s5.r,3)),
    4 * pi / 3 * (pow(s2.r,3) - pow(s3.r,3)),
    4 * pi / 3 * (pow(s1.r,3) - pow(s2.r,3) - pow(s6.r,3)),
    4 * pi / 3 * (pow(s6.r,3)))

carb_sphere   = openmc.Cell(1, fill=carbon, region=-s5, 
                            name = "{:10.4e}".format(vol[0]) + ' | '
                            'R < 0.01 cm sphere')
carb_lsphere  = openmc.Cell(2, fill=carbon, region=-s4,                             
                            name = "{:10.4e}".format(vol[1]) + ' | '
                            'R < 0.1 cm sphere in the center')
inner_sphere  = openmc.Cell(3, fill=fuel,   region=+s4 & -s3 & +s5,                             
                            name = "{:10.4e}".format(vol[2]) + ' | '
                            '0.1 < R < 1. cm layer - sm.sph.')
gap_layer     = openmc.Cell(4, fill=air,    region=+s3 & -s2,                             
                            name = "{:10.4e}".format(vol[3]) + ' | '
                            '1.0 < R < 1.0001 cm layer')
outer_sphere  = openmc.Cell(5, fill=water,  region=+s2 & +s6 & -s1,                             
                            name = "{:10.4e}".format(vol[4]) + ' | '
                            '1.0001 < R < 10 cm layer')
bubble_sphere = openmc.Cell(6, fill=wapor,  region=-s6,                             
                            name = "{:10.4e}".format(vol[5]) + ' | '
                            'R < 0.1 cm sphere on periphery')

cells_list = (carb_sphere, carb_lsphere, inner_sphere, gap_layer, 
              outer_sphere, bubble_sphere)

root = openmc.Universe(0, cells=(
    carb_sphere, carb_lsphere, inner_sphere, gap_layer, outer_sphere, 
    bubble_sphere))

geom = openmc.Geometry(root)
geom.export_to_xml()

lower_left, upper_right = geom.bounding_box
cell_vol = openmc.VolumeCalculation(
    cells_list, 100000000, lower_left, upper_right
    )
cell_vol_rt = openmc.VolumeCalculation(
    cells_list, 30000000, lower_left, upper_right, 'ray'
    )
mat_vol = openmc.VolumeCalculation(
    [fuel, water, wapor, carbon, air],
    100000000, lower_left, upper_right
    )
mat_vol_rt = openmc.VolumeCalculation(
    [fuel, water, wapor, carbon, air],
    30000000, lower_left, upper_right, 'ray'
    )
univ_vol = openmc.VolumeCalculation(
    [root], 100000000, lower_left, upper_right
    )
univ_vol_rt = openmc.VolumeCalculation(
    [root], 30000000, lower_left, upper_right, 'ray'
    )

settings = openmc.Settings()
settings.volume_calculations = [
    cell_vol, mat_vol, univ_vol, cell_vol_rt, mat_vol_rt, univ_vol_rt
    ]
settings.run_mode = 'volume'
settings.export_to_xml()

openmc.calculate_volumes()
