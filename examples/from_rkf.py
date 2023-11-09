from yutility import orbitals
import tcviewer
import pathlib

rkf_dir = pathlib.Path(__file__).parents[1]/'test'/'fixtures'/'NH3BH3'

# load orbitals and choose an MO to draw
orbs = orbitals.Orbitals(rkf_dir/'adf.rkf')
homo = orbs.mos['HOMO']
# generate a cube file
cub = homo.generate_orbital()

with tcviewer.Screen() as scr:
	scr.draw_cub(cub, 0.03, material=tcviewer.materials.orbital_matte)
	