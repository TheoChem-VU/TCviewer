from yutility import orbitals
import tcviewer
import pathlib

rkf_dir = pathlib.Path(__file__).parents[0]/'data'/'NH3BH3'

# load orbitals and choose an MO to draw
orbs = orbitals.Orbitals(rkf_dir/'adf.rkf')
homo = orbs.mos['HOMO']
# generate a cube file
cub = homo.generate_orbital()

with tcviewer.Screen() as scr:
	scr.draw_cub(cub, 0.03, material=tcviewer.materials.orbital_shiny)
	scr.draw_cub(cub, 0, material=tcviewer.materials.orbital_matte, color1=[1,1,1], color2=[1,1,1])
