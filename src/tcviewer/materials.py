import open3d as o3d

orbital_shiny = o3d.visualization.rendering.MaterialRecord()
orbital_shiny.shader = 'defaultLitTransparency'
orbital_shiny.base_color = [1, 1, 1, 0.9]
orbital_shiny.base_roughness = .8
orbital_shiny.base_reflectance = .3
orbital_shiny.base_clearcoat = 1
orbital_shiny.thickness = 1
orbital_shiny.transmission = 0

orbital_matte = o3d.visualization.rendering.MaterialRecord()
orbital_matte.shader = 'defaultLitTransparency'
orbital_matte.base_color = [1, 1, 1, 0.9]
orbital_matte.base_roughness = .8
orbital_matte.base_reflectance = .3
orbital_matte.base_clearcoat = 0
orbital_matte.thickness = 1
orbital_matte.transmission = 0

atom_shiny = o3d.visualization.rendering.MaterialRecord()
# atom_shiny.shader = 'defaultLit'
atom_shiny.base_color = [1, 1, 1, 1]
atom_shiny.base_roughness = .8
atom_shiny.base_reflectance = .3
atom_shiny.base_clearcoat = 1
atom_shiny.thickness = 1
atom_shiny.transmission = 0

atom_matte = o3d.visualization.rendering.MaterialRecord()
# atom_matte.shader = 'defaultLit'
atom_matte.base_color = [1, 1, 1, 1]
atom_matte.base_roughness = 1
atom_matte.base_reflectance = 0
atom_matte.base_clearcoat = 0
atom_matte.thickness = 1
atom_matte.transmission = 0


def orbital_material(material):

	if material == 'shiny':
		return orbital_shiny

	if material == 'matte':
		return orbital_matte

def atom_material(material):

	if material == 'shiny':
		return atom_shiny

	if material == 'matte':
		return atom_matte

