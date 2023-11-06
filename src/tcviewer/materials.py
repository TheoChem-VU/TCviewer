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
