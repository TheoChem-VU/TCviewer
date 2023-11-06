import open3d as o3d
from yutility import atom_data, geometry
from scm import plams
import numpy as np
from tcintegral import grid
import skimage
from tcviewer import materials




class Screen:
    def __init__(self, **kwargs):
        self.screenshot_file = kwargs.get('screenshot_file')
        self.meshes = []

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        if self.screenshot_file is None:
            o3d.visualization.draw(self.meshes, show_skybox=False, title='TCViewer')
            return

        vis = o3d.visualization.Visualizer()
        vis.create_window()
        [vis.add_geometry(mesh) for mesh in self.meshes]
        # vis.update_geometry()
        # vis.update_renderer()
        vis.capture_screen_image(self.screenshot_file, True)
        vis.destroy_window()




    def add_mesh(self, geometry, name=None, material=None):
        self.meshes.append(dict(geometry=geometry, name=name or str(id(geometry)), material=material))

    def draw_molecule(self, mol: str or plams.Molecule, **kwargs):
        if isinstance(mol, str):
            mol = plams.Molecule(mol)
            mol.guess_bonds()

        for atom in mol:
            sphere = o3d.geometry.TriangleMesh.create_sphere(atom_data.radius(atom.symbol)*kwargs.get('atom_scale', .5), resolution=kwargs.get('atom_resolution', 100))
            sphere.translate(atom.coords)
            sphere = sphere.paint_uniform_color(np.array(atom_data.color(atom.symbol))/255)
            sphere.compute_vertex_normals()

            self.add_mesh(sphere, material=kwargs.get('atom_material'))

        for bond in mol.bonds:
            length = bond.length()
            cylinder = o3d.geometry.TriangleMesh.create_cylinder(kwargs.get('bond_width', .07), length)
            cylinder.translate((np.array(bond.atom1.coords) + np.array(bond.atom2.coords))/2)

            bond_vec = np.array(bond.atom1.coords) - np.array(bond.atom2.coords)
            R = geometry.vector_align_rotmat([0, 0, 1], bond_vec)
            cylinder.rotate(R)
            cylinder.compute_vertex_normals()
            cylinder.paint_uniform_color((0, 0, 0))

            self.add_mesh(cylinder, material=kwargs.get('bond_material'))

    def draw_isosurface(self, gridd, thresh=0, color=None, material=None):
        val = gridd.values
        s = gridd.shape
        val = val.reshape(*s)
        verts, faces, normals, values = skimage.measure.marching_cubes(val, thresh, spacing=gridd.spacing*3, method='lorensen')
        verts = o3d.cpu.pybind.utility.Vector3dVector(verts + gridd.origin)
        triangles = o3d.cpu.pybind.utility.Vector3iVector(faces)

        mesh = o3d.geometry.TriangleMesh(verts, triangles)
        if color:
            mesh.paint_uniform_color(color)

        mesh.compute_vertex_normals()
        self.add_mesh(mesh, material=material)

    def draw_orbital(self, orb, gridd=None, isovalue=0.03, color1=[1, 0, 0], color2=[0, 0, 1], material=materials.orbital_shiny):
        if gridd is None:
            gridd = grid.molecule_bounding_box(orb.molecule, spacing=.1)

        gridd.values = orb(gridd.points)
        if gridd.values.min() < isovalue < gridd.values.max():
            self.draw_isosurface(gridd, isovalue, color=color1, material=material)
        gridd.values = orb(gridd.points)
        if gridd.values.min() < -isovalue < gridd.values.max():
            self.draw_isosurface(gridd, -isovalue, color=color2, material=material)

    def draw_axes(self, center=[0, 0, 0], length=1, width=.04, **kwargs):
        arrow_x = o3d.geometry.TriangleMesh.create_arrow(width, width*2, cylinder_height=length, cone_height=length*.2, **kwargs)
        arrow_x.paint_uniform_color([1, 0, 0])
        arrow_x.translate(-np.array(center))
        arrow_x.rotate(geometry.vector_align_rotmat([0, 0, 1], [1, 0, 0]), -np.array(center))

        arrow_y = o3d.geometry.TriangleMesh.create_arrow(width, width*2, cylinder_height=length, cone_height=length*.2, **kwargs)
        arrow_y.paint_uniform_color([0, 1, 0])
        arrow_y.translate(-np.array(center))
        arrow_y.rotate(geometry.vector_align_rotmat([0, 0, 1], [0, 1, 0]), -np.array(center))

        arrow_z = o3d.geometry.TriangleMesh.create_arrow(width, width*2, cylinder_height=length, cone_height=length*.2, **kwargs)
        arrow_z.translate(-np.array(center))
        arrow_z.paint_uniform_color([0, 0, 1])

        self.add_mesh(arrow_x)
        self.add_mesh(arrow_y)
        self.add_mesh(arrow_z)


if __name__ == '__main__':
    from tcintegral import basis_set, MolecularOrbital

    # get a basis set from basis_sets
    bs = basis_set.STO2G

    # read a molecule to use
    mol = plams.Molecule('/Users/yumanhordijk/Desktop/benzene.xyz')
    mol.guess_bonds()

    # get atomic orbitals for this molecule
    # orb_to_get = ['1S', '2S', '3S', '1P:x', '1P:y', '1P:z', '2P:x', '2P:y', '2P:z']  # list of possible AO
    orb_to_get = ['1P:z']  # list of possible AOs
    aos = []
    for atom in mol:
        for orb_name in orb_to_get:
            try:  # try to get the AO for this atom, it will fail if it does not exist
                aos.append(bs.get(atom.symbol, orb_name, atom.coords))
            except IndexError:
                pass

    # build the Hamiltonian
    ionization_potentials = {  # in Hartrees, obtained at the OLYP/QZ4P level in ADF
        'H(1S)': -0.238455,
        'C(1S)': -10.063764,
        'C(2S)': -0.506704,
        'C(1P:x)': -0.188651,
        'C(1P:y)': -0.188651,
        'C(1P:z)': -0.188651,
        'O(1S)': -18.9254,
        'O(2S)': -0.886887,
        'O(1P:x)': -0.329157,
        'O(1P:y)': -0.329157,
        'O(1P:z)': -0.329157,
    }

    K = 1.75
    H = np.zeros((len(aos), len(aos)))
    # build diagonal elements
    for i, ao in enumerate(aos):
        H[i, i] = ionization_potentials.get(ao.name)

    # build off-diagonal elements
    for i, ao1 in enumerate(aos):
        for j, ao2 in enumerate(aos[i+1:], start=i+1):
            H[i, j] = H[j, i] = K * ao1.overlap(ao2) * (H[i, i] + H[j, j]) / 2

    # get MO energies and coeffs
    energies, coefficients = np.linalg.eigh(H)

    for mo_index in range(len(energies)):
        mo = MolecularOrbital(aos, coefficients[:, mo_index], mol)
        with Screen(screenshot_file='test.png') as scr:
            scr.draw_molecule(mol)
            scr.draw_orbital(mo)
            scr.draw_axes()