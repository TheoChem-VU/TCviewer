from PySide6 import *

import vtk
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkAssembly, vtkFollower, vtkPolyDataMapper, vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkLight, vtkCamera
from vtkmodules.vtkFiltersSources import vtkLineSource, vtkSphereSource, vtkRegularPolygonSource
from vtkmodules.vtkFiltersCore import vtkTubeFilter

from tcviewer import mol_widget, settings
import tcutility
import pyfmo
from scm import plams
import numpy as np


class Screen(QtWidgets.QApplication):
    def __post_init__(self):
        self.window = QtWidgets.QMainWindow()
        self.window.layout = QtWidgets.QGridLayout()
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(self.window.layout)
        self.window.setCentralWidget(grid_widget)
        self.window.setWindowTitle("TCViewer 2.0")

        # set up the molecule screen
        self.molview = mol_widget.MoleculeWidget(self)
        self.window.layout.addWidget(self.molview, 0, 0, 1, 3)

        # make a scrollbar for changing the molecule
        self.molviewslider = QtWidgets.QScrollBar()
        self.molviewslider.setOrientation(QtCore.Qt.Horizontal)
        self.molviewslider.setMinimum(0)
        self.molviewslider.setMaximum(0)
        self.molviewslider.valueChanged.connect(self.molview.set_active_mol)
        self.window.layout.addWidget(self.molviewslider, 1, 1)

        self.molviewslider_prevbtn = QtWidgets.QPushButton('<')
        self.molviewslider_prevbtn.clicked.connect(self.molview.previous_mol)
        self.molviewslider_nextbtn = QtWidgets.QPushButton('>')
        self.molviewslider_nextbtn.clicked.connect(self.molview.next_mol)
        self.window.layout.addWidget(self.molviewslider_prevbtn, 1, 0)
        self.window.layout.addWidget(self.molviewslider_nextbtn, 1, 2)

        self.window.layout.setColumnStretch(0, 0)
        self.window.layout.setColumnStretch(1, 1)
        self.window.layout.setColumnStretch(2, 0)

        self.settings = settings.DefaultSettings()
        self.window.layout.addWidget(self.settings, 0, 3, 2, 1)


    def __enter__(self):
        self.__post_init__()
        return self

    def __exit__(self, *args):
        self.window.show()
        self.exec()

    def draw_molecule(self, *args, **kwargs):
        self.molview.draw_molecule(*args, **kwargs)

    def add_molscene(self):
        return self.molview.new_scene()

    def screenshot(self, *args, **kwargs):
        self.molview.screenshot(*args, **kwargs)

    def screenshots(self, *args, **kwargs):
        self.molview.screenshots(*args, **kwargs)

        
if __name__ == '__main__':
    with Screen() as scr:
        with scr.add_molscene() as scene:
            res = tcutility.results.read('/Users/yumanhordijk/PhD/Projects/RadicalAdditionASMEDA/data/DFT/TS_C_O/PyFrag_OLYP_TZ2P/frag_Substrate')
            orbs = pyfmo.orbitals.Orbitals(res.files['adf.rkf'])

            cub = orbs.mos['LUMO'].cube_file()
            mol = cub.molecule
            v_cx = mol.as_array()[1] - mol.as_array()[0]

            T = scene.transform
            T.PostMultiply()
            T.Translate(*(-np.mean(mol.as_array(), axis=0)).tolist())

            R_x = tcutility.geometry.vector_align_rotmat(v_cx, [1, 0, 0])
            angles_x = tcutility.geometry.rotmat_to_angles(R_x)
            T.RotateX(angles_x[0] * 180 / np.pi)
            T.RotateY(angles_x[1] * 180 / np.pi)
            T.RotateZ(angles_x[2] * 180 / np.pi)

            mol_ = tcutility.geometry.apply_rotmat(mol.as_array(), R_x)
            n_cxh = np.cross(mol_[1] - mol_[0], mol_[-1] - mol_[0])
            R_y = tcutility.geometry.vector_align_rotmat(n_cxh, [0, 1, 0])
            angles_y = tcutility.geometry.rotmat_to_angles(R_y)
            T.RotateX(angles_y[0] * 180 / np.pi)
            T.RotateY(angles_y[1] * 180 / np.pi)
            T.RotateZ(angles_y[2] * 180 / np.pi)

            T.RotateX(7)
            # scene.scene_assembly.SetUserTransform(T)

            actor = scene.draw_molecule(mol)
            actor = scene.draw_isosurface(cub, -0.03, [1, 1, 0])
            actor = scene.draw_isosurface(cub,  0.03, [0, 1, 1])

        scr.screenshots(directory='screenshots')