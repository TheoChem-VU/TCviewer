try:
    from PySide6 import *
    has_qt = True
except ImportError:
    has_qt = False

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


class Screen:
    def __new__(cls, headless=False):
        if headless or not has_qt:
            return _HeadlessScreen()
        else:
            return _Screen()


class _HeadlessScreen:
    def __init__(self):
        self.molview = mol_widget.MoleculeWidget(headless=True)
        self.use_parallel_projection = False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # self.exec()
        pass

    def draw_molecule(self, *args, **kwargs):
        self.molview.draw_molecule(*args, **kwargs)

    def add_molscene(self):
        return self.molview.new_scene()

    def screenshot(self, *args, **kwargs):
        self.molview.screenshot(*args, **kwargs)

    def screenshots(self, *args, **kwargs):
        self.molview.screenshots(*args, **kwargs)



if has_qt:
    class _Screen(QtWidgets.QApplication):
        def __post_init__(self):
            self.use_parallel_projection = False
            self.window = QtWidgets.QMainWindow()
            self.window.layout = QtWidgets.QGridLayout()
            grid_widget = QtWidgets.QWidget()
            grid_widget.setLayout(self.window.layout)
            self.window.setCentralWidget(grid_widget)
            self.window.setWindowTitle("TCViewer 2.0")
            
            # set up the settings screen
            self.settings = settings.DefaultSettings()
            self.window.layout.addWidget(self.settings, 0, 3, 2, 1)

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


        def __enter__(self):
            self.__post_init__()
            return self

        def __exit__(self, *args):
            self.window.show()
            self.exec()
            self.shutdown()

        def draw_molecule(self, *args, **kwargs):
            self.molview.draw_molecule(*args, **kwargs)

        def add_molscene(self):
            return self.molview.new_scene()

        def screenshot(self, *args, **kwargs):
            self.molview.screenshot(*args, **kwargs)

        def screenshots(self, *args, **kwargs):
            self.molview.screenshots(*args, **kwargs)

if __name__ == '__main__':
    with Screen(headless=False) as scr:
        ...
        # with scr.add_molscene() as scene:
        #     # res = tcutility.results.read('/Users/yumanhordijk/PhD/Projects/RadicalAdditionASMEDA/data/DFT/TS_C_O/PyFrag_OLYP_TZ2P/frag_Substrate')
        #     orbs = pyfmo.orbitals2.objects.Orbitals('/Users/yumanhordijk/Downloads/FragAnal.adf.rkf')

        #     cub1 = orbs.sfos['1(SOMO)'].cube_file()
        #     cub2 = orbs.sfos['2(21A)'].cube_file()

        #     S = cub1.copy()
        #     S.values = abs(cub1.values * cub2.values)
        #     mol = cub1.molecule

        #     # T = tcutility.geometry.MolTransform(mol)
        #     # T.center(1)  # center on C1
        #     # T.align_to_vector(1, 2, [1, 0, 0])  # put C=O bond on x-axis
        #     # T.align_to_plane(2, 1, 4, [0, 1, 0])  # place the molecule on the xz-plane
        #     # # T.rotate(x=7 * np.pi / 180)
        #     # scene.transform = T.to_vtkTransform()

        #     scene.draw_molecule(cub1.molecule)
        #     scene.draw_isosurface(S,  0.009, [0, 1, 1])
        #     # scene.draw_isosurface(cub2, -0.03, [1, .5, 0])
        #     # scene.draw_isosurface(cub2,  0.03, [0, 1, 1])

        # scr.screenshots(directory='screenshots')
