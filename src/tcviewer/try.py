# coding=utf-8


from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer
# load implementations for rendering and interaction factory classes
import vtkmodules.vtkRenderingOpenGL2
from vtk import vtkLineSource, vtkTubeFilter, vtkInteractorStyleTrackballCamera, vtkLight
import numpy as np

import vtkmodules.qt.QVTKRenderWindowInteractor as QVTK
QVTKRenderWindowInteractor = QVTK.QVTKRenderWindowInteractor

if QVTK.PyQtImpl == 'PySide6':
    from PySide6.QtCore import Qt
    from PySide6 import QtWidgets
    from PySide6.QtWidgets import QApplication, QMainWindow
elif QVTK.PyQtImpl == 'PySide2':
    from PySide2.QtCore import Qt
    from PySide2 import QtWidgets
    from PySide2.QtWidgets import QApplication, QMainWindow
else:
    from PySide.QtCore import Qt
    from PySide import QtWidgets
    from PySide.QtGui import QApplication, QMainWindow

from tcviewer import shapes, widgets


class MolWidget(QVTKRenderWindowInteractor):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        self.renderer = vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)
        self.GetRenderWindow().AddRenderer(self.renderer)

        light = vtkLight()
        light.SetPosition(-1.5, 2, 2)
        light.SetLightTypeToCameraLight()
        self.renderer.AddLight(light)

        self.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.Initialize()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.Start()

    def draw_single_bond(self, p1, p2, color=None, opacity=1, radius=.05):
        p1, p2 = np.array(p1), np.array(p2)
        lineSource = vtkLineSource()
        lineSource.SetPoint1(p1)
        lineSource.SetPoint2(p2)

        tubeFilter = vtkTubeFilter()
        tubeFilter.source = lineSource
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(radius)

        tubeFilter.SetCapping(True)
        tubeFilter.SetNumberOfSides(20)

        tubeMapper = vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

        tubeActor = vtkActor()
        tubeActor.SetMapper(tubeMapper)
        tubeActor.GetProperty().SetColor([0, 0, 0])

        tubeActor.GetProperty().SetOpacity(opacity)

        tubeActor.type = 'bond'
        # a1, a2 = [atom for atom in self.mol if atom.coords == p1][0], [atom for atom in self.mol if atom.coords == p2][0]
        tubeActor.atoms = p1.tolist(), p2.tolist()
        # tubeActor.SetUserTransform(self.transform)
        self.renderer.AddActor(tubeActor)
        return tubeActor

    def draw_atom(self, atom, color=None, opacity=1, ambient=0.65, draw_outline=True, draw_quadrants=True):
        shapes.Atom(self.renderer, atom.symbol, atom.coords)


class TCViewerWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.layout = QtWidgets.QHBoxLayout(self)
        self.molwidget_stack = QtWidgets.QStackedWidget(self)
        # self.layout.addWidget(self.molwidget_stack)
        widgets.PeriodicTable(self)

    def add_molscene(self):
        scene = MolWidget(self)
        self.molwidget_stack.addWidget(scene)
        return scene


class TCViewer(QApplication):
    def __init__(self):
        super().__init__()
        self.window = TCViewerWindow()
        self.window.show()

    def __enter__(self):
        return self.window

    def __exit__(self, *args):
        try:
            self.exec()
        except AttributeError:
            self.exec_()


if __name__ == "__main__":
    from scm import plams

    mol = plams.Molecule(r'D:\Users\Yuman\Desktop\PhD\TCMU\examples\job\asc.xyz')

    with TCViewer() as scr:
        with scr.add_molscene() as scene:
            for atom in mol:
                scene.draw_atom(atom)

            mol.guess_bonds()
            for bond in mol.bonds:
                scene.draw_single_bond(bond.atom1.coords, bond.atom2.coords)
