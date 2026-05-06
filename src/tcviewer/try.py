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
    from PySide6 import QtWidgets, QtGui, QtCore
    from PySide6.QtWidgets import QApplication, QMainWindow
elif QVTK.PyQtImpl == 'PySide2':
    from PySide2.QtCore import Qt
    from PySide2 import QtWidgets, QtGui, QtCore
    from PySide2.QtWidgets import QApplication, QMainWindow
else:
    from PySide.QtCore import Qt
    from PySide import QtWidgets, QtGui, QtCore
    from PySide.QtGui import QApplication, QMainWindow

from tcviewer import shapes, widgets, settings


class MolWidget(QtWidgets.QSplitter):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.scene = MolScene(self)
        self.settings = settings.SettingsWidget()

        self.addWidget(self.scene)
        self.addWidget(self.settings)

    def __enter__(self):
        return self.scene

    def __exit__(self, *args):
        self.scene.renderer.ResetCamera()
        self.scene.Start()


class MolScene(QVTKRenderWindowInteractor):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        self.shapes = []

        self.renderer = vtkRenderer()
        self.renderer.UseFXAAOn()
        self.renderer.SetBackground(1, 1, 1)
        self.GetRenderWindow().AddRenderer(self.renderer)

        light = vtkLight()
        light.SetPosition(-1.5, 2, 2)
        light.SetLightTypeToCameraLight()
        self.renderer.AddLight(light)

        self.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.Initialize()

        self.setFocusPolicy(Qt.NoFocus)

    def draw_single_bond(self, a1, a2, **kwargs):
        idx1 = a1.mol.atoms.index(a1) + 1
        idx2 = a2.mol.atoms.index(a2) + 1

        bond = shapes.Bond(self.renderer, a1, a2, name=f'Bond({a1.symbol}:{idx1} ––– {a2.symbol}:{idx2})', **kwargs)

        self.shapes.append(bond)
        self.parent.settings.add_object(bond)
        return bond

    def draw_atom(self, atom, **kwargs):
        idx = atom.mol.atoms.index(atom) + 1

        atom = shapes.Atom(self.renderer, atom.symbol, atom.coords, name=f'Atom({atom.symbol}:{idx})', **kwargs)
        self.shapes.append(atom)
        self.parent.settings.add_object(atom)
        return atom


class TCViewerWindowKeyFilter(QtCore.QObject):
    def eventFilter(self, widget, event):
        # event 51 is the keyevent
        if event.type() != 51:
            return

        if event.key() == QtCore.Qt.Key_Right:
            self.parent().next_scene()

        if event.key() == QtCore.Qt.Key_Left:
            self.parent().previous_scene()

        # spacebar resets the camera to the right distance
        if event.key() == QtCore.Qt.Key_Space:
            self.parent().current_scene().renderer.ResetCamera()
            self.parent().current_scene().GetRenderWindow().Render()


class TCViewerWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.main_frame = QtWidgets.QFrame(parent=self)
        self.setCentralWidget(self.main_frame)

        self.layout = QtWidgets.QGridLayout(self.main_frame)
        self.molwidget_stack = QtWidgets.QStackedWidget(parent=self)
        self.molwidget_stack.setFocusPolicy(Qt.NoFocus)

        prev_button = QtWidgets.QPushButton('<', parent=self)
        prev_button.clicked.connect(self.previous_scene)

        next_button = QtWidgets.QPushButton('>', parent=self)
        next_button.clicked.connect(self.next_scene)

        self.installEventFilter(TCViewerWindowKeyFilter(parent=self))

        self.scene_slider = QtWidgets.QScrollBar(parent=self)
        self.scene_slider.setOrientation(Qt.Horizontal)
        self.scene_slider.setMinimum(0)
        self.scene_slider.setMaximum(0)
        self.scene_slider.valueChanged.connect(self.molwidget_stack.setCurrentIndex)

        self.layout.addWidget(self.molwidget_stack, 0, 0, 1, 3)
        self.layout.addWidget(prev_button, 1, 0, 1, 1)
        self.layout.addWidget(self.scene_slider, 1, 1, 1, 1)
        self.layout.addWidget(next_button, 1, 2, 1, 1)


    def add_molscene(self):
        scene = MolWidget(self)
        self.molwidget_stack.addWidget(scene)
        self.scene_slider.setMaximum(self.molwidget_stack.count() - 1)
        return scene

    def previous_scene(self):
        idx = (self.molwidget_stack.currentIndex() - 1) % self.molwidget_stack.count()
        self.scene_slider.setSliderPosition(idx)
        self.molwidget_stack.setCurrentIndex(idx)

    def next_scene(self):
        idx = (self.molwidget_stack.currentIndex() + 1) % self.molwidget_stack.count()
        self.scene_slider.setSliderPosition(idx)
        self.molwidget_stack.setCurrentIndex(idx)

    def current_scene(self):
        return self.molwidget_stack.currentWidget().scene


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

    mol = plams.Molecule('/Users/yumanhordijk/PhD/Programs/TheoCheM/TCMU/examples/job/asc.xyz')

    with TCViewer() as scr:
        with scr.add_molscene() as scene:
            for atom in mol:
                a = scene.draw_atom(atom)

            mol.guess_bonds()
            for bond in mol.bonds:
                scene.draw_single_bond(bond.atom1, bond.atom2)

        with scr.add_molscene() as scene:
            for atom in mol:
                a = scene.draw_atom(atom)
                a.set_quadrant_width(0.1)

            mol.guess_bonds()
            for bond in mol.bonds:
                scene.draw_single_bond(bond.atom1, bond.atom2)

        with scr.add_molscene() as scene:
            for atom in mol:
                a = scene.draw_atom(atom)
                a.set_quadrant_width(0.0)

            mol.guess_bonds()
            for bond in mol.bonds:
                scene.draw_single_bond(bond.atom1, bond.atom2)
