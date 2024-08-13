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

from tcviewer import mol_widget
import tcutility
import pyfmo
from scm import plams
import numpy as np


class MoleculeScene:
    def __init__(self, parent):
        self.parent = parent
        self.renderer = vtkRenderer()
        self.renderer.SetBackground([1, 1, 1])
        self.parent.renWin.AddRenderer(self.renderer)

        light = vtkLight()
        light.SetPosition(-1.5, 2, 2)
        light.SetLightTypeToCameraLight()
        self.renderer.AddLight(light)
        self.save_camera()

        self.scene_assembly = vtkAssembly()
        self.transform = vtk.vtkTransform()
        self.camera_followers = []
        # self.renderer.AddActor(self.scene_assembly)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.post_draw()
        pass

    def save_camera(self):
        camera = self.renderer.GetActiveCamera()
        self.cam_settings = {}
        self.cam_settings['position'] = camera.GetPosition()
        self.cam_settings['focal point'] = camera.GetFocalPoint()
        self.cam_settings['view up'] = camera.GetViewUp()
        self.cam_settings['distance'] = camera.GetDistance()
        self.cam_settings['clipping range'] = camera.GetClippingRange()
        self.cam_settings['orientation'] = camera.GetOrientation()


    def load_camera(self):
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(self.cam_settings['position'])
        camera.SetFocalPoint(self.cam_settings['focal point'])
        camera.SetViewUp(self.cam_settings['view up'])
        camera.SetDistance(self.cam_settings['distance'])
        camera.SetClippingRange(self.cam_settings['clipping range'])

    def post_draw(self):
        for follower in self.camera_followers:
            follower['actor'].RotateX(follower['rotatex'])
            follower['actor'].RotateY(follower['rotatey'])
            pos = follower['actor'].GetPosition()
            T = self.scene_assembly.GetUserTransform()
            follower['actor'].SetPosition(T.TransformPoint(pos))

        self.renderer.AddActor(self.scene_assembly)
        self.renderer.ResetCamera()
        self.parent.renWin.Render()
        self.parent.Initialize()
        self.parent.Start()

    def draw_molecule(self, mol):
        for atom in mol:
            self.draw_atom(atom)

        mol.guess_bonds()
        for bond in mol.bonds:
            self.draw_single_bond(bond.atom1.coords, bond.atom2.coords)

    def draw_atom(self, atom):
        # actors = []
        def draw_disk(rotatex, rotatey):
            circle = vtkRegularPolygonSource()
            circle.SetCenter([0, 0, 0])
            circle.SetRadius(tcutility.data.atom.radius(atom.symbol)/2.4 + .02)
            circle.SetNumberOfSides(50)
            circleMapper = vtkPolyDataMapper()
            circleMapper.SetInputConnection(circle.GetOutputPort())
            circleActor = vtkFollower()
            circleActor.SetCamera(self.renderer.GetActiveCamera())
            circleActor.SetPosition(atom.coords)
            circleActor.SetMapper(circleMapper)
            circleActor.GetProperty().SetColor([0, 0, 0])
            circleActor.PickableOff()

            self.camera_followers.append({'actor': circleActor, 'rotatex': rotatex, 'rotatey': rotatey})
            self.renderer.AddActor(circleActor)

        sphere = vtkSphereSource()
        sphere.SetPhiResolution(35)
        sphere.SetThetaResolution(45)
        sphere.SetRadius(tcutility.data.atom.radius(atom.symbol)/2.4)
        sphereMapper = vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())
        sphereActor = vtkActor()
        sphereActor.SetMapper(sphereMapper)
        sphereActor.SetPosition(atom.coords)
        sphereActor.GetProperty().SetAmbient(0.65)
        sphereActor.GetProperty().SetDiffuse(0.5)
        sphereActor.GetProperty().SetSpecular(0.5)
        sphereActor.GetProperty().SetSpecularPower(5.0)
        sphereActor.GetProperty().SetColor([x/255 for x in tcutility.data.atom.color(atom.symbol)])
        sphereActor.type = 'atom'
        self.scene_assembly.AddPart(sphereActor)

        draw_disk(0, 0)
        draw_disk(-65, 0)
        draw_disk(0, -65)

        # return actors

    def draw_single_bond(self, p1, p2):
        p1, p2 = np.array(p1), np.array(p2)
        lineSource = vtkLineSource()
        lineSource.SetPoint1(p1)
        lineSource.SetPoint2(p2)

        tubeFilter = vtkTubeFilter()
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(0.075)
        tubeFilter.SetNumberOfSides(20)

        tubeMapper = vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

        tubeActor = vtkActor()
        tubeActor.SetMapper(tubeMapper)
        tubeActor.GetProperty().SetColor([0, 0, 0])
        tubeActor.type = 'bond'

        self.scene_assembly.AddPart(tubeActor)
        # self.renderer.AddActor(tubeActor)

        # return tubeActor

    def draw_isosurface(self, grid, isovalue=0, color=(1, 1, 0)):
        # vtkImageData is the vtk image volume type
        # this is where the conversion happens
        depthArray = numpy_to_vtk(grid.values.reshape(*grid.shape).ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        imdata = vtk.vtkImageData()
        imdata.SetDimensions(grid.shape)
        imdata.SetSpacing(grid.spacing)
        imdata.SetOrigin(grid.origin)
        imdata.GetPointData().SetScalars(depthArray)
        mcplus = vtk.vtkMarchingCubes()
        mcplus.SetInputData(imdata)
        mcplus.ComputeNormalsOn()
        mcplus.ComputeGradientsOn()
        mcplus.SetValue(0, isovalue)
        mcplus.Update()
        
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(mcplus.GetOutputPort())
        mapper.ScalarVisibilityOff()

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetOpacity(.5)
        actor.GetProperty().SetAmbient(-0)
        actor.GetProperty().SetDiffuse(1)
        actor.GetProperty().SetSpecular(5)
        actor.GetProperty().SetSpecularPower(70)
        actor.GetProperty().SetSpecularColor((1, 1, 1))
        actor.PickableOff()

        # self.renderer.AddActor(actor)
        # self.post_draw()
        self.scene_assembly.AddPart(actor)

        # return actor


class MoleculeWidget(QVTKRenderWindowInteractor):
    def __init__(self, parent):
        super(MoleculeWidget, self).__init__()
        self.scenes = []

        self.parent = parent

        self.molecule_renderers = []
        self.renWin = self.GetRenderWindow()
        self.renWin.BordersOn()
        self.interactor_style = vtk.vtkInteractorStyleTrackballCamera()
        self.SetInteractorStyle(self.interactor_style)
        self._base_ren = vtkRenderer()
        self._base_ren.SetBackground(1, 1, 1)
        self.renWin.AddRenderer(self._base_ren)
        self.SetRenderWindow(self.renWin)
        self.eventFilter = MoleculeWidgetKeyPressFilter(parent=self)
        self.installEventFilter(self.eventFilter)
        self.setAcceptDrops(True)
        self.Initialize()
        self.Start()

        self.selected_actors = []
        self.selected_actor_highlights = {}
        self.picker = vtk.vtkPropPicker()
        self.AddObserver('LeftButtonPressEvent', self.left_mousebutton_press_event)
        # self.AddObserver('LeftButtonPressEvent', self.RenderCallback, -1)
        self.CreateRepeatingTimer(500)

    # def RenderCallback(self, event, interactor):
    #     self.Render()
    #     self.Render()

    # The following three methods set up dragging and dropping for the app
    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dragMoveEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        """
        Drop files directly onto the widget
        File locations are stored in fname
        :param e:
        :return:
        """
        if e.mimeData().hasUrls:
            e.setDropAction(QtCore.Qt.CopyAction)
            e.accept()
            for url in e.mimeData().urls():
                fname = str(url.toLocalFile())
                # try:
                self.draw_molecule(fname)
                # except:
                #     pass
        else:
            e.ignore()

    # The following three methods set up dragging and dropping for the app
    def left_mousebutton_press_event(self, interactor, event):
        pick = self.picker.PickProp(*interactor.GetEventPosition(), list(self.renWin.GetRenderers())[0])
        if pick:
            actor = self.picker.GetViewProp()
            if actor in self.selected_actors:
                self.remove_highlight(actor)
            else:
                self.add_highlight(actor)

    def add_highlight(self, actor):
        ren = list(self.renWin.GetRenderers())[0]
        self.selected_actors.append(actor)
        if actor.type == 'atom':
            actor.SetScale([0.75, 0.75, 0.75])
            highlight = vtkSphereSource()
            highlight.SetPhiResolution(35)
            highlight.SetThetaResolution(45)
            highlight.SetCenter([0, 0, 0])
            radius = actor.GetMapper().GetInputConnection(0, 0).GetProducer().GetRadius()
            highlight.SetRadius(radius)

        elif actor.type == 'bond':
            source = actor.GetMapper().GetInputConnection(0, 0).GetInputConnection(0, 0).GetProducer()
            line = vtkLineSource()
            line.SetPoint1(*source.GetPoint1())
            line.SetPoint2(*source.GetPoint2())
            highlight = vtkTubeFilter()
            highlight.SetInputConnection(line.GetOutputPort())
            highlight.SetRadius(source.GetRadius())
            highlight.SetNumberOfSides(20)
            source.SetRadius(source.GetRadius() / 1.5)

        highlightMapper = vtkPolyDataMapper()
        highlightMapper.SetInputConnection(highlight.GetOutputPort())
        highlightActor = vtkFollower()
        highlightActor.SetCamera(ren.GetActiveCamera())
        highlightActor.SetPosition(actor.GetPosition())
        highlightActor.SetMapper(highlightMapper)
        highlightActor.GetProperty().SetColor([0, 1, 1])
        highlightActor.GetProperty().SetOpacity(.5)
        highlightActor.PickableOff()
        ren.AddActor(highlightActor)
        self.selected_actor_highlights[actor] = highlightActor

    def remove_highlight(self, actor):
        ren = list(self.renWin.GetRenderers())[0]
        self.selected_actors.remove(actor)
        highlightActor = self.selected_actor_highlights.pop(actor)
        actor.SetScale([1, 1, 1])
        ren.RemoveActor(highlightActor)

    def remove_all_highlights(self):
        for actor in self.selected_actors:
            self.remove_highlight(actor)

    def draw_molecule(self, xyz):
        if xyz.endswith('.xyz'):
            mol = plams.Molecule(xyz)
            orbs = False
        else:
            res = tcutility.results.read(xyz)
            orbs = pyfmo.orbitals.Orbitals(res.files['adf.rkf'])
            mol = res.molecule.output

        # disable previous scenes
        [scene.renderer.DrawOff() for scene in self.scenes]

        scene = MoleculeScene(self)
        scene.renderer.SetActiveCamera(self._base_ren.GetActiveCamera())
        scene.draw_molecule(mol)
        if orbs:
            scene.draw_isosurface(tcutility.ensure_list(orbs.mos['LUMO'])[0].cube_file(), 0.03, [0, 1, 1])
            scene.draw_isosurface(tcutility.ensure_list(orbs.mos['LUMO'])[0].cube_file(), -0.03, [1, 1, 0])

        self.scenes.append(scene)
        self.set_active_mol(len(self.scenes) - 1)

    def new_scene(self):
        scene = MoleculeScene(self)
        scene.renderer.SetActiveCamera(self._base_ren.GetActiveCamera())
        [scene.renderer.DrawOff() for scene in self.scenes]
        self.scenes.append(scene)
        self.set_active_mol(len(self.scenes) - 1)
        with scene as scene:
            return scene

    def next_mol(self):
        if len(self.scenes) == 0:
            return 
            
        new_idx = (self.active_mol_index + 1) % len(self.scenes)
        self.set_active_mol(new_idx)

    def previous_mol(self):
        if len(self.scenes) == 0:
            return 

        new_idx = (self.active_mol_index - 1) % len(self.scenes)
        self.set_active_mol(new_idx)


    def set_active_mol(self, index):
        if len(self.scenes) == 0:
            return

        self.scenes[self.active_mol_index].save_camera()

        [scene.renderer.DrawOff() for scene in self.scenes]
        self.scenes[index].renderer.DrawOn()
        self.scenes[index].load_camera()
        self.scenes[index].renderer.Render()
        self.renWin.Initialize()
        self.renWin.Render()
        self.Initialize()
        self.Render()

        self.parent.molviewslider.setMaximum(len(self.scenes) - 1)
        self.parent.molviewslider.setValue(index)
        # self.interactor_style.AutoAdjustCameraClippingRangeOff()


    @property
    def active_mol_index(self):
        if len(self.scenes) == 0:
            return 
        return [i for i, scene in enumerate(self.scenes) if scene.renderer.GetDraw()][0]


class MoleculeWidgetKeyPressFilter(QtCore.QObject):
    def eventFilter(self, widget, event):
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Right:
                self.parent().next_mol()
            if event.key() == QtCore.Qt.Key_Left:
                self.parent().previous_mol()
            if event.key() == QtCore.Qt.Key_C and event.modifiers() == QtCore.Qt.ControlModifier:
                print('copy')
            if event.key() == QtCore.Qt.Key_V and event.modifiers() == QtCore.Qt.ControlModifier:
                print('paste')
        return False