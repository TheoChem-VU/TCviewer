from PySide6 import *

import vtk
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkAssembly, vtkFollower, vtkPolyDataMapper, vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkLight, vtkCamera
from vtkmodules.vtkFiltersSources import vtkLineSource, vtkArcSource, vtkSphereSource, vtkRegularPolygonSource
from vtkmodules.vtkFiltersCore import vtkTubeFilter

from tcviewer import mol_widget
from tcviewer.settings import settings
import tcutility
import pyfmo
from scm import plams
import numpy as np
import os


class MoleculeScene:
    def __init__(self, parent):
        self.parent = parent
        self.renderer = vtkRenderer()
        self.renderer.SetBackground([1, 1, 1])
        self.renderer.UseFXAAOn()
        self.parent.renWin.AddRenderer(self.renderer)

        light = vtkLight()
        light.SetPosition(-1.5, 2, 2)
        light.SetLightTypeToCameraLight()
        self.renderer.AddLight(light)
        self.save_camera()

        self.transform = vtk.vtkTransform()
        self.camera_followers = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.post_draw()
        pass

    def screenshot(self, path: str):
        self.post_draw()
        img_filter = vtk.vtkWindowToImageFilter()
        img_filter.SetInput(self.parent.renWin)
        img_filter.SetScale(2)
        # img_filter.SetMagnification(2)
        img_filter.SetInputBufferTypeToRGBA()
        img_filter.ReadFrontBufferOff()
        img_filter.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(path)
        writer.SetInputConnection(img_filter.GetOutputPort())
        writer.Write()

    #   vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    #   windowToImageFilter->SetInput(renderWindow);
    # #if VTK_MAJOR_VERSION >= 8 || VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 90
    #   windowToImageFilter->SetScale(2); // image quality
    # #else
    #   windowToImageFilter->SetMagnification(2); // image quality
    # #endif
    #   windowToImageFilter->SetInputBufferTypeToRGBA(); // also record the alpha
    #                                                    // (transparency) channel
    #   windowToImageFilter->ReadFrontBufferOff();       // read from the back buffer
    #   windowToImageFilter->Update();

    #   vtkNew<vtkPNGWriter> writer;
    #   writer->SetFileName("screenshot2.png");
    #   writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    #   writer->Write();

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
            follower['actor'].SetPosition(self.transform.TransformPoint(pos))
        
        # for actor in self.renderer.GetActors():
        #     if actor in [follower['actor'] for follower in self.camera_followers]:
        #         continue
        #     actor.SetUserTransform(self.transform)

        self.renderer.ResetCamera()
        self.parent.renWin.Render()
        self.parent.Initialize()
        self.parent.Start()
        self.save_camera()

        self.post_draw = lambda: ...

    def draw_molecule(self, mol):
        for atom in mol:
            self.draw_atom(atom)

        mol.guess_bonds()
        for bond in mol.bonds:
            self.draw_single_bond(bond.atom1, bond.atom2)

    def draw_atom(self, atom):
        def draw_disk(rotatex, rotatey):
            circle = vtkRegularPolygonSource()
            circle.SetCenter([0, 0, 0])
            circle.SetRadius(tcutility.data.atom.radius(atom.symbol) * settings['atom']['size'] + settings['atom']['quadrant_width'])
            circle.SetNumberOfSides(50)
            circleMapper = vtkPolyDataMapper()
            circleMapper.SetInputConnection(circle.GetOutputPort())
            if settings['atom']['quadrant_follow_camera']:
                circleActor = vtkFollower()
                circleActor.SetCamera(self.renderer.GetActiveCamera())
            else:
                circleActor = vtkActor()
            circleActor.SetPosition(atom.coords)
            circleActor.SetMapper(circleMapper)
            circleActor.GetProperty().SetColor([0, 0, 0])
            circleActor.PickableOff()
            # circleActor.SetUserTransform(self.transform)

            self.camera_followers.append({'actor': circleActor, 'rotatex': rotatex, 'rotatey': rotatey})
            self.renderer.AddActor(circleActor)

        sphere = vtkSphereSource()
        sphere.SetPhiResolution(35)
        sphere.SetThetaResolution(45)
        sphere.SetRadius(tcutility.data.atom.radius(atom.symbol) * settings['atom']['size'])
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
        sphereActor.atom = atom
        sphereActor.SetUserTransform(self.transform)
        self.renderer.AddActor(sphereActor)

        if settings['atom']['draw_quadrants']:
            draw_disk(0, 0)
            draw_disk(-65, 0)
            draw_disk(0, -65)

    def draw_single_bond(self, a1, a2):
        p1, p2 = np.array(a1.coords), np.array(a2.coords)
        lineSource = vtkLineSource()
        lineSource.SetPoint1(p1)
        lineSource.SetPoint2(p2)

        tubeFilter = vtkTubeFilter()
        tubeFilter.source = lineSource
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(settings['bond']['radius'])
        tubeFilter.SetNumberOfSides(20)

        tubeMapper = vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

        tubeActor = vtkActor()
        tubeActor.SetMapper(tubeMapper)
        tubeActor.GetProperty().SetColor(settings['bond']['color'])
        tubeActor.type = 'bond'
        tubeActor.atoms = a1, a2
        tubeActor.SetUserTransform(self.transform)
        self.renderer.AddActor(tubeActor)

    def remove_bond(self, a1, a2):
        for actor in self.renderer.GetActors():
            if not hasattr(actor, 'type'):
                continue

            if actor.type != 'bond':
                continue

            if a1 in actor.atoms and a2 in actor.atoms:
                self.renderer.RemoveActor(actor)

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
        actor.type = 'surface'
        actor.SetUserTransform(self.transform)

        self.renderer.AddActor(actor)
        # self.Render()

    def draw_angle(self, a1, a2, a3):
        # c1, c2, c3 = np.array(a1.coords), np.array(a2.coords), np.array(a3.coords)
        c1 = np.array(a1.coords)
        c2 = np.array(a2.coords)
        c3 = np.array(a3.coords)

        arcSource = vtk.vtkArcSource()

        u, v = c1 - c2, c3 - c2
        u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)

        arcSource.SetCenter(c2)
        arcSource.SetPoint1(u * .7 + v * .2 + c2)
        arcSource.SetPoint2(v * .7 + u * .2 + c2)
        arcSource.SetResolution(100)

        tubeFilter = vtkTubeFilter()
        tubeFilter.source = arcSource
        tubeFilter.SetInputConnection(arcSource.GetOutputPort())
        tubeFilter.SetRadius(settings['bond']['radius'] / 5)
        tubeFilter.SetNumberOfSides(20)

        tubeMapper = vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

        tubeActor = vtkActor()
        tubeActor.SetMapper(tubeMapper)
        tubeActor.GetProperty().SetColor([0, 0, 0])
        tubeActor.SetUserTransform(self.transform)
        self.renderer.AddActor(tubeActor)

        # angle = tcutility.geometry.parameter([c1, c2, c3], 0, 1, 2)

        # text = vtk.vtkVectorText()
        # text.SetText(f'{angle: .1f}°')
        # textMapper = vtk.vtkPolyDataMapper()
        # textMapper.SetInputConnection(text.GetOutputPort())
        # textActor = vtk.vtkFollower()
        # textActor.SetMapper(textMapper)
        # # textActor.SetUserTransform(self.transform)
        # textActor.SetCamera(self.renderer.GetActiveCamera())
        # textActor.SetScale(0.2, 0.2, 0.2)
        # textActor.SetPosition(self.transform.TransformPoint(c2 + (u + v) / 1.5))
        # # textActor.SetActiveCamera(self._base_ren.GetActiveCamera())
        # self.camera_followers.append({'actor': textActor, 'rotatex': 0, 'rotatey': 0})
        # self.renderer.AddActor(textActor)


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
        self._base_ren.DrawOff()
        self.renWin.AddRenderer(self._base_ren)
        self.SetRenderWindow(self.renWin)
        self.installEventFilter(MoleculeWidgetKeyPressFilter(parent=self))
        self.setAcceptDrops(True)
        self.Initialize()
        self.Start()

        self.selected_actors = []
        self.selected_actor_highlights = {}
        self.picker = vtk.vtkPropPicker()
        self.AddObserver('EndInteractionEvent', self.highlight_observer)
        # self.GetInteractorStyle().AddObserver(vtk.vtkCommand.LeftButtonReleaseEvent, self.highlight_observer)
        self.AddObserver('LeftButtonPressEvent', self.record_mouse_position)

        self._recording_mouse = False
        self._mouse_pos = None

    def record_mouse_position(self, interactor, event):
        self._recording_mouse = True
        self._mouse_pos = interactor.GetEventPosition()

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
                self.draw_molecule(fname)
        else:
            e.ignore()

    # The following three methods set up dragging and dropping for the app
    def highlight_observer(self, interactor, event):
        # check if an actor was picked by the user
        pick = self.picker.PickProp(*self.GetEventPosition(), self.active_scene.renderer)

        # if we did not pick an actor we check if we clicked (not dragged) outside
        if not pick:
            if self._mouse_pos == self.GetEventPosition():
                self.remove_all_highlights()
            return

        actor = self.picker.GetViewProp()

        # lets toggle the selected actor
        self.toggle_highlight(actor)
        # if the user did not hold the shift key then we deselect the other actors
        if not self.GetShiftKey():
            self.remove_all_highlights(exception=actor)

        if self._recording_mouse:
            self.LeftButtonReleaseEvent()
            self._recording_mouse = False


    def add_highlight(self, actor):
        ren = self.active_scene.renderer
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
            tube = actor.GetMapper().GetInputConnection(0, 0).GetProducer()
            source = tube.source
            line = vtkLineSource()
            line.SetPoint1(*source.GetPoint1())
            line.SetPoint2(*source.GetPoint2())
            highlight = vtkTubeFilter()
            highlight.SetInputConnection(line.GetOutputPort())
            highlight.SetRadius(tube.GetRadius())
            highlight.SetNumberOfSides(20)
            tube.SetRadius(tube.GetRadius() / 1.5)

        highlightMapper = vtkPolyDataMapper()
        highlightMapper.SetInputConnection(highlight.GetOutputPort())
        highlightActor = vtkActor()
        highlightActor.SetPosition(actor.GetPosition())
        highlightActor.SetMapper(highlightMapper)
        highlightActor.GetProperty().SetColor([0, 1, 1])
        highlightActor.GetProperty().SetOpacity(.5)
        highlightActor.PickableOff()
        highlightActor.SetUserTransform(self.active_scene.transform)
        ren.AddActor(highlightActor)
        self.selected_actor_highlights[actor] = highlightActor

    def remove_highlight(self, actor):
        ren = self.active_scene.renderer
        self.selected_actors = [actor_ for actor_ in self.selected_actors if actor_ != actor]
        highlightActor = self.selected_actor_highlights.pop(actor)
        actor.SetScale([1, 1, 1])

        if actor.type == 'bond':
            tube = actor.GetMapper().GetInputConnection(0, 0).GetProducer()
            tube.SetRadius(tube.GetRadius() * 1.5)

        ren.RemoveActor(highlightActor)

    def remove_all_highlights(self, exception=None):
        exception = tcutility.ensure_list(exception)
        for actor in self.selected_actors:
            if actor in exception:
                continue
            self.remove_highlight(actor)

    def toggle_highlight(self, actor):
        if actor in self.selected_actors:
            self.remove_highlight(actor)
        else:
            self.add_highlight(actor)

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

        with self.new_scene() as scene:
            scene.renderer.SetActiveCamera(self._base_ren.GetActiveCamera())
            scene.draw_molecule(mol)
            if orbs:
                scene.draw_isosurface(tcutility.ensure_list(orbs.mos['LUMO'])[0].cube_file(), 0.03, [0, 1, 1])
                scene.draw_isosurface(tcutility.ensure_list(orbs.mos['LUMO'])[0].cube_file(), -0.03, [1, 1, 0])
        self.set_active_mol(-1)

    def new_scene(self):
        scene = MoleculeScene(self)
        scene.renderer.SetActiveCamera(self._base_ren.GetActiveCamera())
        [scene.renderer.DrawOff() for scene in self.scenes]
        self.scenes.append(scene)
        self.set_active_mol(-1)
        return scene

    def next_mol(self):
        self.set_active_mol(self.active_scene_index + 1)

    def previous_mol(self):
        self.set_active_mol(self.active_scene_index - 1)

    def set_active_mol(self, index):
        if len(self.scenes) == 0:
            return

        self.remove_all_highlights()

        self.active_scene.save_camera()

        index = index % len(self.scenes)
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

    @property
    def active_scene_index(self):
        if len(self.scenes) == 0:
            return
        return [i for i, scene in enumerate(self.scenes) if scene.renderer.GetDraw()][0]

    @property
    def active_scene(self):
        return self.scenes[self.active_scene_index]

    @property
    def number_of_scenes(self):
        return len(self.scenes)

    def screenshot(self, path):
        self.active_scene.screenshot(path)

    def screenshots(self, paths=None, directory=None):
        if paths is None and directory is None:
            raise ValueError('You should give either the paths or directory argument')

        if paths is None:
            os.makedirs(directory, exist_ok=True)
            paths = [os.path.join(directory, f'scene{i}.png') for i in range(self.number_of_scenes)]

        for scene, path in zip(self.scenes, paths):
            scene.screenshot(path)



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

            scene = self.parent().scenes[self.parent().active_scene_index]
            bond_selected = len(self.parent().selected_actors) == 1 and self.parent().selected_actors[0].type == 'bond'
            atoms2_selected = len(self.parent().selected_actors) == 2 and self.parent().selected_actors[0].type == 'atom' and self.parent().selected_actors[1].type == 'atom'
            if atoms2_selected:
                actor1, actor2 = self.parent().selected_actors
                if event.key() == QtCore.Qt.Key_1:
                    scene.draw_single_bond(actor1.atom, actor2.atom)
                    self.parent().Render()

                if event.key() in [QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete]:
                    scene.remove_bond(actor1.atom, actor2.atom)
                    self.parent().Render()

            if bond_selected:
                actor = self.parent().selected_actors[0]
                if event.key() in [QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete]:
                    scene.remove_bond(*actor.atoms)
                    self.parent().remove_highlight(actor)
                    self.parent().Render()

        return False
