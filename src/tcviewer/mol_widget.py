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
from vtkmodules.vtkRenderingCore import vtkTextActor, vtkActor, vtkAssembly, vtkFollower, vtkPolyDataMapper, vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkLight, vtkCamera
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
        self.use_parallel_projection = False

        light = vtkLight()
        light.SetPosition(-1.5, 2, 2)
        light.SetLightTypeToCameraLight()
        self.renderer.AddLight(light)
        self.save_camera()

        self.transform = vtk.vtkTransform()
        self.transform.PostMultiply()
        self.camera_followers = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.reset_camera()
        self.post_draw()
        pass

    def screenshot(self, path: str, scale=4):
        self.post_draw()
        img_filter = vtk.vtkWindowToImageFilter()
        img_filter.SetInput(self.parent.renWin)
        img_filter.SetScale(scale)
        img_filter.SetInputBufferTypeToRGBA()
        img_filter.ReadFrontBufferOff()
        img_filter.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(path)
        writer.SetInputConnection(img_filter.GetOutputPort())
        writer.Write()

    def save_camera(self):
        camera = self.renderer.GetActiveCamera()
        self.cam_settings = {}
        self.cam_settings['position'] = camera.GetPosition()
        self.cam_settings['focal point'] = camera.GetFocalPoint()
        self.cam_settings['view up'] = camera.GetViewUp()
        self.cam_settings['distance'] = camera.GetDistance()
        # self.cam_settings['clipping range'] = camera.GetClippingRange()
        self.cam_settings['parallel_projection'] = self.use_parallel_projection
        self.cam_settings['clipping range'] = (0.1, 1000)
        self.cam_settings['orientation'] = camera.GetOrientation()

    def load_camera(self):
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.SetParallelProjection(self.cam_settings['parallel_projection'])
        camera.SetPosition(self.cam_settings['position'])
        camera.SetFocalPoint(self.cam_settings['focal point'])
        camera.SetViewUp(self.cam_settings['view up'])
        camera.SetDistance(self.cam_settings['distance'])
        camera.SetClippingRange(self.cam_settings['clipping range'])

    def reset_camera(self):
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        p = camera.GetPosition()
        plen = np.linalg.norm(p)
        camera.SetPosition([0, 0, plen])
        camera.SetFocalPoint(self.cam_settings['focal point'])
        camera.SetViewUp([0, 1, 0])
        camera.SetClippingRange(self.cam_settings['clipping range'])
        self.save_camera()

    def post_draw(self):
        for follower in self.camera_followers:
            if not follower.get('rotated', False):
                follower['actor'].RotateX(follower['rotatex'])
                follower['actor'].RotateY(follower['rotatey'])
                follower['rotated'] = True

            pos = follower['orig_pos']
            follower['actor'].SetPosition(self.transform.TransformPoint(pos))
        
        self.renderer.ResetCamera()
        self.parent.renWin.Render()
        self.parent.Initialize()
        self.save_camera()

    def draw_molecule(self, mol):
        self.mol = mol

        for atom in mol:
            self.draw_atom(atom)

        mol.guess_bonds()
        for bond in mol.bonds:
            self.draw_single_bond(bond.atom1.coords, bond.atom2.coords)

    def draw_atom(self, atom, color=None):
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
            circleActor.type = f'quadrant_{atom.symbol}'

            self.camera_followers.append({'actor': circleActor, 'rotatex': rotatex, 'rotatey': rotatey, 'cumrotatex': 0, 'cumrotatey': 0, 'orig_pos': atom.coords})
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
        if color is None:
            color = [x/255 for x in tcutility.data.atom.color(atom.symbol)]
        sphereActor.GetProperty().SetColor(color)
        sphereActor.type = f'atom_{atom.symbol}'
        sphereActor.atom = atom
        sphereActor.SetUserTransform(self.transform)
        self.renderer.AddActor(sphereActor)

        if settings['atom']['draw_quadrants']:
            draw_disk(0, 0)
            draw_disk(-65, 0)
            draw_disk(0, -65)

    def draw_single_bond(self, p1, p2):
        p1, p2 = np.array(p1), np.array(p2)
        lineSource = vtkLineSource()
        lineSource.SetPoint1(p1)
        lineSource.SetPoint2(p2)

        tubeFilter = vtkTubeFilter()
        tubeFilter.source = lineSource
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(settings['bond']['radius'])
        tubeFilter.SetCapping(True)
        tubeFilter.SetNumberOfSides(20)

        tubeMapper = vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

        tubeActor = vtkActor()
        tubeActor.SetMapper(tubeMapper)
        tubeActor.GetProperty().SetColor(settings['bond']['color'])
        tubeActor.type = 'bond'
        # a1, a2 = [atom for atom in self.mol if atom.coords == p1][0], [atom for atom in self.mol if atom.coords == p2][0]
        tubeActor.atoms = p1.tolist(), p2.tolist()
        tubeActor.SetUserTransform(self.transform)
        self.renderer.AddActor(tubeActor)
        return tubeActor


    def draw_interaction_bond(self, p1, p2, spacing=.3):
        p1, p2 = np.array(p1), np.array(p2)
        l = np.linalg.norm(p1 - p2)
        direction = (p2 - p1) / l
        n_segments = int(l // spacing)
        # print(n_segments)
        margin = (l - n_segments * spacing) / 2
        for i in range(n_segments):
            pi = p1 + margin * direction + spacing * direction * i
            actor = self.draw_single_bond(pi, pi + spacing / 2 * direction)
            actor.atoms = p1.tolist(), p2.tolist()
            actor.type = 'intbond'


    def remove_bond(self, p1, p2):
        for actor in self.renderer.GetActors():
            if not hasattr(actor, 'type'):
                continue

            if actor.type not in ['bond', 'intbond']:
                continue

            print(actor.atoms, p1, p2)

            if p1 in actor.atoms and p2 in actor.atoms:
                self.renderer.RemoveActor(actor)

    def draw_isosurface(self, grid, isovalue=None, color=(1, 1, 0), opacity=.3, shininess=.0):
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
        actor.GetProperty().SetOpacity(opacity)
        actor.GetProperty().SetAmbient(.3)
        actor.GetProperty().SetDiffuse(1)
        actor.GetProperty().SetSpecular(shininess)
        actor.GetProperty().SetSpecularPower(70)
        actor.GetProperty().SetSpecularColor((1, 1, 1))
        actor.PickableOff()
        if isovalue < 0:
            actor.type = 'surfacem'
        else:
            actor.type = 'surfacep'
        actor.SetUserTransform(self.transform)

        self.renderer.AddActor(actor)

    def draw_dual_isosurface(self, grid, isovalue=None, colorm=(1, 0, 0), colorp=(0, 0, 1), opacity=None, shininess=None):
        if isovalue is None:
            isovalue = self.parent.parent.settings.get_value('Iso Surface', 'Iso Value')
        if opacity is None:
            opacity = self.parent.parent.settings.get_value('Iso Surface', 'Opacity')
        if shininess is None:
            shininess = self.parent.parent.settings.get_value('Iso Surface', 'Shininess')
        self.draw_isosurface(grid, isovalue,  color=colorp, opacity=opacity, shininess=shininess)
        self.draw_isosurface(grid, -isovalue, color=colorm, opacity=opacity, shininess=shininess)

    def draw_axes(self):
        for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
            lineSource = vtkLineSource()
            lineSource.SetPoint1((0, 0, 0))
            lineSource.SetPoint2(v)

            tubeFilter = vtkTubeFilter()
            tubeFilter.source = lineSource
            tubeFilter.SetInputConnection(lineSource.GetOutputPort())
            tubeFilter.SetRadius(0.01)
            tubeFilter.SetNumberOfSides(20)

            tubeMapper = vtkPolyDataMapper()
            tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

            tubeActor = vtkActor()
            tubeActor.SetMapper(tubeMapper)
            tubeActor.GetProperty().SetColor([255*i for i in v])
            # tubeActor.SetUserTransform(self.transform)
            self.renderer.AddActor(tubeActor)

    def draw_contour_lines(self, grid, isovalue=0, transform=None):
        depthArray = numpy_to_vtk(abs(grid.values).reshape(*grid.shape).ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        imdata = vtk.vtkImageData()
        imdata.SetDimensions(grid.shape)
        imdata.SetSpacing(grid.spacing)
        imdata.SetOrigin(grid.origin)
        imdata.GetPointData().SetScalars(depthArray)

        # transform = vtk.vtkTransform()
        planesource = vtk.vtkPlaneSource()
        max_size = 50

        planesource.SetResolution(2000, 2000)
        planesource.SetPoint1((max_size * 2, 0, 0))
        planesource.SetPoint2((0, 0, max_size * 2))
        planesource.SetOrigin((-max_size, 0, -max_size))
        # planesource.SetTransform(transform)

        filt = vtk.vtkTransformPolyDataFilter()
        filt.SetInputConnection(planesource.GetOutputPort())
        filt.SetTransform(transform)

        appendf = vtk.vtkAppendPolyData()
        appendf.AddInputConnection(filt.GetOutputPort())

        probe = vtk.vtkProbeFilter()
        probe.SetSourceData(imdata)
        probe.SetInputConnection(appendf.GetOutputPort())

        contour = vtk.vtkContourFilter()
        contour.SetInputConnection(probe.GetOutputPort())
        n = 20
        contour.GenerateValues(n, (0.001, 0.4))

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(contour.GetOutputPort())
        mapper.SetScalarRange(imdata.GetScalarRange())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(3)
        actor.SetUserTransform(self.transform)

        self.renderer.AddActor(actor)

    def remove_angle(self, a1, a2, a3):
        for act in self.renderer.GetActors():
            if not hasattr(act, 'type') or act.type != 'angle':
                continue

            if all(a in act.atoms for a in (a1, a2, a3)):
                self.renderer.RemoveActor(act)

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
        tubeActor.atoms = (a1, a2, a3)
        tubeActor.type = 'angle'
        self.renderer.AddActor(tubeActor)

    def toggle_angle(self, a1, a2, a3):
        for act in self.renderer.GetActors():
            if not hasattr(act, 'type') or act.type != 'angle':
                continue

            if all(a in act.atoms for a in (a1, a2, a3)):
                self.remove_angle(a1, a2, a3)
                break
        else:
            self.draw_angle(a1, a2, a3)

    def draw_text(self, text: str, location: str = 'bottom right', fontsize=24, color=(255, 255, 255)):
          textActor = vtkTextActor()
          textActor.SetInput(text)
          textActor.SetPosition(50, 50)
          textActor.UseBorderAlignOff()
          textActor.GetTextProperty().SetFontSize(fontsize)
          textActor.GetTextProperty().SetColor(color)
          self.renderer.AddActor2D(textActor)
          
    def use_perspective(self):
        self.use_parallel_projection = False
        camera = self.renderer.GetActiveCamera()
        camera.SetParallelProjection(False)

    def use_parallel(self):
        self.use_parallel_projection = True
        camera = self.renderer.GetActiveCamera()
        camera.SetParallelProjection(True)
        

class MoleculeWidget:
    def __new__(self, parent=None, headless=False):
        if headless or not has_qt:
            return _HeadlessMoleculeWidget()
        else:
            return _MoleculeWidget(parent)


class _HeadlessMoleculeWidget(vtk.vtkRenderWindowInteractor):
    def __init__(self):
        super().__init__()
        self.scenes = []

        self.molecule_renderers = []
        self.renWin = vtk.vtkRenderWindow()
        self.SetRenderWindow(self.renWin)
        # self.renWin.SetMultiSamples(4)
        self.renWin.BordersOn()
        self.interactor_style = vtk.vtkInteractorStyleTrackballCamera()
        # self.GetActiveCamera()
        self.SetInteractorStyle(self.interactor_style)
        self._base_ren = vtkRenderer()
        self._base_ren.SetBackground(1, 1, 1)
        self._base_ren.DrawOff()
        self.renWin.AddRenderer(self._base_ren)

        self.Initialize()
        # self.Start()

        self._recording_mouse = False
        self._mouse_pos = None

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

    def screenshot(self, path, scale=4):
        self.active_scene.screenshot(path, scale=scale)

    def screenshots(self, paths=None, directory=None, scale=4):
        if paths is None and directory is None:
            raise ValueError('You should give either the paths or directory argument')

        if paths is None:
            os.makedirs(directory, exist_ok=True)
            paths = [os.path.join(directory, f'scene{i}.png') for i in range(self.number_of_scenes)]

        for scene, path in zip(self.scenes, paths):
            scene.screenshot(path, scale=scale)


if has_qt:
    class _MoleculeWidget(QVTKRenderWindowInteractor):
        def __init__(self, parent):
            super().__init__()
            self.scenes = []

            self.parent = parent
            self.parent.settings.tabs['Atom'].connect(self._update_atoms)
            self.parent.settings.tabs['Bond'].connect(self._update_bonds)
            self.parent.settings.tabs['Iso Surface'].connect(self._update_isosurfaces)

            self.molecule_renderers = []
            self.renWin = self.GetRenderWindow()
            # self.renWin.SetMultiSamples(4)
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
            if actor.type.startswith('atom'):
                actor.SetScale([0.75, 0.75, 0.75])
                highlight = vtkSphereSource()
                highlight.SetPhiResolution(35)
                highlight.SetThetaResolution(45)
                highlight.SetCenter([0, 0, 0])
                radius = actor.GetMapper().GetInputConnection(0, 0).GetProducer().GetRadius()
                highlight.SetRadius(radius)

            elif actor.type in ['bond', 'intbond']:
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
            highlightActor.GetProperty().SetAmbient(1)
            highlightActor.PickableOff()
            highlightActor.SetUserTransform(self.active_scene.transform)
            highlightActor.type = 'highlight'
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
            if isinstance(xyz, str):
                if xyz.endswith('.xyz'):
                    mol = plams.Molecule(xyz)
                    orbs = False
                else:
                    res = tcutility.results.read(xyz)
                    orbs = pyfmo.orbitals.Orbitals(res.files['adf.rkf'])
                    mol = res.molecule.output
            elif isinstance(xyz, plams.Molecule):
                mol = xyz

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

        def screenshot(self, path, scale=4):
            self.active_scene.screenshot(path, scale=scale)

        def screenshots(self, paths=None, directory=None, scale=4):
            if paths is None and directory is None:
                raise ValueError('You should give either the paths or directory argument')

            if paths is None:
                os.makedirs(os.path.abspath(directory), exist_ok=True)
                paths = [os.path.join(directory, f'scene{i}.png') for i in range(self.number_of_scenes)]

            for scene, path in zip(self.scenes, paths):
                scene.screenshot(path, scale=scale)

        def _update_isosurfaces(self):
            opacity = self.parent.settings.get_value('Iso Surface', 'Opacity')
            shininess = self.parent.settings.get_value('Iso Surface', 'Shininess')

            for scene in self.scenes:
                actors = scene.renderer.GetActors()
                for actor in actors:
                    if not hasattr(actor, 'type'):
                        continue
                    if not actor.type.startswith('surface'):
                        continue

                    isovalue = self.parent.settings.get_value('Iso Surface', 'Iso Value')
                    if actor.type == 'surfacem':
                        isovalue *= -1

                    actor.GetMapper().GetInputConnection(0, 0).GetProducer().SetValue(0, isovalue)
                    actor.GetProperty().SetOpacity(opacity)
                    actor.GetProperty().SetSpecular(shininess)
                    actor.GetMapper().GetInputConnection(0, 0).GetProducer().Update()
            self.renWin.Render()

        def _update_atoms(self):
            draw_quadrants = self.parent.settings.get_value('Atom', 'Draw Quadrants')
            quadrant_width = self.parent.settings.get_value('Atom', 'Quadrant Width')
            size_ratio = self.parent.settings.get_value('Atom', 'Size Ratio')

            for scene in self.scenes:
                actors = scene.renderer.GetActors()
                for actor in actors:
                    if not hasattr(actor, 'type'):
                        continue
                    if not any(actor.type.startswith(typ) for typ in ['atom', 'quadrant']):
                        continue

                    symbol = actor.type.split('_')[1]
                    if actor.type.startswith('atom'):
                        radius = tcutility.data.atom.radius(symbol) * size_ratio
                    else:
                        radius = tcutility.data.atom.radius(symbol) * size_ratio + quadrant_width
                    actor.GetMapper().GetInputConnection(0, 0).GetProducer().SetRadius(radius)

                    if actor.type.startswith('quadrant'):
                        if draw_quadrants and quadrant_width > 0.0:
                            actor.VisibilityOn()
                        else:
                            actor.VisibilityOff()

            self.renWin.Render()

                # tab.add_spinbox('Radius', 0.08, minimum=0, decimals=3, suffix_text='Å')
                # tab.add_color('Color', (0, 0, 0))
                # tab.add_spinbox('Dashed Bond Radius', default=0.04, minimum=0, decimals=3, suffix_text='Å')

        def _update_bonds(self):
            radius = self.parent.settings.get_value('Bond', 'Radius')
            color = self.parent.settings.get_value('Bond', 'Color')
            dashbond_radius = self.parent.settings.get_value('Bond', 'Dashed Bond Radius')

            for scene in self.scenes:
                actors = scene.renderer.GetActors()
                for actor in actors:
                    if not hasattr(actor, 'type'):
                        continue
                    if not any(actor.type.startswith(typ) for typ in ['bond', 'intbond']):
                        continue

                    if actor.type == 'bond':
                        rad = radius
                    else:
                        rad = dashbond_radius
                    actor.GetMapper().GetInputConnection(0, 0).GetProducer().SetRadius(rad)
                    actor.GetProperty().SetColor(color)


            self.renWin.Render()

    class MoleculeWidgetKeyPressFilter(QtCore.QObject):
        def eventFilter(self, widget, event):
            if event.type() == QtCore.QEvent.KeyPress:
                scene = self.parent().scenes[self.parent().active_scene_index]
                if event.key() == QtCore.Qt.Key_Right:
                    self.parent().next_mol()
                if event.key() == QtCore.Qt.Key_Left:
                    self.parent().previous_mol()
                if event.key() == QtCore.Qt.Key_C and event.modifiers() == QtCore.Qt.ControlModifier:
                    print('copy')
                if event.key() == QtCore.Qt.Key_V and event.modifiers() == QtCore.Qt.ControlModifier:
                    print('paste')
                if event.key() == QtCore.Qt.Key_Space:
                    scene.reset_camera()
                    scene.post_draw()

                if event.key() == QtCore.Qt.Key_P:
                    atoms = [actor.atom for actor in self.parent().selected_actors if actor.type.startswith('atom')]
                    if len(atoms) == 3:
                        scene.toggle_angle(*atoms)
                        scene.post_draw()

                if event.key() == QtCore.Qt.Key_S:
                    filename, t = QtWidgets.QFileDialog.getSaveFileName()
                    scene.screenshot(filename.removesuffix('.png') + '.png')

                bond_selected = len(self.parent().selected_actors) == 1 and self.parent().selected_actors[0].type == 'bond'
                atoms2_selected = len(self.parent().selected_actors) == 2 and self.parent().selected_actors[0].type.startswith('atom') and self.parent().selected_actors[1].type.startswith('atom')
                atom_and_bond_selected = len(self.parent().selected_actors) == 2 and any(actor.type == 'bond' for actor in self.parent().selected_actors) and any(actor.type.startswith('atom') for actor in self.parent().selected_actors)
                if atoms2_selected:
                    actor1, actor2 = self.parent().selected_actors
                    if event.key() == QtCore.Qt.Key_1:
                        scene.draw_single_bond(actor1.atom.coords, actor2.atom.coords)
                        self.parent().Render()

                    if event.key() in [QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete]:
                        scene.remove_bond(actor1.atom.coords, actor2.atom.coords)
                        self.parent().Render()

                    if event.key() == QtCore.Qt.Key_I:
                        scene.draw_interaction_bond(actor1.atom.coords, actor2.atom.coords)
                        self.parent().Render()
                        # scene.dr


                if bond_selected:
                    actor = self.parent().selected_actors[0]
                    if event.key() in [QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete]:
                        print(actor.atoms)
                        scene.remove_bond(*actor.atoms)
                        self.parent().remove_highlight(actor)
                        self.parent().Render()


                if atom_and_bond_selected:
                    if event.key() == QtCore.Qt.Key_I:
                        actor1, actor2 = self.parent().selected_actors
                        if actor1.type != 'bond':
                            actor1, actor2 = actor2, actor1
                        p1, p2 = actor1.atoms
                        scene.draw_interaction_bond((np.array(p1) + np.array(p2))/2, actor2.atom.coords)
                        self.parent().Render()

                if event.key() == QtCore.Qt.Key_A and event.modifiers() == QtCore.Qt.ControlModifier:
                    atoms = [actor.atom for actor in self.parent().selected_actors if actor.type == 'atom']
                    bonds = [actor.atoms for actor in self.parent().selected_actors if actor.type == 'bond']
                    
                    T = scene.transform
                    if len(atoms) == 0:
                        # if the user selected a single bond we align it to the x-axis
                        if len(bonds) == 1:
                            a1, a2 = bonds[0]
                            c1, c2 = np.array(scene.transform.TransformPoint(a1.coords)), np.array(scene.transform.TransformPoint(a2.coords))
                            R = tcutility.geometry.vector_align_rotmat(c1 - c2, [1, 0, 0])
                            angles = tcutility.geometry.rotmat_to_angles(R)
                            T.Translate(-c1)
                            T.RotateX(angles[0] * 180 / np.pi)
                            T.RotateY(angles[1] * 180 / np.pi)
                            T.RotateZ(angles[2] * 180 / np.pi)

                            # self.parent().renWin.Render()
                            self.parent().Render()
                            scene.reset_camera()
                            scene.post_draw()

                        # if the user selected a single bond we align it to the x-axis
                        if len(bonds) == 2:
                            a1, a2, a3, a4 = [*bonds[0], *bonds[1]]
                            c1 = np.array(scene.transform.TransformPoint(a1.coords))
                            c2 = np.array(scene.transform.TransformPoint(a2.coords))

                            R = tcutility.geometry.vector_align_rotmat(c1 - c2, [1, 0, 0])
                            angles = tcutility.geometry.rotmat_to_angles(R)
                            T.Translate(-c1)
                            T.RotateX(angles[0] * 180 / np.pi)
                            T.RotateY(angles[1] * 180 / np.pi)
                            T.RotateZ(angles[2] * 180 / np.pi)


                            c1 = np.array(scene.transform.TransformPoint(a1.coords))
                            c2 = np.array(scene.transform.TransformPoint(a2.coords))
                            c3 = np.array(scene.transform.TransformPoint(a3.coords))
                            c4 = np.array(scene.transform.TransformPoint(a4.coords))

                            n = np.cross(c1 - c2, c3 - c4)

                            R = tcutility.geometry.vector_align_rotmat(n, [0, 1, 0])
                            angles = tcutility.geometry.rotmat_to_angles(R)
                            T.RotateX(angles[0] * 180 / np.pi)
                            T.RotateY(angles[1] * 180 / np.pi)
                            T.RotateZ(angles[2] * 180 / np.pi)

                            # self.parent().renWin.Render()
                            self.parent().Render()
                            scene.reset_camera()
                            scene.post_draw()

                    print('align!', atoms)


            return False
