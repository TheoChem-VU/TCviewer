from typing import List, Union
import vtk
import tcmu
from tcviewer import shapes, register_setting


class Atom(shapes.EditableShape):
    def __init__(
        self, 
        renderer: vtk.vtkRenderer, 
        atom_type: Union[str, int], 
        coords: List[float],
        name: str = None,

        radius: float = None,
        scale: float = 0.5,
        quadrant_width: float = 0.02,
        quadrants_enabled: bool = True,

        opacity: float = 1,
        ambient: float = 0.65,
        diffuse: float = 0.5,
        specular: float = 0.5,
        specular_power: float = 5.0,
    ):
        self.renderer = renderer
        self.symbol = tcmu.data.atom.symbol(atom_type)
        self.coords = coords
        self.name = name

        # settings holds all properties related to the atom shape
        self.settings = {}

        if radius is None:
            self.settings['radius'] = tcmu.data.atom.radius(self.symbol)
        else:
            self.settings['radius'] = radius
        self.settings['scale'] = scale
        self.settings['quadrant_width'] = quadrant_width
        self.settings['quadrants_enabled'] = quadrants_enabled

        self.settings['opacity'] = opacity
        self.settings['ambient'] = ambient
        self.settings['diffuse'] = diffuse
        self.settings['specular'] = specular
        self.settings['specular_power'] = specular_power

        def draw_disk(rotatex, rotatey):
            circle = vtk.vtkRegularPolygonSource()
            circle.SetCenter([0, 0, 0])
            circle.SetNumberOfSides(50)
            circleMapper = vtk.vtkPolyDataMapper()
            circleMapper.SetInputConnection(circle.GetOutputPort())
            circleActor = vtk.vtkFollower()
            circleActor.SetCamera(renderer.GetActiveCamera())
            circleActor.SetPosition(coords)
            circleActor.SetMapper(circleMapper)
            circleActor.PickableOff()
            circleActor.RotateX(rotatex)
            circleActor.RotateY(rotatey)
            renderer.AddActor(circleActor)

            return circle, circleActor

        self.sphere = vtk.vtkSphereSource()
        self.sphere.SetPhiResolution(35)
        self.sphere.SetThetaResolution(45)
        self.sphere.SetRadius(self.settings['radius'] * self.settings['scale'])
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInputConnection(self.sphere.GetOutputPort())
        self.sphere_actor = vtk.vtkActor()
        self.sphere_actor.SetMapper(sphereMapper)
        self.sphere_actor.SetPosition(coords)
        self.settings['color'] = [x/255 for x in tcmu.data.atom.color(self.symbol)]
        self.settings['outline_color'] = [0, 0, 0]
        renderer.AddActor(self.sphere_actor)

        self.outline_disk, self.outline_disk_actor = draw_disk(0, 0)
        self.quadrant_disk1, self.quadrant_disk1_actor = draw_disk(-65, 0)
        self.quadrant_disk2, self.quadrant_disk2_actor = draw_disk(0, -65)

        self._reset_properties()

    def __str__(self):
        return self.name if self.name is not None else repr(self)

    @register_setting('Geometry', 'Radius', 'radius', limits={'radius': (0, float('inf'))})
    def set_radius(self, radius: float = None):
        if radius is None:
            self.settings['radius'] = tcmu.data.atom.radius(self.symbol)
        else:
            self.settings['radius'] = radius
        self._reset_properties()

    @register_setting('Geometry', 'Scale', 'scale', limits={'scale': (0, float('inf'))})
    def set_scale(self, scale: float = 0.5):
        self.settings['scale'] = scale
        self._reset_properties()

    @register_setting('Geometry', 'Quadrant Width', 'quadrant_width', limits={'width': (0, float('inf'))})
    def set_quadrant_width(self, width: float = 0.02):
        self.settings['quadrant_width'] = width
        self._reset_properties()

    @register_setting('Geometry', 'Draw Quadrants', 'quadrants_enabled')
    def set_quadrants_enabled(self, enable: bool = True):
        self.settings['quadrants_enabled'] = enable
        self._reset_properties()

    @register_setting('Material', 'Opacity', 'opacity', limits={'opacity': (0, 1)})
    def set_opacity(self, opacity: float = 1):
        self.settings['opacity'] = opacity
        self._reset_properties()

    @register_setting('Material', 'Ambient', 'ambient', limits={'ambient': (0, 1)})
    def set_ambient(self, ambient: float = 0.65):
        self.settings['ambient'] = ambient
        self._reset_properties()

    @register_setting('Material', 'Diffuse', 'diffuse', limits={'diffuse': (0, 1)})
    def set_diffuse(self, diffuse: float = 0.5):
        self.settings['diffuse'] = diffuse
        self._reset_properties()

    @register_setting('Material', 'Specular', 'specular', limits={'specular': (0, 1)})
    def set_specular(self, specular: float = 0.5):
        self.settings['specular'] = specular
        self._reset_properties()

    @register_setting('Material', 'Specular Power', 'specular_power', limits={'specular_power': (0, float('inf'))})
    def set_specular_power(self, specular_power: float = 5.0):
        self.settings['specular_power'] = specular_power
        self._reset_properties()

    @register_setting('Material', 'Color', 'color', limits={'R': (0, 1), 'G': (0, 1), 'B': (0, 1)})
    def set_color(self, R: float = None, G: float = None, B: float = None):
        if R is not None:
            self.settings['color'][0] = R
        if G is not None:
            self.settings['color'][1] = G
        if B is not None:
            self.settings['color'][2] = B
        self._reset_properties()

    @register_setting('Material', 'Outline Color', 'outline_color', limits={'R': (0, 1), 'G': (0, 1), 'B': (0, 1)})
    def set_outline_color(self, R: float = None, G: float = None, B: float = None):
        if R is not None:
            self.settings['outline_color'][0] = R
        if G is not None:
            self.settings['outline_color'][1] = G
        if B is not None:
            self.settings['outline_color'][2] = B
        self._reset_properties()

    def _reset_properties(self):
        self.sphere.SetRadius(self.settings['radius'] * self.settings['scale'])
        if self.settings['quadrants_enabled']:
            self.outline_disk.SetRadius(self.settings['radius'] * self.settings['scale'] + self.settings['quadrant_width'])
            self.quadrant_disk1.SetRadius(self.settings['radius'] * self.settings['scale'] + self.settings['quadrant_width'])
            self.quadrant_disk2.SetRadius(self.settings['radius'] * self.settings['scale'] + self.settings['quadrant_width'])
        else:
            self.outline_disk.SetRadius(0)
            self.quadrant_disk1.SetRadius(0)
            self.quadrant_disk2.SetRadius(0)

        for actor in self.renderer.GetActors():
            actor.GetProperty().SetOpacity(self.settings['opacity'])
            actor.GetProperty().SetAmbient(self.settings['ambient'])
            actor.GetProperty().SetDiffuse(self.settings['diffuse'])
            actor.GetProperty().SetSpecular(self.settings['specular'])
            actor.GetProperty().SetSpecularPower(self.settings['specular_power'])

        self.sphere_actor.GetProperty().SetColor(self.settings['color'])
        self.outline_disk_actor.GetProperty().SetColor(self.settings['outline_color'])
        self.quadrant_disk1_actor.GetProperty().SetColor(self.settings['outline_color'])
        self.quadrant_disk2_actor.GetProperty().SetColor(self.settings['outline_color'])



