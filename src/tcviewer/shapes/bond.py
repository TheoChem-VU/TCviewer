from typing import List, Union
import vtk
import tcmu
from tcviewer import shapes, register_setting
import numpy as np
from scm import plams


class Bond(shapes.EditableShape):
    def __init__(self, 
            renderer: vtk.vtkRenderer, 
            a1: plams.Atom, 
            a2: plams.Atom,
            name: str = None,

            radius: float = 0.05,
            color: List[float] = (0, 0, 0),

            opacity: float = 1,
            ambient: float = 0.65,
            diffuse: float = 0.5,
            specular: float = 0.5,
            specular_power: float = 5.0,
            ):
        self.renderer = renderer
        self.a1 = a1
        self.a2 = a2
        self.name = name

        # settings holds all properties related to the atom shape
        self.settings = {}
        self.settings['radius'] = radius
        self.settings['opacity'] = opacity
        self.settings['ambient'] = ambient
        self.settings['diffuse'] = diffuse
        self.settings['specular'] = specular
        self.settings['specular_power'] = specular_power
        self.settings['color'] = color


        a1, a2 = np.array(a1.coords), np.array(a2.coords)
        self.line_source = vtk.vtkLineSource()
        self.line_source.SetPoint1(a1)
        self.line_source.SetPoint2(a2)

        self.tube_filter = vtk.vtkTubeFilter()
        self.tube_filter.source = self.line_source
        self.tube_filter.SetInputConnection(self.line_source.GetOutputPort())
        self.tube_filter.SetRadius(self.settings['radius'])

        self.tube_filter.SetNumberOfSides(10)

        tube_mapper = vtk.vtkPolyDataMapper()
        tube_mapper.SetInputConnection(self.tube_filter.GetOutputPort())

        self.tube_actor = vtk.vtkActor()
        self.tube_actor.SetMapper(tube_mapper)
        renderer.AddActor(self.tube_actor)

        self._reset_properties()

    def __str__(self):
        return self.name if self.name is not None else repr(self)

    @register_setting('Geometry', 'Radius', 'radius', limits={'radius': (0, float('inf'))})
    def set_radius(self, radius: float = 0.05):
        self.settings['radius'] = radius
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

    def _reset_properties(self):
        self.tube_filter.SetRadius(self.settings['radius'])
        self.tube_actor.GetProperty().SetColor(self.settings['color'])
        self.tube_actor.GetProperty().SetOpacity(self.settings['opacity'])
        self.tube_actor.GetProperty().SetAmbient(self.settings['ambient'])
        self.tube_actor.GetProperty().SetDiffuse(self.settings['diffuse'])
        self.tube_actor.GetProperty().SetSpecular(self.settings['specular'])
        self.tube_actor.GetProperty().SetSpecularPower(self.settings['specular_power'])
