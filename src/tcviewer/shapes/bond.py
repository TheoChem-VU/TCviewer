from typing import List, Union
import vtk
import tcmu
from tcviewer import shapes, register_setting
import numpy as np
from scm import plams


class Bond(shapes.EditableShape):
    def __init__(
        self, 
        renderer: vtk.vtkRenderer, 
        a1: plams.Atom, 
        a2: plams.Atom,
        name: str = None,

        radius: float = 0.1,
        outline_radius: float = 0.02,
        color: List[float] = [0, 0, 0],
        use_atom_colors: bool = False,

        opacity: float = 1,
        ambient: float = 0.65,
        diffuse: float = 0.5,
        specular: float = 0.5,
        specular_power: float = 5.0,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.a1 = a1
        self.a2 = a2
        self.name = name

        self.renderer = renderer

        # settings holds all properties related to the bond
        self.settings = {}
        self.settings['radius'] = radius
        self.settings['opacity'] = opacity
        self.settings['ambient'] = ambient
        self.settings['diffuse'] = diffuse
        self.settings['specular'] = specular
        self.settings['specular_power'] = specular_power
        self.settings['color'] = color
        self.settings['use_atom_colors'] = use_atom_colors
        self.settings['outline_radius'] = outline_radius

        c1 = np.array(a1.coords)
        c2 = np.array(a2.coords)

        self.line_source1 = vtk.vtkLineSource()
        self.line_source1.SetPoint1(c1)
        self.line_source1.SetPoint2((c1 + c2) / 2)

        self.line_source2 = vtk.vtkLineSource()
        self.line_source2.SetPoint1((c1 + c2) / 2)
        self.line_source2.SetPoint2(c2)

        self.tube_filter1 = vtk.vtkTubeFilter()
        self.tube_filter1.source = self.line_source1
        self.tube_filter1.SetInputConnection(self.line_source1.GetOutputPort())
        self.tube_filter1.SetRadius(self.settings['radius'])
        self.tube_filter1.SetNumberOfSides(10)

        self.tube_filter2 = vtk.vtkTubeFilter()
        self.tube_filter2.source = self.line_source2
        self.tube_filter2.SetInputConnection(self.line_source2.GetOutputPort())
        self.tube_filter2.SetRadius(self.settings['radius'])
        self.tube_filter2.SetNumberOfSides(10)

        tube_mapper1 = vtk.vtkPolyDataMapper()
        tube_mapper1.SetInputConnection(self.tube_filter1.GetOutputPort())

        tube_mapper2 = vtk.vtkPolyDataMapper()
        tube_mapper2.SetInputConnection(self.tube_filter2.GetOutputPort())

        self.tube_actor1 = vtk.vtkActor()
        self.tube_actor1.SetMapper(tube_mapper1)
        renderer.AddActor(self.tube_actor1)

        self.tube_actor2 = vtk.vtkActor()
        self.tube_actor2.SetMapper(tube_mapper2)
        renderer.AddActor(self.tube_actor2)

        # draw the outline
        self.line_source_outline = vtk.vtkLineSource()
        self.line_source_outline.SetPoint1(c1)
        self.line_source_outline.SetPoint2(c2)

        self.tube_filter_outline = vtk.vtkTubeFilter()
        self.tube_filter_outline.SetInputConnection(self.line_source_outline.GetOutputPort())
        self.tube_filter_outline.SetRadius(self.settings['radius'] + self.settings['outline_radius'])
        self.tube_filter_outline.SetNumberOfSides(10)
        self.tube_filter_outline.CappingOn()

        outlineMapper = vtk.vtkPolyDataMapper()
        outlineMapper.SetInputConnection(self.tube_filter_outline.GetOutputPort())

        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(outlineMapper)
        outlineActor.GetProperty().FrontfaceCullingOn()
        outlineActor.GetProperty().SetColor(0.0, 0.0, 0.0)

        renderer.AddActor(outlineActor)
        self._reset_properties()

        # this actor is used for drawing an outline around the atom
        self.outline_actor = vtk.vtkActor()
        self.outline_actor.SetMapper(outlineMapper)
        self.outline_actor.GetProperty().LightingOff()
        self.outline_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d('Magenta'))
        # self.outline_actor.SetPosition(self.coords)

    def __str__(self):
        return self.name if self.name is not None else repr(self)

    @register_setting('Geometry', 'Radius', 'radius', limits={'radius': (0, float('inf'))})
    def set_radius(self, radius: float = 0.1):
        self.settings['radius'] = radius
        self._reset_properties()

    @register_setting('Geometry', 'Outline Radius', 'outline_radius', limits={'outline_radius': (0, float('inf'))})
    def set_outline_radius(self, outline_radius: float = 0.02):
        self.settings['outline_radius'] = outline_radius
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

    @register_setting('Material', 'Use Atom Colors', 'use_atom_colors')
    def set_use_atom_colors(self, enable: bool = False):
        self.settings['use_atom_colors'] = enable
        self._reset_properties()

    def _reset_properties(self):
        self.tube_filter1.SetRadius(self.settings['radius'])
        self.tube_filter2.SetRadius(self.settings['radius'])
        self.tube_filter_outline.SetRadius(self.settings['radius'] + self.settings['outline_radius'])

        if self.settings['use_atom_colors']:
            self.tube_actor1.GetProperty().SetColor([x/255 for x in tcmu.data.atom.color(self.a1.symbol)])
            self.tube_actor2.GetProperty().SetColor([x/255 for x in tcmu.data.atom.color(self.a2.symbol)])
        else:
            self.tube_actor1.GetProperty().SetColor(self.settings['color'])
            self.tube_actor2.GetProperty().SetColor(self.settings['color'])

        self.tube_actor1.GetProperty().SetOpacity(self.settings['opacity'])
        self.tube_actor1.GetProperty().SetAmbient(self.settings['ambient'])
        self.tube_actor1.GetProperty().SetDiffuse(self.settings['diffuse'])
        self.tube_actor1.GetProperty().SetSpecular(self.settings['specular'])
        self.tube_actor1.GetProperty().SetSpecularPower(self.settings['specular_power'])

        self.tube_actor2.GetProperty().SetOpacity(self.settings['opacity'])
        self.tube_actor2.GetProperty().SetAmbient(self.settings['ambient'])
        self.tube_actor2.GetProperty().SetDiffuse(self.settings['diffuse'])
        self.tube_actor2.GetProperty().SetSpecular(self.settings['specular'])
        self.tube_actor2.GetProperty().SetSpecularPower(self.settings['specular_power'])
        # renderer.GetRenderWindow().Render()
