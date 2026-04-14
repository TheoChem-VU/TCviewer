from typing import List, Union
import vtk
import tcmu
from tcviewer import shapes


class Atom(shapes.EditableShape):
    def __init__(self, renderer: vtk.vtkRenderer, atom_type: Union[str, int], coords: List[float]):
        self.renderer = renderer
        self.symbol = tcmu.data.atom.symbol(atom_type)
        self.coords = coords

        def draw_disk(rotatex, rotatey, name):
            circle = vtk.vtkRegularPolygonSource()
            circle.SetCenter([0, 0, 0])
            circle.SetRadius(tcmu.data.atom.radius(self.symbol) * 0.5 + 0.02)
            circle.SetNumberOfSides(50)
            circleMapper = vtk.vtkPolyDataMapper()
            circleMapper.SetInputConnection(circle.GetOutputPort())
            circleActor = vtk.vtkFollower()
            circleActor.SetCamera(renderer.GetActiveCamera())
            circleActor.GetProperty().SetOpacity(1)
            circleActor.SetPosition(coords)
            circleActor.SetMapper(circleMapper)
            circleActor.GetProperty().SetColor([0, 0, 0])
            circleActor.PickableOff()
            circleActor.RotateX(rotatex)
            circleActor.RotateY(rotatey)
            renderer.AddActor(circleActor)

        sphere = vtk.vtkSphereSource()
        sphere.SetPhiResolution(35)
        sphere.SetThetaResolution(45)
        sphere.SetRadius(tcmu.data.atom.radius(self.symbol) * 0.5)
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())
        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)
        sphereActor.SetPosition(coords)
        sphereActor.GetProperty().SetOpacity(1)
        sphereActor.GetProperty().SetAmbient(0.65)
        sphereActor.GetProperty().SetDiffuse(0.5)
        sphereActor.GetProperty().SetSpecular(0.5)
        sphereActor.GetProperty().SetSpecularPower(5.0)
        color = [x/255 for x in tcmu.data.atom.color(self.symbol)]
        sphereActor.GetProperty().SetColor(color)
        # sphereActor.type = f'atom_{atom.symbol}'
        # sphereActor.atom = atom
        # sphereActor.SetUserTransform(self.transform)
        renderer.AddActor(sphereActor)

        draw_disk(0, 0, 'outline')
        draw_disk(-65, 0, 'quadrant')
        draw_disk(0, -65, 'quadrant')
