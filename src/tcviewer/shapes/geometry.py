from tcviewer import shapes


class Geometry(shapes.EditableShape):
    def __init__(self) -> None:
        super().__init__()

        self.register('Geometry', 'Opacity', self.set_opacity, shapes.option.Float(0, 1, 1))

    def set_opacity(self, val: float):
        for actor in self.renderer.GetActors():
            actor.GetProperty().SetOpacity(val)

    def set_opacity(self, val: float):
        for actor in self.renderer.GetActors():
            actor.GetProperty().SetOpacity(val)
