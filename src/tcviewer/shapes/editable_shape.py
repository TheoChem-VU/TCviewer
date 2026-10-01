from typing import Callable
from tcviewer.shapes import option


class EditableShape:
    def __init__(self, outline_pass = None, **kwargs) -> None:
        self.edit_functions = []
        self.outline_pass = outline_pass

    def register(self, category: str, display_name: str, function: Callable, options: option.Option):
        self.edit_functions.append((category, display_name, function, options))

    def enable_outline(self):
        if hasattr(self, 'outline_pass'):
            self.outline_pass.add_actor(self.outline_actor)

    def disable_outline(self):
        if hasattr(self, 'outline_pass'):
            self.outline_pass.remove_actor(self.outline_actor)
