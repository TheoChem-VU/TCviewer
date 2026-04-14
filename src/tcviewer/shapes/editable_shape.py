from typing import Callable
from tcviewer.shapes import option


class EditableShape:
    def __init__(self) -> None:
        self.edit_functions = []

    def register(self, category: str, display_name: str, function: Callable, options: option.Option):
        self.edit_functions.append((category, display_name, function, options))
