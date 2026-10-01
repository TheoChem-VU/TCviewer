from PySide6 import QtWidgets


class Option(QtWidgets.QFrame):
    def __init__(self) -> None:
        super().__init__()
        self.layout = QtWidgets.HBoxLayout()
        self.setLayout(self.layout)


class Float(Option):
    def __init__(self, min_val: float = None, max_val: float = None, default_val: float = None) -> None:
        super().__init__()

        self.min_val = min_val
        self.max_val = max_val

        self.layout.addWidget(QtWidgets.QLineEdit())

