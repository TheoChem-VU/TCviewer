from typing import Any, Tuple, List

try:
    from PySide6 import *
    has_qt = True
except ImportError:
    has_qt = False

settings = {
    'atom': {
        'size': 1/2.4,
        'draw_quadrants': True,
        'quadrant_follow_camera': True,
        'quadrant_width': 0.02,
    },
    'bond': {
    	'radius': 0.075,
    	'color': [0, 0, 0],
    },
}

if has_qt:
	class SettingsTab(QtWidgets.QWidget):
		def __init__(self):
			super().__init__()
			self.layout = QtWidgets.QGridLayout()
			self.setLayout(self.layout)

		def __enter__(self):
			return self

		def __exit__(self, *args):
			...

		def add(self, section: str = None, name: str = None, default: Any = None, range: Tuple[float] = None, options: List[str] = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)

		def add_checkbox(self, section: str = None, name: str = None, default: bool = True):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			widg = QtWidgets.QCheckBox()
			widg.setCheckState(QtCore.Qt.CheckState.Checked if default else QtCore.Qt.CheckState.Unchecked)
			self.layout.addWidget(widg, row, 1)

		def add_spinbox(self, section: str = None, name: str = None, default: float = None, minimum: float = None, maximum: float = None, step: float = None, decimals: int = 2, suffix_text: str = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			widg = QtWidgets.QDoubleSpinBox()
			if default:
				widg.setValue(default)
			if minimum:
				widg.setMinimum(minimum)
			if maximum:
				widg.setMaximum(maximum)
			if step:
				widg.setSingleStep(step)
			if decimals:
				widg.setDecimals(decimals)

			self.layout.addWidget(widg, row, 1)

			if suffix_text:
				self.layout.addWidget(QtWidgets.QLabel(suffix_text), row, 2)



	class SettingsWidget(QtWidgets.QTabWidget):
		def __init__(self):
			super().__init__()

		def add_tab(self, tab_name: str):
			widg = SettingsTab()
			self.addTab(widg, tab_name)
			return widg


	class DefaultSettings(SettingsWidget):
		def __init__(self):
			super().__init__()

			with self.add_tab('Molecule') as tab:
				tab.add_spinbox('Atom', 'Size Ratio', 1/2.4, minimum=0, step=0.05, decimals=3)
				tab.add_checkbox('Atom', 'Draw Quadrants', True)
				tab.add_spinbox('Atom', 'Quadrant Width', default=0.02, minimum=0, step=0.001, decimals=3, suffix_text='Å')
			# self.add_tab('Bond')


if __name__ == '__main__':
	sett = SettingsWidget()
	with sett.add_tab('Atom') as tab:
		...
