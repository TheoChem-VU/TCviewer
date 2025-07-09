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
    	'radius': 0.08,
    	'color': [0, 0, 0],
    	'dashed_bond_ratio': .5,
    },
}

if has_qt:
	class SettingsTab(QtWidgets.QWidget):
		def __init__(self):
			super().__init__()
			self.layout = QtWidgets.QGridLayout()
			self.setLayout(self.layout)
			self.widget_value_getters = {}
			self.widgets = {}

		def __enter__(self):
			return self

		def __exit__(self, *args):
			...

		def add(self, name: str = None, default: Any = None, range: Tuple[float] = None, options: List[str] = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)

		def add_checkbox(self, name: str = None, default: bool = True):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			widg = QtWidgets.QCheckBox()
			widg.setCheckState(QtCore.Qt.CheckState.Checked if default else QtCore.Qt.CheckState.Unchecked)
			self.layout.addWidget(widg, row, 1)
			self.widgets[name] = widg
			self.widget_value_getters[name] = widg.isChecked

		def add_spinbox(self, name: str = None, default: float = None, minimum: float = None, maximum: float = None, step: float = None, decimals: int = 2, suffix_text: str = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			widg = QtWidgets.QDoubleSpinBox()
			if default:
				widg.setValue(default)
			if minimum:
				widg.setMinimum(minimum)
			if maximum:
				widg.setMaximum(maximum)
			if decimals:
				widg.setDecimals(decimals)
				if step is None:
					step = 10**-decimals
			if step:
				widg.setSingleStep(step)

			self.layout.addWidget(widg, row, 1)

			if suffix_text:
				self.layout.addWidget(QtWidgets.QLabel(suffix_text), row, 2)

			self.widgets[name] = widg
			self.widget_value_getters[name] = widg.value

		def add_color(self, name: str = None, default: tuple[float] = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			sublayout = QtWidgets.QGridLayout()
			self.layout.addLayout(sublayout, row, 1)
			widgR = QtWidgets.QSpinBox()
			widgG = QtWidgets.QSpinBox()
			widgB = QtWidgets.QSpinBox()
			if default:
				widgR.setValue(default[0])
				widgG.setValue(default[1])
				widgB.setValue(default[2])

			widgR.setMinimum(0)
			widgG.setMinimum(0)
			widgB.setMinimum(0)
			widgR.setMaximum(255)
			widgG.setMaximum(255)
			widgB.setMaximum(255)

			sublayout.addWidget(QtWidgets.QLabel('R:'), 0, 1)
			sublayout.addWidget(widgR, 0, 2)
			sublayout.addWidget(QtWidgets.QLabel('G:'), 0, 3)
			sublayout.addWidget(widgG, 0, 4)
			sublayout.addWidget(QtWidgets.QLabel('B:'), 0, 5)
			sublayout.addWidget(widgB, 0, 6)

			class widg:
				def value(self):
					return (widgR.value()/255, widgG.value()/255, widgB.value()/255)

				def connect(self, func):
					widgR.valueChanged.connect(func)
					widgG.valueChanged.connect(func)
					widgB.valueChanged.connect(func)

			self.widgets[name] = widg()
			self.widget_value_getters[name] = self.widgets[name].value

		def add_slider(self, name: str = None, default: float = None, minimum: float = None, maximum: float = None, step: float = None, decimals: int = 2, suffix_text: str = None):
			row = self.layout.rowCount()
			self.layout.addWidget(QtWidgets.QLabel('    ' + name + ':'), row, 0)
			widg = QtWidgets.QDoubleSpinBox()

			if minimum:
				widg.setMinimum(minimum)
			if maximum:
				widg.setMaximum(maximum)
			if default:
				widg.setValue(default)
			if step:
				widg.setSingleStep(step)
			if decimals:
				widg.setDecimals(decimals)

			self.layout.addWidget(widg, row, 1)

			if suffix_text:
				self.layout.addWidget(QtWidgets.QLabel(suffix_text), row, 2)

			self.widgets[name] = widg
			self.widget_value_getters[name] = widg

		def get_value(self, name: str):
			return self.widget_value_getters[name]()

		def connect(self, func: callable):
			for widg in self.widgets.values():
				if hasattr(widg, 'valueChanged'):
					widg.valueChanged.connect(func)
				elif hasattr(widg, 'stateChanged'):
					widg.stateChanged.connect(func)
				else:
					widg.connect(func)


	class SettingsWidget(QtWidgets.QTabWidget):
		def __init__(self):
			super().__init__()
			self.tabs = {}

		def add_tab(self, tab_name: str):
			widg = SettingsTab()
			self.addTab(widg, tab_name)
			self.tabs[tab_name] = widg
			return widg

		def get_value(self, section: str, name: str):
			return self.tabs[section].get_value(name)


	class DefaultSettings(SettingsWidget):
		def __init__(self):
			super().__init__()

			with self.add_tab('Atom') as tab:
				tab.add_spinbox('Size Ratio', 1/2.4, minimum=0, step=0.05, decimals=3)
				tab.add_checkbox('Draw Quadrants', True)
				tab.add_spinbox('Quadrant Width', default=0.02, minimum=0, step=0.001, decimals=3, suffix_text='Å')

			with self.add_tab('Bond') as tab:
				tab.add_spinbox('Radius', 0.08, minimum=0, decimals=3, suffix_text='Å')
				tab.add_color('Color', (0, 0, 0))
				tab.add_spinbox('Dashed Bond Radius', default=0.04, minimum=0, decimals=3, suffix_text='Å')

			with self.add_tab('Iso Surface') as tab:
				tab.add_spinbox('Iso Value', 0.03, minimum=0, maximum=1, decimals=3)
				tab.add_spinbox('Opacity', 0.15, minimum=0, maximum=1, decimals=2, step=.1)
				tab.add_spinbox('Shininess', 0.0, minimum=0, maximum=1, decimals=2, step=.1)


if __name__ == '__main__':
	sett = SettingsWidget()
	with sett.add_tab('Atom') as tab:
		...
