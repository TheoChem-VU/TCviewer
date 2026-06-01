from typing import Any, Tuple, List, get_type_hints
import inspect

try:
    from PySide6 import QtWidgets, QtCore
    has_qt = True
except ImportError:
    has_qt = False

settings = {
    'atom': {
        'size': 1/2.4,
        'draw_quadrants': True,
        'draw_outline': True,
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
            self.layout.setRowStretch(self.layout.rowCount(), 1)

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
                tab.add_checkbox('Draw Outline', True)
                tab.add_spinbox('Quadrant Width', default=0.02, minimum=0, step=0.001, decimals=3, suffix_text='Å')

            with self.add_tab('Bond') as tab:
                tab.add_spinbox('Radius', 0.08, minimum=0, decimals=3, suffix_text='Å')
                tab.add_color('Color', (0, 0, 0))
                tab.add_spinbox('Dashed Bond Radius', default=0.04, minimum=0, decimals=3, suffix_text='Å')

            with self.add_tab('Iso Surface') as tab:
                tab.add_spinbox('Iso Value', 0.03, minimum=0, maximum=1, decimals=3)
                tab.add_spinbox('Opacity', 0.15, minimum=0, maximum=1, decimals=2, step=.1)
                tab.add_spinbox('Shininess', 0.0, minimum=0, maximum=1, decimals=2, step=.1)
                tab.add_checkbox('Switch Phase Colors', False)

import random
_registered_settings = []


class SettingWidget(QtWidgets.QFrame):
    def __init__(self, objs, name, internal_variable, func, limits):
        super().__init__()
        self.setLayout(QtWidgets.QHBoxLayout())
        self.layout().addWidget(QtWidgets.QLabel(name))

        self.func = func
        self.objs = objs

        annots = func.__annotations__
        # internal_values = [ for obj in objs]
        internal_values = []
        for obj in objs:
            v = obj.settings[internal_variable]
            if not isinstance(v, list) and not isinstance(v, tuple):
                v = [v]

            internal_values.append(v)

        # check if the internal_values are all the same
        # if there is more than 1 different value we set
        # the value to None
        def all_same(vals):
            diffs = []
            for val in vals:
                if val not in diffs:
                    diffs.append(val)
            return len(diffs) == 1

        if all_same(internal_values):
            internal_values = internal_values[0]
        else:
            internal_values = [None] * len(internal_values[0])

        self.value_getters = []
        for val, (name, typ) in zip(internal_values, annots.items()):
            if typ is float:
                sb = QtWidgets.QDoubleSpinBox()
                if val is not None:
                    sb.setValue(val)
                self.layout().addWidget(sb)
                self.value_getters.append(sb.value)
                sb.valueChanged.connect(self.apply_setting)

            if typ is bool:
                sb = QtWidgets.QCheckBox()
                if val is not None:
                    sb.setChecked(val)
                self.layout().addWidget(sb)
                self.value_getters.append(sb.isChecked)
                sb.stateChanged.connect(self.apply_setting)


    def get_value(self):
        return [vg() for vg in self.value_getters]

    def apply_setting(self):
        v = self.get_value()
        renderers = set()
        for obj in self.objs:
            self.func(obj, *v)
            renderers.add(obj.renderer)

        for renderer in renderers:
            renderer.GetRenderWindow().Render()



class SettingsWidget(QtWidgets.QSplitter):
    def __init__(self):
        super().__init__()
        self.setOrientation(QtCore.Qt.Vertical)
        self.object_list = QtWidgets.QListWidget()
        self.object_list.setSelectionMode(QtWidgets.QListWidget.ExtendedSelection)
        self.object_list.itemSelectionChanged.connect(self.set_tabs)
        self.addWidget(self.object_list)
        self.settings_tabs = QtWidgets.QTabWidget()
        self.addWidget(self.settings_tabs)
        self.objects = {}

    def add_object(self, obj):
        item = QtWidgets.QListWidgetItem(str(obj))
        self.object_list.addItem(item)
        self.objects[str(item)] = obj
        # cls = obj.__class__.__name__
        # setts = [row for row in _registered_settings if row[-1] == cls]
        # if cls == 'Atom':
        #     setts[0][4](obj, random.random())
        # else:
        #     setts[0][4](obj, random.random()/10)

    def set_tabs(self):
        selected = self.object_list.selectedItems()

        self.settings_tabs.clear()
        objs = []
        for sel in selected:
            obj = self.objects[str(sel)]
            objs.append(obj)

        classes = set(obj.__class__.__name__ for obj in objs)

        # for obj in objs:
        for cls in classes:
            setts = [row for row in _registered_settings if row[-1] == cls]

            tabs = {}

            for sett in setts:
                tab, variable, internal_variable, limits, func, _ = sett
                if tab not in tabs:
                    frame = QtWidgets.QFrame()
                    frame.setLayout(QtWidgets.QVBoxLayout())
                    tabs[tab] = frame

                widg = SettingWidget(objs, variable, internal_variable, func, limits)
                tabs[tab].layout().addWidget(widg)

        for tab_name, tab_frame in tabs.items():
            self.settings_tabs.addTab(tab_frame, tab_name)
            tab_frame.layout().addStretch()


def register_setting(category: str, variable_name: str, internal_name: str, limits: dict = None):
    def inner_dec(func):
        cls = func.__qualname__.split('.')[-2]
        _registered_settings.append((category, variable_name, internal_name, limits, func, cls))
        return func
    return inner_dec


if __name__ == '__main__':
    sett = SettingsWidget()
    with sett.add_tab('Atom') as tab:
        ...
