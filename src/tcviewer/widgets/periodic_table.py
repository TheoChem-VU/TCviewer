from PySide6 import QtWidgets


elements_main = [
    ['H',  '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   'He'],
    ['Li', 'Be', '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   'B',  'C',  'N',  'O',  'F',  'Ne'],
    ['Na', 'Mg', '',   '',   '',   '',   '',   '',   '',   '',   '',   '',   'Al', 'Si', 'P',  'S',  'Cl', 'Ar'],
    ['K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    ['Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe'],
    ['Cs', 'Ba', 'La', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'],
    ['Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
]

elements_extra = [
    ['Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'],
    ['Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
]


class PeriodicTable(QtWidgets.QFrame):
    def __init__(self, parent) -> None:
        super().__init__(parent)
        self.make_buttons()

    def make_buttons(self):
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel('hey'), 2, 1)

        for j, row in enumerate(elements_main):
            for i, el in enumerate(row):
                if el == '':
                    continue

                btn = QtWidgets.QPushButton(el)
                layout.addWidget(btn, i, j)
