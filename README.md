# TCviewer
Quickly display molecules and surfaces with our TCviewer program. TCviewer is easy to use and yields beautiful figures.

## Installation <a name=installation></a>
*Currently, the installation can be done by cloning this repository and having all necessary dependencies.*

The following is for people who would like to install the repository themselves. For example, to edit and/or contribute code to the project.

First clone this repository:
```
git clone git@github.com:TheoChem-VU/TCviewer.git
```

Then move into the new directory and install the package:

```
cd TCviewer
python -m pip install --upgrade build 
python -m build 
python -m pip install -e .
```

To get new updates, simply run:
```
git pull
```
## Usage <a name=usage></a>
TCviewer exports the `Screen` object which can be used to quickly build a scene. All the heavy work, such as setting up the camera, building and coloring meshes, will be done in the background. Please see the examples section for some simple example scripts for you to try out.

## Examples <a name=examples></a>
To draw a molecule, simply run:

```python
from tcviewer import Screen
from scm import plams

molecule = 'path/to/molecule.xyz'  # molecules can be given as paths to xyz files
molecule = plams.Molecule('path/to/molecule.xyz')  # or as plams.Molecule objects
molecule.guess_bonds()
with Screen() as scr:
  scr.draw_molecule(molecule)
```

To load and draw the HOMO:

```python
from yutility import orbitals
from tcviewer import Screen, materials
import pathlib

rkf_dir = pathlib.Path(__file__).parents[0]/'data'/'NH3BH3'

# load orbitals and choose an MO to draw
orbs = orbitals.Orbitals(rkf_dir/'adf.rkf')
homo = orbs.mos['HOMO']
# generate a cube file
cub = homo.generate_orbital()

with Screen() as scr:
	scr.draw_cub(cub, 0.03, material=materials.orbital_matte)
```
