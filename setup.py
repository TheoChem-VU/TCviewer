from setuptools import setup

APP = ['src/tcviewer/screen.py']
DATA_FILES = []
OPTIONS = {'argv_emulation': True}

setup(
   app=APP,
   data_files=DATA_FILES,
   options={'py2app': OPTIONS},
   setup_requires=['py2app'],
)

if __name__ == "__main__":
    setup()
