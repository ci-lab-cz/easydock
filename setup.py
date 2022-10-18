import setuptools
from os import path
import moldock

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="moldock",
    version=moldock.__version__,
    author="Pavel Polishchuk",
    author_email="pavel_polishchuk@ukr.net",
    description="Python moldock to facilitate molecular docking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ci-lab-cz/docking-scripts",
    packages=['moldock'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    extras_require={
        'rdkit': ['rdkit>=2017.09'],
    },
    entry_points={'console_scripts':
                      ['vina_dock = moldock.vina_dock:main',
                       'gnina_dock = moldock.gnina_dock:main',
                       'get_sdf_from_dock_db = moldock.get_sdf_from_dock_db:main']}
)
