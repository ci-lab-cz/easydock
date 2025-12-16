import setuptools
from os import path
import easydock

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="easydock",
    version=easydock.__version__,
    author="Pavel Polishchuk, Guzel Minibaeva, Aleksandra Ivanova",
    author_email="pavel_polishchuk@ukr.net",
    description="EasyDock Python module to facilitate molecular docking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ci-lab-cz/easydock",
    packages=['easydock'],
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
                      ['run_dock = easydock.run_dock_obsolete:main',
                       'easydock = easydock.run_dock:main',
                       'get_sdf_from_easydock = easydock.get_sdf_from_dock_db:main',
                       'make_clean_copy = easydock.make_clean_copy:main',
                       'easydock_plif = easydock.easydock_plif:main',]}
)
