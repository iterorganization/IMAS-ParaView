# Author: Panchumarti Jaswant EXT
# Description: Install vtkggdtools on any system with pip.

import setuptools
import os
import versioneer

from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install

from vtkggdtools._version import __version__
from vtkggdtools.generator import generate

base_dir = os.path.dirname(__file__)
input_dir = os.path.join(os.path.join(base_dir, "vtkggdtools"), "plugins/templates")
output_dir = os.path.join(os.path.join(base_dir, "vtkggdtools"), "plugins")
generate(input_dir, output_dir)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vtkggdtools",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    # setuptools_git_versioning={
    #     "enabled": True,
    # },
    # setup_requires=[
    #     "setuptools-git-versioning"
    # ],
    # version_config={
    #     #"version_callback": __version__,
    #     "template": "{tag}",
    #     "dirty_template": "{tag}.{ccount}.{sha}",
    # },
    author="Panchumarti Jaswant EXT",
    author_email="jaswant.panchumarti@iter.org",
    description="Tools to expose GGD readers/writers as ParaView plugins.",
    long_description=long_description,
    url="https://git.iter.org/projects/IMEX/repos/ggd-vtk/browse",
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords="VTK-GGD, VTK, Tools",
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.8",
    install_requires=["jinja2 >= 3.0.3", "numpy >= 1.19.4", "vtk >= 9.1.0"],
)
