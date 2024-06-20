from setuptools import setup

import versioneer

# Generate the plugin using jinja:
import os
from vtkggdtools.generator import generate
base_dir = os.path.dirname(__file__)
input_dir = os.path.join(os.path.join(base_dir, "vtkggdtools"), "plugins/templates")
output_dir = os.path.join(os.path.join(base_dir, "vtkggdtools"), "plugins")
generate(input_dir, output_dir)

setup(version=versioneer.get_version(), cmdclass=versioneer.get_cmdclass())
