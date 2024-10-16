#!/bin/bash
# Bamboo CI script to install vtkggdtools and run all benchmarks
# Note: this script should be run from the root of the git repository

# Debuggging:
set -e -o pipefail
echo "Loading modules:" $@

# Set up environment such that module files can be loaded
source /etc/profile.d/modules.sh
module purge
# Modules are supplied as arguments in the CI job:
module load $@

# Debuggging:
echo "Done loading modules"
set -x

# Export current PYTHONPATH so ASV benchmarks can import imas
export ASV_PYTHONPATH="$PYTHONPATH"

# Set up the testing venv
rm -rf venv  # Environment should be clean, but remove directory to be sure
python -m venv venv
source venv/bin/activate

# Install asv and vtkggd
pip install --upgrade pip setuptools wheel
pip install virtualenv .[benchmark]

# Copy previous results (if any)
mkdir -p /mnt/bamboo_deploy/vtkggdtools/benchmarks/results
mkdir -p .asv
cp -rf /mnt/bamboo_deploy/vtkggdtools/benchmarks/results .asv/

# Ensure there is a machine configuration
asv machine --yes

# Run ASV for the current commit, develop and main
asv run -v --skip-existing-successful HEAD^!
asv run -v --skip-existing-successful develop^!

# TODO: Turn on benchmarking for master branch as well. Currently there is an 
# issue with the master branch, where it is unable to find the 
# 'vtkggdtools.generator' module. As this is not used anymore since the IMASPy 
# migration. We can turn on benchmarking after a new master version.
# asv run --skip-existing-successful master^!

# Compare results
if [ `git rev-parse --abbrev-ref HEAD` == develop ]
then
    # asv compare master develop
    echo ""
else
    asv compare develop HEAD
fi

# Publish results
asv publish

# And persistently store them
cp -rf .asv/{results,html} /mnt/bamboo_deploy/vtkggdtools/benchmarks/

