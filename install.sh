# Source this script with an installation prefix.
# For example:
# 	$ source install.sh /home/ITER/panchuj/work/Debug/install
# vtkggdtools will be installed in /home/ITER/panchuj/work/Debug/install/lib/pythonx.y/site-packages

install_dir=$1
python_majmin_ver=$(python -c 'import sys; print(str(sys.version_info[0])+"."+str(sys.version_info[1]))')
pip install --prefix=${install_dir} --no-deps .

echo "Installed into ${install_dir}/lib/python${python_majmin_ver}/site-packages/vtkggdtools"

export PYTHONPATH=${PYTHONPATH}:${install_dir}/lib/python${python_majmin_ver}/site-packages
export PYTHONPATH=$(python3 -c 'import os; from collections import OrderedDict; \
    l=os.environ["PYTHONPATH"].split(":"); print(":".join(OrderedDict.fromkeys(l)))' )

export PV_PLUGIN_PATH=${PV_PLUGIN_PATH}:${install_dir}/lib/python${python_majmin_ver}/site-packages/vtkggdtools/plugins
export PV_PLUGIN_PATH=$(python3 -c 'import os; from collections import OrderedDict; \
    l=os.environ["PV_PLUGIN_PATH"].split(":"); print(":".join(OrderedDict.fromkeys(l)))' )
