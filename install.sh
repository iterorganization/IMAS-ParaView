# Source this script with an installation prefix.
# For example:
# 	$ source install.sh /home/ITER/panchuj/work/Debug/install
# vtkggdtools will be installed in /home/ITER/panchuj/work/Debug/install/lib/pythonx.y/site-packages
# the paraview plugin will go in /home/ITER/panchuj/work/Debug/install/bin/plugins

python_majmin_ver=$(python -c 'import sys; print(str(sys.version_info[0])+"."+str(sys.version_info[1]))')
install_dir=$1
if [ ! -d ${install_dir} ]; then
  install_dir=$(pwd)
else
  install_dir=$(cd ${install_dir}; pwd)
fi

library_install_dir="${install_dir}/lib/python${python_majmin_ver}/site-packages/vtkggdtools"
plugin_install_dir="${install_dir}/bin/plugins"

pip install --prefix="${install_dir}" --no-deps .

echo "Installed into ${library_install_dir}"

[ -d ${plugin_install_dir} ] || mkdir -p ${plugin_install_dir}
mv "${library_install_dir}/plugins/VTKGGDTools.py" "${plugin_install_dir}/"

export PYTHONPATH=${PYTHONPATH}:"${install_dir}/lib/python${python_majmin_ver}/site-packages"
export PYTHONPATH=$(python3 -c 'import os; from collections import OrderedDict; \
    l=os.environ["PYTHONPATH"].split(":"); print(":".join(OrderedDict.fromkeys(l)))' )

export PV_PLUGIN_PATH=${PV_PLUGIN_PATH}:"${plugin_install_dir}"
export PV_PLUGIN_PATH=$(python3 -c 'import os; from collections import OrderedDict; \
    l=os.environ["PV_PLUGIN_PATH"].split(":"); print(":".join(OrderedDict.fromkeys(l)))' )

unset python_majmin_ver
unset library_install_dir
unset plugin_install_dir
unset install_dir
