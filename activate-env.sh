#!/bin/bash
# create and activate python environment for CgP-saarekkeet project
# script follows loosely https://google.github.io/styleguide/shell.xml
# about virtual envs: https://docs.python.org/3/library/venv.html


ENV_NAME=env-saarekkeet


#######################################
# Add created virtual env folder to 
# .gitignore file if not already present.
# Globals:
#   ENV_NAME
# Arguments:
#   None
# Returns:
#   None
#######################################
git_ignore_env() {
    if ! grep -q ${ENV_NAME} .gitignore; then
        echo $ENV_NAME >> .gitignore
    fi
}


#######################################
# Activate the virtual env.
# Globals:
#   ENV_NAME
# Arguments:
#   None
# Returns:
#   None
#######################################
activate_env() {
    . ${ENV_NAME}/bin/activate
    echo "* Virtual environment ${ENV_NAME} activated."
}
 

#######################################
# Install needed modules to created env.
#
# Installed modules:
# http://biopython.org/
# https://github.com/PyCQA/pycodestyle
# https://pypi.python.org/pypi/autopep8/
# https://www.pylint.org/
#
# Globals:
#   ENV_NAME
# Arguments:
#   None
# Returns:
#   None
#######################################
install_modules() {
    echo "* Installing modules for ${ENV_NAME}..."
    pip install biopython
    pip install pycodestyle 
    pip install --upgrade autopep8
    pip install pylint
}


#######################################
# Create virtual env using python3.
# Globals:
#   ENV_NAME
# Arguments:
#   None
# Returns:
#   None
#######################################
create_env() {
    echo "* Creating $ENV_NAME..."
    #python3 -m venv ${ENV_NAME}
    virtualenv ${ENV_NAME}
}


# the "main"
if [ -d ${ENV_NAME} ]; then
    activate_env
else
    create_env
   
    if [ -d ${ENV_NAME} ]; then
        activate_env
        install_modules
        git_ignore_env
    else
        echo "** Failed to create environment!"
    fi
fi
