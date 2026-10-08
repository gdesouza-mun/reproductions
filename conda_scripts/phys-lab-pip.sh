#!/bin/bash

ENV_NAME="physics-lab"
ENV_DIR="$HOME/.ENVS/$ENV_NAME"

mkdir -p "$HOME/.ENVS"

python3 -m venv "$ENV_DIR"

source "$ENV_DIR/bin/activate"

python -m pip install --upgrade pip

pip install \
    jupyterlab \
    numpy \
    scipy \
    pandas \
    matplotlib \
    uncertainties \
    odrpack

echo "Created environment: $ENV_DIR"
echo "Activate with: source $ENV_DIR/bin/activate"
echo "Run JupyterLab with: jupyter lab"
