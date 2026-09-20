---
icon: lucide/download
tags:
  - usage
---

# Installation

`dynamic_characterization` is a Python package, available on [PyPI](https://pypi.org/project/dynamic-characterization/) and on [conda](https://anaconda.org/diepers/dynamic_characterization). It needs Python 3.9 or newer.

## Installing with `uv` (recommended)

[uv](https://docs.astral.sh/uv/) is a fast, modern Python package manager written in Rust. It's significantly faster than pip and handles dependency resolution more reliably. Installing uv is a [one-liner](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer).

=== "New Project"

    Create a new project with uv:

    ```bash
    uv init my-dynchar-project
    cd my-dynchar-project
    uv add dynamic-characterization
    ```

    Run your scripts with:

    ```bash
    uv run python my_script.py
    ```

    Or directly get started in [JupyterLab](https://jupyter.org):

    ```bash
    uv run --with jupyter jupyter lab
    ```

    By default, `jupyter lab` will start the server at [http://localhost:8888/lab](http://localhost:8888/lab).

=== "New Environment"

    Create and activate a new virtual environment:

    ```bash
    uv venv .venv
    source .venv/bin/activate  # Linux/macOS
    # or: .venv\Scripts\activate  # Windows
    ```

    Install `dynamic_characterization`:

    ```bash
    uv pip install dynamic-characterization
    ```

=== "Existing Environment"

    If you already have an activated virtual environment:

    ```bash
    uv pip install dynamic-characterization
    ```

## Installing with `pip`

1. Install `python` from [the website](https://www.python.org/downloads/), your system package manager, or [Homebrew](https://docs.brew.sh/Homebrew-and-Python).

2. In a console or terminal window, create and activate a new virtual environment:

    ```bash
    python -m venv .venv
    source .venv/bin/activate  # Linux/macOS
    # or: .venv\Scripts\activate  # Windows
    ```

3. Install `dynamic_characterization`:

    ```bash
    pip install dynamic-characterization
    ```

You can also use pip to install useful libraries like `jupyterlab`.

## Installing with `conda` or `mamba`

!!! important "Prerequisites"

    1. A working installation of [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [`mamba`](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). If you are using `conda`, we recommend installing the [libmamba solver](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community).
    2. Basic knowledge of [Conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

1. Create a new Conda environment with `dynamic_characterization`:

    ```bash
    conda create -n dynchar -c conda-forge -c diepers dynamic_characterization
    ```

2. Activate the environment:

    ```bash
    conda activate dynchar
    ```

3. (Optional but recommended) You can also use conda to install useful libraries like `jupyterlab`:

    ```bash
    conda install -c conda-forge jupyterlab
    ```
