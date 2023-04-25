# Build docker for LightVegeManager with Singularity

In this project, we had to try the tool on distant servers, therefore we had to provide a working
environment. We choose Singularity to build our container.

The recipe `plantmix_env.def` is based on an existing singularity container with a pre-installed debian system with miniconda 3 and python 3.7.

The dependencies file `spec-openalea_py37.txt` stores a conda environment tested with LightVegeManager.

Building command:

```bash
sudo singularity build mycontainer.sif plantmix_env.def
```

Run command:

```bash
singularity run mycontainer.sif mypythonscript.py options
```
