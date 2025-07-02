# TC-DTD
<p align="justify">
A glacier's surface debris layer significantly modifies underlying ice melt dependent on the thermal resistance of the debris cover, with thermal resistance being a function of debris thickness and effective thermal conductivity.
The most commonly used method to calculate effective thermal conductivity of supraglacial debris layers applies heat diffusion principles to a vertical array of temperature measurements through the supraglacial debris cover combined with an estimate of volumetric heat capacity of the debris as presented by <a href="https://iahs.info/uploads/dms/11685.89-97-264-Conway.pdf">Conway and Rasmussen (2000)</a>.
Application of this approach is only appropriate if the temperature data indicate that the system is predominantly conductive and, even in the case of a pure conductive system, if the method necessarily introduces numerical errors that can impact the derived values.
The sampling strategies used in published applications of this method vary in sensor precision and spatiotemporal temperature sampling strategies, hampering inter-site comparisons of the derived values and their usage at unmeasured sites.
Here we provide an interactive analysis tool to help researchers investigate the effect of the sampling interval on calculated sub-debris ice melt and plan future measurement campaigns. 
</p>

### Jupyter notebooks

This repository consists of two jupyter notebooks to interactively generate and analyze thermal diffusivity datasets:

**DTDmodel.ipynb** - In this jupyter notebook you can create your own debris thermal diffusivity data. You can select between a sinusoidal and a skewed sine forcing.

**DTDstd.ipynb** - In this jupyter notebook you can investigate the sampling interval dependency of debris thermal diffusivity calculations. 

### .py
The following python files contain the functions to run the jupyter notebooks:

    analysis.py
    equations.py
    functions.py
    plotfunction.py


### Python Libraries required:

    NumPy (v1.20.0)
    Matplotlib (v3.4.1)
    Statsmodels (v0.12.2)
    Pandas (v1.2.4)
    DateTime (v4.3)
    SciPy (v1.6.2)
    TQDM (v4.59.0)
    Natsort (v7.1.1)

The program code was tested with the stated library versions.
