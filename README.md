# Alhazen-Ptolemy

The most up-to-date version of the associated paper, entitled *Solving the Alhazen-Ptolemy Problem: 
Determining Specular Points on Spherical Surfaces for Radiative Transfer of Titan's Seas*, can be 
found [here](https://iopscience.iop.org/article/10.3847/PSJ/abe4dd).

Implementations of the solutions for determining the specular point on a spherical surface given an 
arbitrary spatial configuration of source and observer (known as the Alhazen-Ptolemy problem).

`Alhazen-Ptolemy.c++` has three use modes: benchmarking, one-finite, and both-finite. 

The benchmarking usage is 

```
alhazen-ptolemy
```

The one-finite usage is 

```
alhazen-ptolemy <observer_angle> <c>
```

And the both-finite usage is

```
alhazen-ptolemy <observer_angle> <c> <b>
```

`observer_angle` corresponds to <i>&theta;</i><sub>obs</sub> from the paper (or 
<i>&theta;</i><sub>src</sub>, if the observer was rotated to &pi; / 2 instead). Note that the 
branch-deductions assume the point with the smaller radius has been rotated to &pi; / 2. 
`observer_angle` is specified in degrees as the angle from the positive <i>x</i>-axis. `b` is the 
ratio between the radius of the sphere and the radius of either the observer or the source -- 
whichever is larger -- and `c` is the other ratio. In the associated paper, `b` is the ratio between 
the radius of the sphere and the radius of the observer, i.e. <i>R</i><sub>sph</sub> / 
<i>R</i><sub>obs</sub>, and `c` is the ratio between the radius of the sphere and the radius of the 
source, i.e. <i>R</i><sub>sph</sub> / <i>R</i><sub>src</sub>.

`alhazenptolemy.py` has the same usage as `Alhazen-Ptolemy.c++` except it does not feature any 
benchmarking. This directory is also an installable python module. Installing with

```
python -m pip install Alhazen-Ptolemy
```

or

```
pip install Alhazen-Ptolemy
```

will make the `alhazenptolemy.py` module importable via

```
import alhazenptolemy
```


