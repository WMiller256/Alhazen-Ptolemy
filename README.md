# Alhazen-Ptolemy
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

`observer_angle` corresponds to <i>&theta;</i><sub>obs</sub> and is specified in degrees as the angle from the positive <i>x</i>-axis, 
`c` is the ratio between the radius of the sphere and the radius of the observer, i.e. <i>R</i><sub>sph</sub> / <i>R</i><sub>obs</sub>,
and `b` is the ratio between the radius of the sphere and the radius of the source, i.e. <i>R</i><sub>sph</sub> / <i>R</i><sub>src</sub> .

`Alhazen-Ptolemy.py` has the same usage as `Alhazen-Ptolemy.c++` except it does not feature any benchmarking.
