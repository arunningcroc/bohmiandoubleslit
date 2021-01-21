## Bohmian particles in a double slit

This code shows a Bohmian particle in a double slit experiment.
To use it, you need numpy, matplotlib and cython. Pull the code and do

python setup.py build_ext --inplace

and it should work. Inspired heavily by
[this blog post of Thomas Bronzwaer](https://thomasbronzwaer.wordpress.com/2016/03/25/numerical-quantum-mechanics-the-time-dependent-schrodinger-equation-i/)

[I also wrote a blog post on this](https://runningcrocodile.fi/articles/bohmiandoubleslit.html). 

The walls are a bunch of gaussians next to each other to avoid the numerical problems of sharp-edge potentials.