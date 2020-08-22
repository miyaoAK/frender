# frender
Fractal renderer.

**NOTE**
This project is very buggy atm, please don't run it on your own PC. This repository merely serves as a backup atm.

I have given up on the boundary-tracing algorithm. This is simple because I misunderstood the M-set. It turns out that most of the bands that are formed by the escape-time
algorithm contain very few points (1-3 points). The version that I currently implemented would color the M-set correctly, if a few additions were made (rn it only colors the borders, and it doesn't skip already-traced bands). However, the fact that there are so many micro bands discouraged me to invest more time into this algorithm.

At a resolution of 1920x1080, using the border-tracing algorithm in combination with an optimized escape-time algorithm, I have to compute ~90.000 borders, which takes
approximately 2.5s. At the cost of visual fidelity, one could cut down on the number of borders, if X consecutive iteration bands share the same color.  

The optimized approach that uses periodicity- and cardoid checking renders the same image in ~0.25s. 
Neither approach uses multi-threading, but both use auto-vectorization.

Right now, I simply scan each pixel row-wise. Even a simple interleave algorithm would offer significant performance gains. I think, an interleave algorithm in combination
with multi-threading is able to render the M-set in less than 20ms. That's my goal at least.
