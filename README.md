# An extended exercise in yak shaving

This repository contains some example code that shows how to get fast 
binomial convolutions if you really really care.  It's largely 
illustrative, and isn't particularly intended for Real Use.

The code was somewhat born out of frustration in trying to use R's 
gsDesign package on largish problems, but is frankly overkill by a 
couple of orders of magnitude.  I have almost certainly spent more time 
optimizing this than I will ever make back in using it.

If you want to play with this stuff you'll need FFTW3 handy.
