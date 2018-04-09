# wiener
recursion algorithm for Wiener Hopf equation, Simpson sideways recursion, least-squares
For Wiener filter, we usually need to solve the Wiener-Hopf equation.It is prohibitive to use the standard matrix operation to get 
the solution, specially for large number of filter coefficients. A series of recursion algorithm making using of the Toeplitz 
structure, can efficiently solve the equation.

We can find different filter coefficients corresponding to different the time lag between input data and output data. It provides 
a relatively inexpensive way find a filter of a fixed length having the optimum lag.Simpson sideways recursion is more efficient 
than simple brute-force method.

The names of the functions are

impuls  (impulse)
eureka  (general Toeplitz recursion)
peo     (auxiliary Toepliz recursion)
invtop  (inverse Toepliz matrix)
shape   (shaping filter)
spike   (spiking filter)
side    (Simpson sideways iteration)
spiker  (spiking filter, more efficient than spike)
shaper  (shaping filter for optimum positioning) 
