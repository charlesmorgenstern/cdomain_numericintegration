# cdomain_numericintegration
Numeric integration methods for a "c" shaped domain with Green's Theorem and Monte Carlo

The domain is defined as four arrays of xy points which define paths that form the boundary.

getpaths(n1,n2,n3,n4) loads the paths. ni define the  number of points used for the ith path.

plotdomain()  plots the domain and shows the indexing for the paths.

greensintegral(n1,n2,n3,n4) uses Green's theorem to approximate the area of the domain by turning the double integral into a bounday integral. ni are the number of points used for numeric integration on the ith path. Relative error is given.

montecarlo_cdom(n) uses Monte Carlo to approximate the area with n random points. The points are plotted on the domain showing which points are in and out. Relative error is given.

montecarlo_noplot(n) is the same but does not produce a plot.
