import Interpolation  # importing the necessary classes

x = [20, 21, 23, 24, 25, 27, 29]  # x values
y = [346, 362, 343, 339, 347, 346, 339]  # y values

# Lagrange Interpolation
lagrange = Interpolation.Lagrange(x, y)  # creating the lagrange interpolation class
lagrange.interpolate(26)  # interpolating to the given value
lagrange.plot()  # plotting the graph

# Newton's Divided Difference Polynomial Interpolation
nddp = Interpolation.NDDP(x, y)  # creating the Newton's divided differences polynomial interpolation class
nddp.interpolate(26)  # interpolating to the given value
nddp.plot()  # plotting the graph

# Direct Method of Polynomial Interpolation
directMethod = Interpolation.DirectMethod(x, y)  # creating the direct method of polynomial interpolation class
directMethod.interpolate(26)  # interpolating to the given value
directMethod.plot()  # plotting the graph
