import Interpolation  # importing the necessary classes

x = [20, 21, 23, 24, 25, 27, 29]  # x values
y = [346, 362, 343, 339, 347, 346, 339]  # y values
x_value = 26  # x-axis value of desired approximation

# Lagrange Interpolation
lagrange = Interpolation.Lagrange(x, y)  # creating the lagrange interpolation class
lagrange.interpolate(x_value)  # interpolating to the given value
lagrange.plot()  # plotting the graph

# Newton's Divided Difference Polynomial Interpolation
nddp = Interpolation.NDDP(x, y)  # creating the Newton's divided differences polynomial interpolation class
nddp.interpolate(x_value)  # interpolating to the given value
nddp.plot()  # plotting the graph

# Direct Method of Polynomial Interpolation
directMethod = Interpolation.DirectMethod(x, y)  # creating the direct method of polynomial interpolation class
directMethod.interpolate(x_value)  # interpolating to the given value
directMethod.plot()  # plotting the graph
