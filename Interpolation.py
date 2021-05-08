import numpy as np  # for required calculations
import matplotlib.pyplot as plt  # for plotting the graphs


class Interpolation():

    def __init__(self, x, y):  # default constructor
        self.x = x  # will be used in calculations and plotting graph
        self.y = y  # will be used in calculations and plotting graph

    def interpolate(self, xn):  # will be overridden in sub classes
        print("Not defined interpolation")

    def plot(self):  # all the sub classes will use this basic scheme to plot a graph
        x_axis = np.msort(self.x)  # x values must be sorted
        y_axis = np.msort(self.y)  # y values must be sorted

        plt.plot(x_axis, y_axis)  # plotting the graph according to x and y values
        plt.title(self.__class__.__name__ + " Interpolation")  # adding the title of interpolation
        plt.xlabel("x")  # labeling the x-axis
        plt.ylabel("y")  # labeling the y-axis
        plt.savefig(self.__class__.__name__ + "_interpolation.png", dpi=500)  # saving the graph before showing it
        plt.show()  # showing the graph


class Lagrange(Interpolation):

    def interpolate(self, xn):  # overridden interpolate method to calculate lagrange interpolation
        x = np.array(self.x)  # getting the x values as numpy array
        y = np.array(self.y)  # getting the y values as numpy array

        m = len(x)  # length of values will be used in calculations
        n = m - 1  # basically, n-valued interpolation means (n-1)th degree polynomial

        yn = 0  # for calculating f(xn) = p0(xn)*y0 + p1(xn)*y1 +...+ pi(xn)*yi = yn
        for i in range(n + 1):
            p = 1
            for j in range(n + 1):
                if j != i:
                    p *= (xn - x[j]) / (x[i] - x[j])  # pi(xn)
            yn += y[i] * p  # pi(xn)*yn

        self.x = np.append(self.x, xn)  # saving the xn to include into the graph
        self.y = np.append(self.y, yn)  # saving the yn to include into the graph

        print("Lagrange Interpolation has been done, for x=%2.f, y=%f" % (xn, yn))


class NDDP(Interpolation):

    def interpolate(self, xn):  # overridden interpolate method to calculate nddp interpolation

        x = np.array(self.x)  # getting the x values as numpy array
        y = np.array(self.y)  # getting the y values as numpy array

        c_vector = self.get_ndd_coefficients(x, y)  # to get ndd coefficients

        final_pol = np.polynomial.Polynomial([0.])  # target polynomial
        n = c_vector.shape[0]  # to get number of coefficients
        for i in range(n):
            p = np.polynomial.Polynomial([1.])
            for j in range(i):
                p_temp = np.polynomial.Polynomial([-x[j], 1.])  # x - xj
                p = np.polymul(p, p_temp)  # multiplying dummy with expression
            p *= c_vector[i]  # applying coefficient
            final_pol = np.polyadd(final_pol, p)  # adding to target polynomial
        p = np.flip(final_pol[0].coef)  # coefficients

        yn = np.polyval(p, xn)  # adding to final polynomial

        self.x = np.append(self.x, xn)  # saving the xn to include into the graph
        self.y = np.append(self.y, yn)  # saving the yn to include into the graph

        print("Newton's Divided Differences Polynomial Interpolation has been done, for x=%2.f, y=%f" % (xn, yn))

    def get_ndd_coefficients(self, x, y):
        n = np.shape(y)[0]
        pyramid = np.zeros([n, n])  # creating a square matrix to hold pyramid
        pyramid[::, 0] = y  # first column is y
        for j in range(1, n):
            for i in range(n - j):
                pyramid[i][j] = (pyramid[i + 1][j - 1] - pyramid[i][j - 1]) / (x[i + j] - x[i])
        return pyramid[0]


class DirectMethod(Interpolation):

    def interpolate(self, xn):
        n = len(self.x)  # how many inputs we'll work on -> (n-1)th degree polynomial

        x = np.array(self.x)
        y = np.array(self.y).reshape(n, 1)

        c_matrix = np.column_stack([x ** e for e in range(0, n)])  # to get coefficients matrix
        s_matrix = np.linalg.solve(c_matrix, y)  # we are solving CS = Y with respect to S

        yn = self.get_value(s_matrix, xn)  # calculating the value of y using b0, b1,..., bn

        self.x = np.append(self.x, xn)  # saving the xn to include into the graph
        self.y = np.append(self.y, yn)  # saving the yn to include into the graph

        print("Direct Method Interpolation has been done, for x=%2.f, y=%f" % (xn, yn))

    def get_value(self, s_matrix, xn):  # direct method coefficients b0, b1,..., bn
        yn = 0  # will be used to calculate f(xn) = b0 + b1*x + b2*x^2 +...+ bn*x^n = yn
        for i in range(len(s_matrix)):
            yn += (s_matrix[i]) * (pow(xn, i))  # to add every single different term into yn to calculate f(xn)
        return yn
