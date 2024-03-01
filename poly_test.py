from poly import Polynomial, Basis, InterpolationPoly
from curve import Scalar

def poly_test():
    vals = [1, 2, 3, 4]
    vals_scalar = [Scalar(int(x)) for x in vals]
    roots_of_unity = Scalar.roots_of_unity(4)

    poly_lag = Polynomial(vals_scalar, Basis.LAGRANGE)
    poly_coeff = poly_lag.ifft()
    points = roots_of_unity + [Scalar(2), Scalar(3), Scalar(4)]
    for i in range(len(points)):
      point = points[i]
      eval_lag = poly_lag.barycentric_eval(point)
      coeff_eval = poly_coeff.coeff_eval(point)
      assert eval_lag == coeff_eval

    quo = poly_coeff / Scalar(2)
    print("quo: ", quo.values)


def lagrange_poly_test():
    x_vals = [2, 3, 6]
    y_vals = [3, 4, 8]
    X = [Scalar(x) for x in x_vals]
    Y = [Scalar(y) for y in y_vals]
    lagrange_poly = InterpolationPoly(X, Y)
    v_poly = lagrange_poly.vanishing_poly()
    assert v_poly.values[0] == -Scalar(36)
    assert v_poly.values[1] == 36
    assert v_poly.values[2] == -Scalar(11)
    assert v_poly.values[3] == 1
    print("v_poly.values: ", v_poly.values)
    d_poly = lagrange_poly.vanishing_poly_diff()
    assert d_poly.values[0] == 36
    assert d_poly.values[1] == -Scalar(22)
    assert d_poly.values[2] == 3
    print("d_poly.values: ", d_poly.values)
    l_poly = lagrange_poly.lagrange_poly(0)
    print("l_poly.values: ", l_poly.values)
    poly = lagrange_poly.poly()
    print("poly.values: ", poly.values)
    for i in range(len(X)):
      assert poly.coeff_eval(X[i]) == Y[i]

if __name__ == "__main__":
    print("===========> Beginning Polynomial test <===========")
    poly_test()
    lagrange_poly_test()
    print("===========> Polynomial Test success <===========")
