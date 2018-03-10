import numpy as np
import matplotlib.pyplot as plt


def BasisFunc(s, j_spline, knots, spline_order):

    s_min = knots[j_spline]
    s_max = knots[j_spline+1]

    if spline_order == 0:
        if (s <= s_max and s >= s_min ):
            B = 1
        else:
            B = 0
    else:
        d1 = knots[j_spline+spline_order]-knots[j_spline]

        if d1 == 0:
            mult1 = 0
        else:
            mult1 = (s-knots[j_spline])/d1

        d2 = knots[j_spline+spline_order+1]-knots[j_spline+1]

        if d2 == 0:
            mult2 = 0
        else:
            mult2 = (knots[j_spline+spline_order+1] - s)/d2

        B = mult1 * BasisFunc(s, j_spline, knots, spline_order - 1) + mult2 * BasisFunc(s, j_spline + 1, knots, spline_order - 1)

    return B


# Used to compute the first derivative of B-splines
def BasisFunc_der(s, j_spline, knots, spline_order):

	d1 = knots[j_spline+spline_order]-knots[j_spline]

	if d1 == 0:
		mult1 = 0
	else:
		mult1 = spline_order/d1

	d2 = knots[j_spline+spline_order+1]-knots[j_spline+1]

	if d2 == 0:
		mult2 = 0
	else:
		mult2 = spline_order/d2

	B_der = mult1*BasisFunc(s, j_spline, knots, spline_order-1) - mult2*BasisFunc(s, j_spline+1, knots, spline_order-1)

	return B_der
