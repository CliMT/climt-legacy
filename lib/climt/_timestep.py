def leapfrog(afc, xinc, x, xold):

    alpha = 0.5
    #print 'AFC: ', afc
    # Step forward in time
    xnew = xold + xinc
    #xnew = x + xinc

    # Do asselin filtering
    d = afc*(xold + xnew - 2*x)
    xold_out = x + alpha*d
    #xold_out = x

    #update data
    x_out = xnew + (alpha-1)*d

    return x_out, xold_out


def asselin(afc, xnew, x, xold):


    # Do asselin filtering
    xold_out = x + afc*(xold + xnew - 2*x)

    #update data
    x_out = xnew

    return x_out, xold_out


c1 = 1.91666666667 #23/12
c2 = 1.33333333333 #16/12
c3 = 0.41666666667 #5/12

def adams_bashforth(xinc, xinc_old, xinc_older, x):

    return x + (23.*xinc - 16.*xinc_old + 5.*xinc_older)/12.
