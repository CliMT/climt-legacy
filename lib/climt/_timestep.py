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

