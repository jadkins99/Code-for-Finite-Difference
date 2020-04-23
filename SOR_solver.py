

def SOR_solver(ud,bp,ap,omega,eps):

    loops = 0
    error = 0.000001

    num = len(ud)
    
    while(error > eps):
        error = 0

        for np in range(1,num-1):
            y = bp[np]*ap*(ud[np-1]+ud[np+1])/(1+2*ap)
            y = ud[np]+omega*(y-ud[np])
            error += (ud[np] - y)**2
            ud[np] = y
        loops += 1
    return ud,loops