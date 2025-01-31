############################################################
#                                                          #
#               General Utility Functions                  #
#               Author: Wyatt Giroux                       #
#               Date: 9/28/23                              #
#                                                          #
############################################################
def bisection(Lb, Rb, RHSfunc, LHS, inputs=[], direction=True, iterLim=5000, errAllow=1e-10, verbose=True):
    '''
    Bisection algorithm to solve an equation LHS = RHS('test') as a function of an independent
    variable defined below, 'test'.
    
    Lb        : left (lower) bound of the solution space
    Rb        : right (upper) bound of the solution space
    RHSfunc   : function returning the value of RHS given a 'test' value
    LHS       : constant representing the desired value of RHS
    inputs    : additional inputs required by RHSfunc. Empty by default
    direction : variable controlling the assumed derivative sign of RHSfunc
    iterLim   : maximum allowed iterations. 500 by default
    errAllow  : allowable error between RHS and LHS. 1e-6 by default
    verbose   : controls whether debug strings are printed to console
    
    returns   : critical solution value to LHS = RHS. Returns None if no solution is found within iterLim
    '''
    crit = None
    
    # For iterLim iterations...
    for i in range(iterLim):
        test = Lb + (Rb - Lb)/2 # Assume a test value at the midpoint of the solution space
        
        # Get the RHS value for a given test value. Add additional inputs if there are any
        if inputs == []: 
            RHS = RHSfunc(test)
        else:
            RHS = RHSfunc(test, *inputs)
        
        # print(f'{RHS} == {LHS}')

        # Check error of RHS relative to LHS. If within error tolerance, return the solution
        if abs(LHS - RHS)/LHS < errAllow:
            crit = test
            break
        
        # Assuming a solution was not found, adjust the solution space accordingly
        if direction:
            if RHS > LHS:
                Lb = test
            else:
                Rb = test
        else:
            if RHS < LHS:
                Lb = test
            else:
                Rb = test
                
    # Print debug strings if desired
    if verbose:
        if crit is None:
            print(f'Solver failed to converge in {iterLim} iterations')
            return None
        else:
            print(f'Solver converged at {round(crit,5)}; took {i} iterations')
    
    return crit
