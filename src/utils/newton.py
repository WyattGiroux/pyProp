import numpy as np

def check_bounds(x, xlow, xhigh):
    for i in range(len(x)):
        if not xlow[i] <= x[i] <= xhigh[i]:
            return False
    return True


def newton_base(x0, f, J, itermax=500, reltol=1e-6, abstol=1e-6, verbose=False):
    iter = 1
    xold = x0
    xnew = np.zeros(len(x0))
    xs = [x0]
    dx = [[0,0]]
    ns = [np.inf]
    
    
    while iter < itermax:
        if np.shape(xold) == np.shape(np.array([0])):
            deltax = np.array([f(xold)[0] / J(xold)[0]])
        else:
            deltax = np.linalg.solve(J(xold), f(xold))
        xnew = xold - deltax
        wn = reltol * xnew + abstol
        
        norm = np.sum((deltax / wn)**2)
        xs.append(xnew)
        dx.append(deltax)
        ns.append(norm)
        
        if verbose:
            print(f'Iteration {iter}:\n\
> x: {np.round(xnew,3)}\n\
> dx: {np.round(deltax,3)}\n\
> norm: {np.round(norm,3)}\n')
        
        if norm < 1:
            iter_data = {
                'xs': np.array(xs),
                'dx': np.array(dx, dtype=object),
                'ns': np.array(ns)}
            if verbose: print(f'Solution Found:\n\
> Solution: {xnew}\n\
> Iterations: {iter}\n\
> Jacobians Required: {iter}\n')
            return xnew, deltax, norm, iter, iter_data
        
        iter += 1
        xold = xnew
    print('Reached Maximum Number of Iterations')
    return None, None, None, None, None
        
        
def newton_relax(x0, f, J, xLow, xHigh, itermax=500, recalcLim=5, trialMax=10, reltol=1e-6, abstol=1e-6, verbose=False, dRelax=np.sqrt(2)):
    # Create iteration counter and xk, xk+1
    iter = 1
    xold = x0
    xnew = np.zeros(len(x0))
    dxOld = np.zeros(len(x0))
    dxNew = np.zeros(len(x0))
    
    # Debug lists for recording iteration data
    xs = [x0]
    dx = [[0,0]]
    ns = [np.inf]
    rel = [1]
    nJ = 1
    
    # Relaxation factor starts at 1
    rFactor = 1
    
    # Initial Jacobian calculation
    Jcurr = J(xold)
    
    # Begin Newton iteration loop
    while iter < itermax:
        trials = 0
        rFactor = 1
        recalc = False # Flag for recalculating the jacobian
        
        # Get initial correction vector
        if np.shape(xold) == np.shape(np.array([0])):
            dxOld = np.array([f(xold)[0] / Jcurr[0]])
        else:
            dxOld = np.linalg.solve(Jcurr, f(xold))
        
        # Determine relaxation factor
        while trials < trialMax:      
            xnew = xold - rFactor * dxOld
            
            # Relax to keep solution in trust region determined by xLow, xHigh
            while not check_bounds(xnew, xLow, xHigh):
                trials += 1
                if trials > trialMax or recalc: break
                if trials % recalcLim == 0:
                    recalc = True
                    break
                            
                rFactor = rFactor / dRelax
                xnew = xold - rFactor * dxOld
            
            # Calculate trust-relaxed xk+1 and correction vector
            xnew = xold - rFactor * dxOld 
            if np.shape(xold) == np.shape(np.array([0])):
                dxNew = np.array([f(xnew)[0] / Jcurr[0]])
            else:
                dxNew = np.linalg.solve(Jcurr, f(xold))
                
            # Relax to ensure dxk+1 < dxk
            while np.linalg.norm(dxNew) > np.linalg.norm(dxOld):
                trials += 1
                if trials > trialMax or recalc: break
                if trials % recalcLim == 0:
                    recalc = True
                    break
            
                rFactor = rFactor / dRelax
                
                xnew = xold - rFactor * dxOld
                if np.shape(xold) == np.shape(np.array([0])):
                    dxNew = np.array([f(xnew)[0] / Jcurr[0]])
                else:
                    dxNew = np.linalg.solve(Jcurr, f(xold))
                
            # If the loop broke because the magnitude condition was met, the step is successful
            if np.linalg.norm(dxNew) <= np.linalg.norm(dxOld): 
                if verbose: print(f'Trials Required: {trials}')
                break                
            
            # Recalculate the jacobian if necessary
            if recalc:
                Jcurr = J(xold)
                nJ += 1
                if np.shape(xold) == np.shape(np.array([0])):
                    dxOld = np.array([f(xold)[0] / Jcurr[0]])
                else:
                    dxOld = np.linalg.solve(Jcurr, f(xold))
                rFactor = 1
                recalc = False
        
        if trials >= trialMax: print('Exceeded Number of Trials; Taking Step and Moving to Next Iteration')
        
        # Calculate weighted norm
        wn = reltol * xnew + abstol
        norm = np.sum((dxNew / wn)**2)
        
        # Update iteration info lists
        xs.append(xnew)
        dx.append(dxOld)
        ns.append(norm)
        rel.append(rFactor)
        
        if verbose:
            print(f'Iteration {iter}:\n\
> x: {np.round(xnew,3)}\n\
> dx: {np.round(dxOld,3)}\n\
> norm: {np.round(norm,3)}\n\
> relax: {np.round(rFactor, 3)}\n')
        
        # If the weighted norm is < 1, solver is converged
        if norm < 1:
            iter_data = {
                'xs': np.array(xs),
                'dx': np.array(dx, dtype=object),
                'ns': np.array(ns),
                'rel': np.array(rel)}
            if verbose: print(f'Solution Found:\n\
> Solution: {xnew}\n\
> Iterations: {iter}\n\
> Jacobians Required: {nJ}\n')
            
            return xnew, dxOld, norm, iter, iter_data
        
        # Increase iteration, update xk, and reset relaxation factor
        iter += 1
        xold = xnew
        rFactor = 1
    print('Reached Maximum Number of Iterations')
    return None, None, None, None, None