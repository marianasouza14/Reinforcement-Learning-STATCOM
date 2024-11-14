import numpy as np

# Main Program
def main():
    # Read system data
    PowerFlowsData()
    YR, YI = YBus(tlsend, tlrec, tlresis, tlreac, tlsuscep, tlcond, shbus, shresis, shreac, ntl, nbb, nsh)
    VM, VA, it = NewtonRaphson(nmax, tol, itmax, ngn, nld, nbb, bustype, genbus, loadbus, PGEN, QGEN, QMAX, QMIN, PLOAD, QLOAD, YR, YI, VM, VA)
    PQsend, PQrec, PQloss, PQbus = PQflows(nbb, ngn, ntl, nld, genbus, loadbus, tlsend, tlrec, tlresis, tlreac, tlcond, tlsuscep, PLOAD, QLOAD, VM, VA)
    
    print(f'Number of iterations: {it}')
    print(f'Voltage magnitude at bus: {VM}')
    print(f'Voltage angle at bus (Degrees): {VA * 180 / np.pi}')
    print(f'Sending active and reactive powers: {PQsend}')
    print(f'Receiving active and reactive powers: {PQrec}')
    print(f'Sending active and reactive powers (WVAR): {PQsend * 100}')
    print(f'Receiving active and reactive powers (WVAR): {PQrec * 100}')

# YBus Function
def YBus(tlsend, tlrec, tlresis, tlreac, tlsuscep, tlcond, shbus, shresis, shreac, ntl, nbb, nsh):
    YR = np.zeros((nbb, nbb))
    YI = np.zeros((nbb, nbb))
    
    # Transmission lines contribution
    for kk in range(ntl):
        ii = tlsend[kk]
        jj = tlrec[kk]
        denom = tlresis[kk]**2 + tlreac[kk]**2
        YR[ii, ii] += tlresis[kk] / denom + 0.5 * tlcond[kk]
        YI[ii, ii] -= tlreac[kk] / denom + 0.5 * tlsuscep[kk]
        YR[ii, jj] -= tlresis[kk] / denom
        YI[ii, jj] += tlreac[kk] / denom
        YR[jj, ii] -= tlresis[kk] / denom
        YI[jj, ii] += tlreac[kk] / denom
        YR[jj, jj] += tlresis[kk] / denom + 0.5 * tlcond[kk]
        YI[jj, jj] -= tlreac[kk] / denom + 0.5 * tlsuscep[kk]

    # Shunt elements contribution
    for kk in range(nsh):
        ii = shbus[kk]
        denom = shresis[kk]**2 + shreac[kk]**2
        YR[ii, ii] += shresis[kk] / denom
        YI[ii, ii] -= shreac[kk] / denom
    
    return YR, YI

# Newton-Raphson Method
def NewtonRaphson(nmax, tol, itmax, ngn, nld, nbb, bustype, genbus, loadbus, PGEN, QGEN, QMAX, QMIN, PLOAD, QLOAD, YR, YI, VM, VA):
    D = np.zeros(nmax)
    flag = 0
    it = 1

    # Calculate net powers
    PNET, QNET = NetPowers(nbb, ngn, nld, genbus, loadbus, PGEN, QGEN, PLOAD, QLOAD)
    
    while it < itmax and flag == 0:
        # Calculated powers
        PCAL, QCAL = CalculatedPowers(nbb, VM, VA, YR, YI)
        
        # Check generator's reactive powers limits violations
        QNET, bustype = GeneratorsLimits(ngn, genbus, bustype, QGEN, QMAX, QMIN, QCAL, QNET, QLOAD, it, VM, nld, loadbus)
        
        # Power mismatches
        DPQ, DP, DQ, flag = PowerMismatches(nmax, nbb, tol, bustype, flag, PNET, QNET, PCAL, QCAL)
        
        # Jacobian formation
        JAC = NewtonRaphsonJacobian(nmax, nbb, bustype, PCAL, QCAL, VM, VA, YR, YI)
        
        # Solve for state variables
        D = np.linalg.solve(JAC, DPQ)
        
        # Update state variables
        VA, VM = StateVariablesUpdates(nbb, D, VA, VM)
        it += 1
    
    return VM, VA, it

# Calculate net scheduled powers
def NetPowers(nbb, ngn, nld, genbus, loadbus, PGEN, QGEN, PLOAD, QLOAD):
    PNET = np.zeros(nbb)
    QNET = np.zeros(nbb)
    
    for ii in range(ngn):
        PNET[genbus[ii]] += PGEN[ii]
        QNET[genbus[ii]] += QGEN[ii]
    
    for ii in range(nld):
        PNET[loadbus[ii]] -= PLOAD[ii]
        QNET[loadbus[ii]] -= QLOAD[ii]
    
    return PNET, QNET

# Calculate injected bus powers
def CalculatedPowers(nbb, VM, VA, YR, YI):
    PCAL = np.zeros(nbb)
    QCAL = np.zeros(nbb)

    for ii in range(nbb):
        PSUM = 0
        QSUM = 0
        for jj in range(nbb):
            PSUM += VM[ii] * VM[jj] * (YR[ii, jj] * np.cos(VA[ii] - VA[jj]) + YI[ii, jj] * np.sin(VA[ii] - VA[jj]))
            QSUM += VM[ii] * VM[jj] * (YR[ii, jj] * np.sin(VA[ii] - VA[jj]) - YI[ii, jj] * np.cos(VA[ii] - VA[jj]))
        PCAL[ii] = PSUM
        QCAL[ii] = QSUM

    return PCAL, QCAL

# Check generator's reactive powers limits violations
def GeneratorsLimits(ngn, genbus, bustype, QGEN, QMAX, QMIN, QCAL, QNET, QLOAD, it, VM, nld, loadbus):
    if it > 2:
        flag2 = 0
        for ii in range(ngn):
            jj = genbus[ii]
            if bustype[jj] == 2:
                if QCAL[jj] > QMAX[ii]:
                    QNET[genbus[ii]] = QMAX[ii]
                    bustype[jj] = 3
                    flag2 = 1
                elif QCAL[jj] < QMIN[ii]:
                    QNET[genbus[ii]] = QMIN[ii]
                    bustype[jj] = 3
                    flag2 = 1

                if flag2 == 1:
                    for ii in range(nld):
                        if loadbus[ii] == jj:
                            QNET[loadbus[ii]] -= QLOAD[ii]

    return QNET, bustype

# Compute power mismatches
def PowerMismatches(nmax, nbb, tol, bustype, flag, PNET, QNET, PCAL, QCAL):
    DPQ = np.zeros(nmax)
    DP = np.zeros(nbb)
    DQ = np.zeros(nbb)
    
    DP = PNET - PCAL
    DQ = QNET - QCAL

    for ii in range(nbb):
        if bustype[ii] == 1:
            DP[ii] = 0
            DQ[ii] = 0
        elif bustype[ii] == 2:
            DQ[ii] = 0
    
    kk = 0
    for ii in range(nbb):
        DPQ[kk] = DP[ii]
        DPQ[kk + 1] = DQ[ii]
        kk += 2

    for ii in range(nbb * 2):
        if np.abs(DPQ[ii]) < tol:
            flag = 1

    return DPQ, DP, DQ, flag

# Build the Jacobian matrix
def NewtonRaphsonJacobian(nmax, nbb, bustype, PCAL, QCAL, VM, VA, YR, YI):
    JAC = np.zeros((nmax, nmax))
    iii = 0
    
    for ii in range(nbb):
        jjj = 0
        for jj in range(nbb):
            if ii == jj:
                JAC[iii, jjj] = -QCAL[ii] - VM[ii]**2 * YI[ii, ii]
                JAC[iii, jjj + 1] = PCAL[ii] + VM[ii]**2 * YR[ii, ii]
