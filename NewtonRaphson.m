clc
clear all
%***- - - Programa Principal
PowerFlowsData; %Read system data
[YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,shbus,shresis,shreac,ntl,nbb,nsh);
[VM,VA,it] = NewtonRaphson(nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA);
[PQsend,PQrec,PQloss,PQbus] = PQflows(nbb,ngn,ntl,nld,genbus,loadbus,tlsend,tlrec,tlresis,tlreac,tlcond,tlsuscep,PLOAD,QLOAD,VM,VA);
it %Número de iterações
VM %Magnitude de tensão no nó (p.u.)
VA = VA*180/pi % Ângulo de tensão no nó (Deg)
PQsend %Sending active and reactive powers (p.u.)
PQrec %Receiving active and reactive powers (p.u.)
PQsend_WVAR = PQsend*100
PQrec_WVAR = PQrec*100
%End Programa Principal

% Matriz admitância
function [YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,shbus,shresis,shreac,ntl,nbb,nsh);
YR=zeros(nbb,nbb);
YI=zeros(nbb,nbb);
% Transmission lines contribution
for kk = 1:ntl
    ii = tlsend(kk);
    jj = tlrec(kk);
    denom = tlresis(kk)^2+tlreac(kk)^2;
    YR(ii,ii) = YR(ii,ii) + tlresis(kk)/denom + 0.5*tlcond(kk);
    YI(ii,ii) = YI(ii,ii) - tlreac(kk)/denom + 0.5*tlsuscep(kk);
    YR(ii,jj) = YR(ii,jj) - tlresis(kk)/denom;
    YI(ii,jj) = YI(ii,jj) + tlreac(kk)/denom;
    YR(jj,ii) = YR(jj,ii) - tlresis(kk)/denom;
    YI(jj,ii) = YI(jj,ii) + tlreac(kk)/denom;
    YR(jj,jj) = YR(jj,jj) + tlresis(kk)/denom + 0.5*tlcond(kk);
    YI(jj,jj) = YI(jj,jj) - tlreac(kk)/denom + 0.5*tlsuscep(kk);
end
% Shunt elements contribution
for kk = 1: nsh
    ii = shbus(kk);
    denom = shresis(kk)^2+shreac(kk)^2;
    YR(ii,ii) = YR(ii,ii) + shresis(kk)/denom;
    YI(ii,ii) = YI(ii,ii) - shreac(kk)/denom;
end
end
% End of function YBus

%Carry out iterative solution using the Newton-Raphson method
function [VM,VA,it] = NewtonRaphson(nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA);
% GENERAL SETTINGS
D = zeros(1,nmax);
flag = 0;
it = 1;
% CALCULATE NET POWERS
[PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD);
while ( it < itmax && flag==0 )
    % CALCULATED POWERS
    [PCAL,QCAL] = CalculatedPowers(nbb,VM,VA,YR,YI);
    % CHECK FOR POSSIBLE GENERATOR’S REACTIVE POWERS LIMITS VIOLATIONS
    [QNET,bustype] = GeneratorsLimits(ngn,genbus,bustype,QGEN,QMAX,QMIN,QCAL,QNET, QLOAD, it, VM, nld, loadbus);
    % POWER MISMATCHES
    [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL);
    % JACOBIAN FORMATION
    [JAC] = NewtonRaphsonJacobian(nmax,nbb,bustype,PCAL,QCAL,VM,VA,YR,YI);
    % SOLVE FOR THE STATE VARIABLES VECTOR
    D = JAC\DPQ';
    % UPDATE STATE VARIABLES
    [VA,VM] = StateVariablesUpdates(nbb,D,VA,VM);
    it = it + 1;
end
end
% End function Newton-Raphson

% Function to calculate the net scheduled powers
function [PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD);
% CALCULATE NET POWERS
PNET = zeros(1,nbb);
QNET = zeros(1,nbb);
for ii = 1: ngn
    PNET(genbus(ii)) = PNET(genbus(ii)) + PGEN(ii);
    QNET(genbus(ii)) = QNET(genbus(ii)) + QGEN(ii);
end
for ii = 1: nld
    PNET(loadbus(ii)) = PNET(loadbus(ii)) - PLOAD(ii);
    QNET(loadbus(ii)) = QNET(loadbus(ii)) - QLOAD(ii);
end
end
%End function NetPowers

%Function to calculate injected bus powers
function [PCAL,QCAL] = CalculatedPowers(nbb,VM,VA,YR,YI)
% Include all entries
PCAL = zeros(1,nbb);
QCAL = zeros(1,nbb);
for ii = 1: nbb
    PSUM = 0;
    QSUM = 0;
    for jj = 1: nbb
        PSUM = PSUM + VM(ii)*VM(jj)*(YR(ii,jj)*cos(VA(ii)-VA(jj)) + YI(ii,jj)*sin(VA(ii)-VA(jj)));
        QSUM = QSUM + VM(ii)*VM(jj)*(YR(ii,jj)*sin(VA(ii)-VA(jj)) - YI(ii,jj)*cos(VA(ii)-VA(jj)));
    end
    PCAL(ii) = PSUM;
    QCAL(ii) = QSUM;
end
end
%End of functionCalculatePowers

%Function to check whether or not solution is within generators limits
function [QNET,bustype] = GeneratorsLimits(ngn,genbus,bustype,QGEN,QMAX,QMIN,QCAL,QNET, QLOAD, it, VM, nld, loadbus);
% CHECK FOR POSSIBLE GENERATOR’S REACTIVE POWERS LIMITS VIOLATIONS
if it > 2
    flag2 = 0;
    for ii = 1: ngn
        jj = genbus(ii);
        if (bustype(jj) == 2)
            if ( QCAL(jj) > QMAX(ii) )
                QNET(genbus(ii)) = QMAX(ii);
                bustype(jj) = 3;
                flag2 = 1;
            elseif ( QCAL(jj) < QMIN(ii) )
                QNET(genbus(ii)) = QMIN(ii);
                bustype(jj) = 3;
                flag2 = 1;
            end
        if flag2 == 1
            for ii = 1:nld
                if loadbus(ii) == jj
                    QNET(loadbus(ii) == QNET(loadbus(ii)) - QLOAD(ii))
                end
            end
        end
        end
    end
end
end
%End function Generatorslimits

%Function to compute power mismatches
function [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL);
% POWER MISMATCHES
DPQ = zeros(1,nmax);
DP = zeros(1,nbb);
DQ = zeros(1,nbb);
DP = PNET - PCAL;
DQ = QNET - QCAL;
% To remove the active and reactive powers contributions of the slack bus and reactive power of all PV buses
for ii = 1: nbb
    if (bustype(ii) == 1 )
        DP(ii) = 0;
        DQ(ii) = 0;
    elseif (bustype(ii) == 2 )
        DQ(ii) = 0;
    end
end
% Re-arrange mismatch entries
kk = 1;
for ii = 1: nbb
    DPQ(kk) = DP(ii);
    DPQ(kk+1) = DQ(ii);
    kk = kk + 2;
end
% Check for convergence
for ii = 1: nbb*2
    if ( abs(DPQ) < tol)
        flag = 1;
    end
end
end
%End function PowerMismatches

%Function to built the Jacobian matrix
function [JAC] = NewtonRaphsonJacobian(nmax,nbb,bustype,PCAL,QCAL,VM,VA,YR,YI);
% JACOBIAN FORMATION
% Include all entries
JAC = zeros(nmax,nmax);
iii = 1;
for ii = 1: nbb
    jjj = 1;
    for jj = 1: nbb
        if ii == jj
            JAC(iii,jjj) = -QCAL(ii) - VM(ii)^2*YI(ii,ii);
            JAC(iii,jjj+1) = PCAL(ii) + VM(ii)^2*YR(ii,ii);
            JAC(iii+1,jjj) = PCAL(ii) - VM(ii)^2*YR(ii,ii);
            JAC(iii+1,jjj+1) = QCAL(ii) - VM(ii)^2*YI(ii,ii);
        else
            JAC(iii,jjj) = VM(ii)*VM(jj)*(YR(ii,jj)*sin(VA(ii)-VA(jj))-YI(ii,jj)*cos(VA(ii)-VA(jj)));
            JAC(iii+1,jjj) = -VM(ii)*VM(jj)*(YI(ii,jj)*sin(VA(ii)-VA(jj))+YR(ii,jj)*cos(VA(ii)-VA(jj)));
            JAC(iii,jjj+1) = -JAC(iii+1,jjj);
            JAC(iii+1,jjj+1) = JAC(iii,jjj);
        end
        jjj = jjj + 2;
    end
iii = iii + 2;
end
% Delete the voltage magnitude and phase angle equations of the slack bus and voltage magnitude equations corresponding to PV buses
for kk = 1: nbb
    if (bustype(kk) == 1)
        ii = kk*2-1;
        for jj = 1: 2*nbb
            if ii == jj
                JAC(ii,ii) = 1;
            else
                JAC(ii,jj) = 0;
                JAC(jj,ii) = 0;
            end
        end
    end
    if (bustype(kk) == 1) || (bustype(kk) == 2)
        ii = kk*2;
        for jj = 1: 2*nbb
            if ii == jj
                JAC(ii,ii) = 1;
            else
                JAC(ii,jj) = 0;
                JAC(jj,ii) = 0;
            end
        end
    end
end
end
%End of function NewtonRaphsonJacobian

%Function to update state variables
function [VA,VM] = StateVariablesUpdates(nbb,D,VA,VM)
iii = 1;
for ii = 1: nbb
    VA(ii) = VA(ii) + D(iii);
    VM(ii) = VM(ii) + D(iii+1)*VM(ii);
    iii = iii + 2;
end
end
%End function StateVariableUpdating

%Function to calculate the power flows
function [PQsend,PQrec,PQloss,PQbus] = PQflows(nbb,ngn,ntl,nld,genbus,loadbus,tlsend,tlrec,tlresis,tlreac,tlcond,tlsuscep,PLOAD,QLOAD,VM,VA);
PQsend = zeros(1,ntl);
PQrec = zeros(1,ntl);
% Calculate active and reactive powers at the sending and receiving
% ends of tranmsission lines
for ii = 1: ntl
    Vsend = ( VM(tlsend(ii))*cos(VA(tlsend(ii))) + VM(tlsend(ii))*sin(VA(tlsend(ii)))*i );
    Vrec = ( VM(tlrec(ii))*cos(VA(tlrec(ii))) + VM(tlrec(ii))*sin(VA(tlrec(ii)))*i );
    tlimped = tlresis(ii) + tlreac(ii)*i;
    current =(Vsend - Vrec) / tlimped + Vsend*( tlcond(ii) + tlsuscep(ii)*i )*0.5 ;
    PQsend(ii) = Vsend*conj(current);
    current =(Vrec - Vsend) / tlimped + Vrec*( tlcond(ii) + tlsuscep(ii)*i )*0.5 ;
    PQrec(ii) = Vrec*conj(current);
    PQloss(ii) = PQsend(ii) + PQrec(ii);
end
% Calculate active and reactive powers injections at buses
PQbus = zeros(1,nbb);
for ii = 1: ntl
    PQbus(tlsend(ii)) = PQbus(tlsend(ii)) + PQsend(ii);
    PQbus(tlrec(ii)) = PQbus(tlrec(ii)) + PQrec(ii);
end
% Make corrections at generator buses, where there is load, in order to get correct generators contributions
for ii = 1: nld
    jj = loadbus(ii);
    for kk = 1: ngn
        ll = genbus(kk);
        if jj == ll
            PQbus(jj) = PQbus(jj) + ( PLOAD(ii) + QLOAD(ii)*i );
        end
    end
end
end
%End function PQflows

