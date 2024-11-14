clear all
clc
% - - - Main STATCOM Program
PowerFlowsData; %Function to read network data
SSCData; %Function to read the STATCOM data
[YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,ntl,nbb);
[VM,VA,it,Vvr,Tvr] = SSCNewtonRaphson(tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA,NSSC,SSCsend,Xvr,TarVol,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr,VvrHi,VvrLo);
[PQsend,PQrec,PQloss,PQbus] = PQflows(nbb,ngn,ntl,nld,genbus,loadbus,tlsend,tlrec,tlresis,tlreac,tlcond,tlsuscep,PLOAD,QLOAD,VM,VA);
[Psend,Qsend,PSSC,QSSC] = SSCPQPowers(VM,VA,NSSC,SSCsend,Xvr,Vvr,Tvr);
%Print results
it %Number of iterations
VM %Nodal voltage magnitude (p.u)
VA=VA*180/pi %Nodal voltage phase angles (Deg)
Vvr %Final voltage magnitude source (p.u.)
Tvr=Tvr*180/pi %Final voltage phase angle source (Deg)
PQsend=Psend + j*Qsend; %Active and reactive powers in bus (p.u.)
PQSSC=PSSC + j*QSSC; %Active and reactive powers in STATCOM (p.u.)
% End of MAIN STATCOM PROGRAM

%Carry out iterative solution using the Newton–Raphson method
function [VM,VA,it,Vvr,Tvr] = SSCNewtonRaphson(tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA,NSSC,SSCsend,Xvr,TarVol,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr,VvrHi,VvrLo)
% GENERAL SETTINGS
flag = 0;
it = 1;
% CALCULATE NET POWERS
[PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD);
while ( it < itmax && flag==0 )
    % CALCULATED POWERS
    [PCAL,QCAL] = CalculatedPowers(nbb,VM,VA,YR,YI);
    %STATCOM CALCULATED POWERS
    [PCAL,QCAL,PSSC,QSSC] = SSCCalculatePowers(PCAL,QCAL,VM,VA,NSSC,SSCsend,Xvr,Vvr,Tvr);
    %POWER MISMATCHES
    [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL);
    %STATCOM POWER MISMATCHES
    [DPQ] = SSCMismatches(DPQ,nbb,VM,VA,NSSC,SSCsend,Xvr,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr);
    % JACOBIAN FORMATION
    [JAC] = NewtonRaphsonJacobian(nbb,bustype,PCAL,QCAL,DPQ,VM,VA,YR,YI);
    % STATCOM JACOBIAN
    [JAC] = SSCJacobian(nbb,JAC,VM,VA,NSSC,SSCsend,Xvr,TarVol,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr);
    % SOLVE FOR THE STATE VARIABLES VECTOR
    D = JAC\DPQ';
    % UPDATE THE STATE VARIABLES VALUES, WITH TRUNCATED CORRECTIONS IF NECESSARY (VM increments < +-0.1 p.u. and VA inrements < +- 5 deg)
    [VA,VM] = StateVariablesUpdating(nbb,D,VA,VM,it);
    %UPDATE STATCOM STATE VARIABLES
    [VM,Vvr,Tvr] = SSCUpdating(nbb,D,VM,VA,NSSC,SSCsend,TarVol,VSta,Psp,Vvr,Tvr);
    %CHECK VOLTAGE MAGNITUDE LIMITS
    [Vvr] = SSCLimits(NSSC,Vvr,VvrHi,VvrLo);
    it = it + 1;
end
end


%Function to calculate injected bus powers by the STATCOM
function [PCAL,QCAL,PSSC,QSSC] = SSCCalculatePowers(PCAL,QCAL,VM,VA,NSSC,SSCsend,Xvr,Vvr,Tvr);
for ii = 1 : NSSC
B(ii)=1/Xvr(ii);
A1 = Tvr(ii)-VA(SSCsend(ii));
A2 = VA(SSCsend(ii))-Tvr(ii);
PCAL(SSCsend(ii)) = PCAL(SSCsend(ii)) + VM(SSCsend(ii))*Vvr(ii)*(B(ii)*sin(A2));
QCAL(SSCsend(ii)) = QCAL(SSCsend(ii)) + VM(SSCsend(ii))^2*B(ii) - Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A2));
PSSC(ii) = Vvr(ii)*VM(SSCsend(ii))*(B(ii)*sin(A1));
QSSC(ii) = - Vvr(ii)^2*B(ii) + Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A1));
end
end

%Function to compute power mismatches for the STATCOM
function [DPQ] = SSCMismatches(DPQ,nbb,VM,VA,NSSC,SSCsend,Xvr,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr);
for ii = 1 : NSSC
    B(ii)=1/Xvr(ii);
    A1 = Tvr(ii)-VA(SSCsend(ii));
    A2 = VA(SSCsend(ii))-Tvr(ii);
    Pcal = VM(SSCsend(ii))*Vvr(ii)*(B(ii)*sin(A2));
    Qcal = - VM(SSCsend(ii))^2*B(ii) + Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A2));
    DPQ(2*(nbb + ii)-1) = Pcal - Psp(ii);
    if (QSta(ii) == 1)
    DPQ(2*(nbb + ii)) = Qcal - Qsp(ii);
    else
    DPQ(2*(nbb + ii)) = 0;
    end
end
end

%Function to add the STATCOM elements to Jacobian matrix
function [JAC] = SSCJacobian(nbb,JAC,VM,VA,NSSC,SSCsend,Xvr,TarVol,VSta,Psp,PSta,Qsp,QSta,Vvr,Tvr);
for ii = 1 : NSSC
    B(ii)=1/Xvr(ii);
    if VSta(ii) == 1
    JAC(: , 2*SSCsend(ii) )=0;
    end
    JAC(2*(nbb + ii)-1,2*(nbb + ii)-1) = 1;
    JAC(2*(nbb + ii),2*(nbb + ii)) = 1;
    A1 = Tvr(ii)-VA(SSCsend(ii));
    A2 = VA(SSCsend(ii))-Tvr(ii);
    Pcal = - VM(SSCsend(ii))*Vvr(ii)*( + B(ii)*sin(A2));
    DQcal = Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A2));
    Pssc = - Vvr(ii)*VM(SSCsend(ii))*(B(ii)*sin(A1));
    DQssc = Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A1));
    JAC(2*SSCsend(ii)-1,2*SSCsend(ii)-1) = JAC(2*SSCsend(ii)-1,2*SSCsend(ii)-1) + VM(SSCsend(ii))^2*B(ii);
    JAC(2*SSCsend(ii),2*SSCsend(ii)-1) = JAC(2*SSCsend(ii),2*SSCsend(ii)-1) - Pcal;
    if (QSta(ii) == 1 )
    JAC(2*SSCsend(ii)-1,2*SSCsend(ii)) = JAC(2*SSCsend(ii)-1,2*SSCsend(ii)) - Pcal;
    JAC(2*SSCsend(ii),2*SSCsend(ii)) = JAC(2*SSCsend(ii),2*SSCsend(ii)) + VM(SSCsend(ii))^2*B(ii);
    else
    JAC(2*SSCsend(ii)-1,2*SSCsend(ii)) = JAC(2*SSCsend(ii)-1,2*SSCsend(ii)) - Pssc;
    JAC(2*SSCsend(ii),2*SSCsend(ii)) = JAC(2*SSCsend(ii),2*SSCsend(ii)) - DQssc;
    end
if (PSta(ii) == 1)
JAC(2*(nbb + ii)-1,2*SSCsend(ii)-1) = JAC(2*(nbb + ii)-1, 2*SSCsend(ii)-1) + DQcal;
JAC(2*SSCsend(ii)-1,2*(nbb + ii)-1) = JAC(2*SSCsend(ii)-1,2*(nbb + ii)-1) - DQssc;
JAC(2*SSCsend(ii),2*(nbb + ii)-1) = JAC(2*SSCsend(ii),2*(nbb + ii)-1) - Pssc;
JAC(2*(nbb + ii)-1,2*(nbb + ii)-1) = - DQssc;
if (QSta == 1)
JAC(2*(nbb + ii),2*(nbb + ii)-1) = JAC(2*(nbb + ii),2*(nbb + ii)-1) - Pssc;
JAC(2*(nbb + ii)-1,2*SSCsend(ii)) = JAC(2*(nbb + ii)-1,2*SSCsend(ii)) - Pcal;
else
JAC(2*(nbb + ii),2*(nbb + ii)-1) = 0.0;
JAC(2*(nbb + ii)-1,2*SSCsend(ii)) = JAC(2*(nbb + ii)-1,2*SSCsend(ii)) + Pssc;
end
else
JAC(2*(nbb + ii)-1,2*(nbb + ii)-1) = 1.0;
end
if (QSta(ii) == 1)
JAC(2*(nbb + ii),2*SSCsend(ii)-1) = JAC(2*(nbb + ii),2*SSCsend(ii)-1)- Pcal;
JAC(2*(nbb + ii),2*SSCsend(ii)) = JAC(2*(nbb + ii),2*SSCsend(ii)) + DQcal;
JAC(2*SSCsend(ii)-1,2*(nbb + ii)) = JAC(2*SSCsend(ii)-1,2*(nbb + ii)) + Pssc;
JAC(2*SSCsend(ii),2*(nbb + ii)) = JAC((nbb + ii),2*(nbb + ii)) - DQcal;
JAC(2*(nbb + ii),2*(nbb + ii)) = -2*Vvr(ii)^2*B(ii) + DQssc;
if (PSta(ii) == 1)
    JAC(2*(nbb + ii)-1,2*(nbb + ii)) = JAC(2*(nbb + ii)-1,2*(nbb + ii)) - Pssc;
else
JAC(2*(nbb + ii)-1,2*(nbb + ii)) = 0.0;
end
else
JAC(2*(nbb + ii),2*(nbb + ii)) = 1.0;
end
end
end

%Function to update STATCOM state variable
function [VM,Vvr,Tvr] = SSCUpdating(nbb,D,VM,VA,NSSC,SSCsend,TarVol,VSta,Psp,Vvr,Tvr);
for ii = 1 : NSSC
if (VSta(ii) == 1)
% Adjust the Volatge Magnitud target
Vvr(ii) = Vvr(ii) + Vvr(ii)*D(2*SSCsend(ii));
VM(SSCsend(ii)) = TarVol(ii);
if (Psp(ii) == 0)
Tvr(ii) = VA(SSCsend(ii));
else
Tvr(ii) = Tvr(ii) + D(2*(nbb + ii)-1);
end
else
Vvr(ii) = Vvr(ii) + Vvr(ii)*D(2*(nbb + ii));
Tvr(ii) = VA(SSCsend(ii));
end
end
end

%Function to check source voltages limits in the STATCOM
function [Vvr] = SSCLimits(NSSC,Vvr,VvrHi,VvrLo);
for ii = 1 : NSSC
%Check STATCOM Vvr Limits
if (Vvr(ii) > VvrHi(ii))
Vvr(ii) = VvrHi(ii);
elseif (Vvr(ii) < VvrLo(ii))
Vvr(ii) = VvrLo(ii);
end
end
end

%Function to calculate the power flows in the STATCOM
function [Psend,Qsend,PSSC,QSSC] = SSCPQPowers(VM,VA,NSSC,SSCsend,Xvr,Vvr,Tvr);
for ii = 1 : NSSC
B(ii)=1/Xvr;
A1 = Tvr(ii)-VA(SSCsend(ii));
A2 = VA(SSCsend(ii))-Tvr(ii);
Psend(ii) = VM(SSCsend(ii))*Vvr(ii)*(B(ii)*sin(A2));
Qsend(ii) = - VM(SSCsend(ii))^2*B(ii) + Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A2));
PSSC(ii) = Vvr(ii)*VM(SSCsend(ii))*(B(ii)*sin(A1));
QSSC(ii) = - Vvr(ii)^2*B(ii) + Vvr(ii)*VM(SSCsend(ii))*(B(ii)*cos(A1));
end
end

% Matriz admitância
function [YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,ntl,nbb)
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
end
% End of function YBus


%Function to calculate the power flows
function [PQsend,PQrec,PQloss,PQbus] = PQflows(nbb,ngn,ntl,nld,genbus,loadbus,tlsend,tlrec,tlresis,tlreac,tlcond,tlsuscep,PLOAD,QLOAD,VM,VA)
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

%Function to compute power mismatches
function [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL)
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

%Function to built the Jacobian matrix
function [JAC] = NewtonRaphsonJacobian(nmax,nbb,bustype,PCAL,QCAL,VM,VA,YR,YI)
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


%This function is used exclusively to enter data for:
% STATIC SYNCHRONOUS COMPENSATOR (STATCOM)
% NSSC : Number of STATCOM’s
% SSCsend: STATCOM’s bus
% Xvr : Converter’s reactance (p.u.)
% TarVol: Target nodal voltage magnitude (p.u.)
% VSta : Indicate the control status over nodal voltage magnitude: 1 is
% on; 0 is off
% Psp : Target active power flow (p.u.)
% PSta : Indicate the control status over active power: 1 is on; 0 is off
% Qsp : Target reactive power flow (p.u.)
% QSta : Indicate the control status over reactive power:1 is on; 0 is off
% Vvr : Initial condition for the source voltage magnitude (p.u.)
% Tvr : Initial condition for the source voltage angle (deg)
% VvrHi : Lower limit source voltage magnitude (p.u.)
% VvrLo : higher limit source voltage magnitude (p.u.)
NSSC = 1;
SSCsend(1)=3; Xvr(1)=10; TarVol(1)=1.0; VSta(1)=1;
Psp(1)=0.0; PSta(1)=1; Qsp(1)=0.0; QSta(1)=0;
Vvr(1)=1.0; Tvr(1)=0.0; VvrHi(1)=1.1; VvrLo(1)=0.9;


% Function to calculate the net scheduled powers
function [PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD)
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