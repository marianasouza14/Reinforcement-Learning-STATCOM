%The following convention is used for the four types of buses available in conventional power flow studies:
%bustype = 1 is slack or swing bus
%bustype = 2 is generator PV bus
%bustype = 3 is load PQ bus
%bustype = 4 is generator PQ bus

%The five buses in the network shown in Figure 4.6 are numbered for the purpose of the power flow solution, as follows:
%North = 1
%South = 2
%Lake = 3
%Main = 4
%Elm = 5

%Bus data
%nbb = number of buses
%bustype = type of bus
%VM = nodal voltage magnitude
%VA = nodal voltage phase angle
nbb = 5 ;
bustype(1) = 1 ; VM(1) = 1.06 ; VA(1) =0 ;
bustype(2) = 2 ; VM(2) = 1 ; VA(2) =0 ;
bustype(3) = 3 ; VM(3) = 1 ; VA(3) =0 ;
bustype(4) = 3 ; VM(4) = 1 ; VA(4) =0 ;
bustype(5) = 3 ; VM(5) = 1 ; VA(5) =0 ;

%Generator data
%ngn = number of generators
%genbus = generator bus number
%PGEN = scheduled active power contributed by the generator
%QGEN = scheduled reactive power contributed by the generator
%QMAX = generator reactive power upper limit
%QMIN = generator reactive power lower limit
ngn = 2 ;
genbus(1) = 1 ; PGEN(1) = 0 ; QGEN(1) = 0 ; QMAX(1) = 5 ; QMIN(1) = -5 ;
genbus(2) = 2 ; PGEN(2) = 0.4 ; QGEN(2) = 0 ; QMAX(2) = 3 ; QMIN(2) = -3 ;

%Transmission line data
%ntl = number of transmission lines
%tlsend = sending end of transmission line
%tlrec = receiving end of transmission line
%tlresis = series resistance of transmission line
%tlreac = series reactance of transmission line
%tlcond = shunt conductance of transmission line
%tlsuscep = shunt susceptance of transmission line
ntl = 7 ;
tlsend(1) = 1 ; tlrec(1) = 2 ; tlresis(1) = 0.02 ; tlreac(1) = 0.06 ;
tlcond(1) = 0 ; tlsuscep(1) = 0.06 ;
tlsend(2) = 1 ; tlrec(2) = 3 ; tlresis(2) = 0.08 ; tlreac(2) = 0.24 ;
tlcond(2) = 0 ; tlsuscep(2) = 0.05 ;
tlsend(3) = 2 ; tlrec(3) = 3 ; tlresis(3) = 0.06 ; tlreac(3) = 0.18 ;
tlcond(3) = 0 ; tlsuscep(3) = 0.04 ;
tlsend(4) = 2 ; tlrec(4) = 4 ; tlresis(4) = 0.06 ; tlreac(4) = 0.18 ;
tlcond(4) = 0 ; tlsuscep(4) = 0.04 ;
tlsend(5) = 2 ; tlrec(5) = 5 ; tlresis(5) = 0.04 ; tlreac(5) = 0.12 ;
tlcond(5) = 0 ; tlsuscep(5) = 0.03 ;
tlsend(6) = 3 ; tlrec(6) = 4 ; tlresis(6) = 0.01 ; tlreac(6) = 0.03 ;
tlcond(6) = 0 ; tlsuscep(6) = 0.02 ;
tlsend(7) = 4 ; tlrec(7) = 5 ; tlresis(7) = 0.08 ; tlreac(7) = 0.24 ;
tlcond(7) = 0 ; tlsuscep(7) = 0.05 ;

%Shunt data
%nsh = number of shunt elements
%shbus = shunt element bus number
%shresis = resistance of shunt element
%shreac = reactance of shunt element:
%+ve for inductive reactance and –ve for capacitive reactance
nsh = 0 ;
shbus(1) = 0 ; shresis(1) = 0 ; shreac(1) = 0 ;

%Load data
%nld = number of load elements
%loadbus = load element bus number
%PLOAD = scheduled active power consumed at the bus
%QLOAD = scheduled reactive power consumed at the bus
nld = 4 ;
loadbus(1) = 2 ; PLOAD(1) = 0.2 ; QLOAD(1) = 0.1 ;
loadbus(2) = 3 ; PLOAD(2) = 0.45 ; QLOAD(2) = 0.15 ;
loadbus(3) = 4 ; PLOAD(3) = 0.4 ; QLOAD(3) = 0.05 ;
loadbus(4) = 5 ; PLOAD(4) = 0.6 ; QLOAD(4) = 0.1 ;
%General parameters
%itmax = maximum number of iterations permitted before the iterative process is terminated – protection against infinite iterative loops
%tol = criterion tolerance to be met before the iterative solution is successfully brought to an end
itmax = 100;
tol = 1e-12;
nmax = 2*nbb;
%End of function PowerFlowsData