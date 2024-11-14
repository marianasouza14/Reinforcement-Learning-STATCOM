clc
clear all

%- - - Programa Principal ---
PowerFlowsData; % Ler informações do sistema de energia (LTs)
SSCData; % Ler informações do STATCOM
[YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,shbus,shresis,shreac,ntl,nbb,nsh); % Montar Matriz admitância
passo = 0.01;
P_base = 100; % Base é 100 MVA/MW/MVAr
[T_reforco, max_estados, linha_inicial, potencia_ref] = funcao_reforco(nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA,limStatcom,passo,SSCsend);

% Aprendizado por reforço
function [T_reforco, max_estados, linha_inicial, potencia_ref] = funcao_reforco (nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA,limStatcom,passo,SSCsend)
tolerancia = 0.0005;
valor_um = false;
max_estados = 2*limStatcom/passo + 1; % Usado como contador para o loop que determina os Q-estados - A soma de 1 é porque temos a "linha zero"
Qaux = QLOAD(2); % Variável auxiliar
limStatcomaux = limStatcom; % Variável auxiliar
Qc = zeros(1,max_estados);
Ql = zeros(1,max_estados);
for i=1:1:max_estados
    QLOAD(2) = Qaux + limStatcomaux;
    [VM,VA] = NewtonRaphson(nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA);
    % X(i) = round(VM(SSCsend(1)),3);
    X(i) = VM(SSCsend(1));
    if (i == 1)
        Qc(i)=-0.01; Ql(i)=-1;
    end
    if (i == max_estados)
        Qc(i)=-1; Ql(i)=-0.01;
    end
    if (i ~= 1) && (i ~= max_estados)
        if X(i)<1
            Qc(i)=-0.01; Ql(i)=-0.5;
        end
        if X(i)>1
            Qc(i)=-0.5; Ql(i)=-0.01;
        end
        % As condições abaixo analisam se eu tenho uma tensão igual a 1. Se não tiver, vai considerar uma aproximação de 1 com base na tolerância (tol) definida
        % anteriormente. A partir do momento queeu tenha alguma dessas condições, altero a variável auxiliar valor_um para evitar duas linhas com reforço igual a 1
        if X(i) == 1 && ~valor_um
            Qc(i)=1; Ql(i)=1;
            linha_inicial = i;
            valor_um = true;
            potencia_ref = limStatcomaux;
        elseif abs(X(i) - 1) < tolerancia && ~valor_um
            Qc(i)=1; Ql(i)=1;
            linha_inicial = i;
            valor_um = true;
            potencia_ref = limStatcomaux;
        end
    end
    limStatcomaux = limStatcomaux - passo;
end
% Criar uma matriz com variáveis X, Qc e Ql
QLOAD(2) = Qaux;
T_reforco = [X', Qc', Ql']; %'VariableNames', 'Estado (V)', 'Qc', 'Ql'; %Tabela de reforço
end

%Tabela Q-Valores
% estou em um estado - ou eu injeto reativo ou capacitivo
T_reforco_aux = T_reforco(:, [2,3]); % Tabela sem considerar a primeira coluna - Facilitar processos
Q_atual = T_reforco_aux;
% tol: tolerância para o critério de parada
% gamma: fator de desconto
Q_novo = T_reforco_aux;

% Variável para verificar se houve mudança significativa
convergiu = false;

% Critério de parada baseado na diferença entre iterações
iteracao = 0;
tol = 0.001;

% Montar tabela de Q-estados
while ~convergiu && iteracao < 1000
    % Marca que a tabela Q não foi atualizada
    delta = 0;
    % Primeiro, atualiza a tabela "andando" para cima a partir do estado cuja tensão é igual a 1
    for estado = (linha_inicial-1):-1:1
        for acao = 1:1:2
            acao_oposta = 3 - acao; % Ação oposta no estado atual
            % Ver o maior Q-valor no estado anterior (abaixo)
            valor_anterior1 = Q_atual(estado+1, acao);
            valor_anterior2 = Q_atual(estado+1, acao_oposta);
            valor_anterior = max(valor_anterior1, valor_anterior2);
            % Se não for o primeiro estado, pegará valor "próximo" (que seria o de cima) - Como não há linha acima do primeiro estado, tem essa condição
            if estado ~= 1
                valor_proximo1 = Q_atual(estado-1, acao);
                valor_proximo2 = Q_atual(estado-1, acao_oposta);
                valor_proximo = max(valor_proximo1, valor_proximo2);
            end
            % Se for o primeiro estado, analisará qual dos Q-estados é igual a (-1)
            if estado == 1
                if Q_atual(estado,acao) == -1
                    continue
                else
                    Q_novo(estado, acao) = T_reforco_aux(estado, acao) + max(valor_anterior);
                end
            else
            % Reforço recebido + desconto máximo entre os 3 valores comparados
            Q_novo(estado, acao) = T_reforco_aux(estado, acao) + max(valor_anterior, valor_proximo);
            end
            % Verifica a diferença absoluta entre o novo valor de Q e o antigo
            delta = max(delta, abs(Q_novo(estado, acao) - Q_atual(estado, acao)));
        end
    end
    for estado = (linha_inicial+1):1:max_estados
        for acao = 1:1:2
            acao_oposta = 3 - acao;
            % Ação do estado anterior (mesma ação)
            valor_anterior1 = Q_atual(estado-1, acao);
            valor_anterior2 = Q_atual(estado-1, acao_oposta);
            valor_anterior = max(valor_anterior1, valor_anterior2);
            % Se não for o último estado, pegará valor "próximo" (que seria o de baixo) - Como não há linha abaixo do último estado, tem essa condição
            if estado ~= max_estados
                valor_proximo1 = Q_atual(estado + 1, acao);
                valor_proximo2 = Q_atual(estado + 1, acao_oposta);
                valor_proximo = max(valor_proximo1, valor_proximo2);
            end
            % Se for o último estado, analisará qual dos Q-estados é igual a (-1)
            if estado == max_estados
                if Q_atual(estado,acao) == -1
                    continue
                else
                    Q_novo(estado, acao) = T_reforco_aux(estado, acao) + max(valor_anterior);
                end
            else
            % Reforço recebido + desconto máximo entre os 3 valores comparados
            Q_novo(estado, acao) = T_reforco_aux(estado, acao) + max(valor_anterior, valor_proximo);
            end
            if estado == 1
                continue
            else
            % Reforço recebido + desconto máximo entre os 3 valores comparados
            Q_novo(estado, acao) = T_reforco_aux(estado, acao) + max(valor_anterior, valor_proximo);
            end
            % Verifica a diferença absoluta entre o novo valor de Q e o antigo
            delta = max(delta, abs(Q_novo(estado, acao) - Q_atual(estado, acao)));
        end
    end
    % Atualiza o critério de parada com a diferença máxima observada
    if delta < tol
        convergiu = true;
    end
    Q_atual = round(Q_novo,4);
    % Incrementa o contador de iteração
    iteracao = iteracao + 1;
end
% Exibe o número de iterações necessárias para convergir
fprintf('A tabela Q convergiu após %d iterações.\n', iteracao);
Q_valor = [T_reforco(:,1), Q_atual];

% Critério de parada - encontrar a linha onde a recompensa é 1
% Variáveis de entrada
referencia = 0.99; % Exemplo de valor de referência
soma_passos = 0;  % Acumular o valor somado dos passos

% Encontrar a linha correspondente ao valor da referência
[~, linha_referencia] = min(abs(Q_valor(:, 1) - referencia));
% Verificar se a referência foi encontrada na tabela
if isempty(linha_referencia)
    error('Valor de referência não encontrado na tabela!');
end
% Ponteiro movendo de acordo com a referência
if referencia > 1
    % Caso a referência seja maior que 1, desce a tabela
    for i = linha_referencia:-1:1
        if Q_valor(i, 2) == 1 || Q_valor(i, 3) == 1
            % Critério de parada (reforço = 1)
            break;
        end
        soma_passos = soma_passos + 0.01;  % Soma -0.01 a cada passo
    end
elseif referencia < 1
    % Caso a referência seja menor que 1, sobe a tabela
    for i = linha_referencia:1:size(Q_valor, 1)
        if Q_valor(i, 2) == 1 || Q_valor(i, 3) == 1
            % Critério de parada (reforço = 1)
            break;
        end
        soma_passos = soma_passos + 0.01;  % Soma +0.01 a cada passo
    end
else
    disp('A referência já é 1, sem necessidade de andar.');
end


% Exibir resultado da soma dos passos
% fprintf('A soma dos passos é: %.2f\n', soma_passos);
% if soma_passos < 0 
%     fprintf('Deve-se injetar %.2f MVAr capacitivo (ou absorver %.2f MVAr indutivo\n', P_base*soma_passos);
% else
%    fprintf('Deve-se injetar %.2f MVAr indutivo\n', P_base*soma_passos);
% end


% [VM,VA,it] = NewtonRaphson(nmax,tol,itmax,ngn,nld,nbb,bustype,genbus,loadbus,PGEN,QGEN,QMAX,QMIN,PLOAD,QLOAD,YR,YI,VM,VA);
% [PQsend,PQrec,PQloss,PQbus] = PQflows(nbb,ngn,ntl,nld,genbus,loadbus,tlsend,tlrec,tlresis,tlreac,tlcond,tlsuscep,PLOAD,QLOAD,VM,VA);
% % it %Número de iterações
% VM %Magnitude de tensão no nó (p.u.)
% VA = VA*180/pi % Ângulo de tensão no nó (Deg)
% PQsend %Sending active and reactive powers (p.u.)
% PQrec %Receiving active and reactive powers (p.u.)
% PQsend_WVAR = PQsend*100
% PQrec_WVAR = PQrec*100

%End Programa Principal

% Matriz admitância
function [YR,YI] = YBus(tlsend,tlrec,tlresis,tlreac,tlsuscep,tlcond,shbus,shresis,shreac,ntl,nbb,nsh);
YR=zeros(nbb,nbb);
YI=zeros(nbb,nbb);
% Transmission lines contribution
for kk = 1:ntl
%kk - laço para cada linha de transmissão
    ii = tlsend(kk);
    jj = tlrec(kk);
%tlsend - de onde está saindo o fluxo - equivalente ao k do material estudado
%tlrec - onde está entrando o fluxo - equivalente ao m do material estudado
    denom = tlresis(kk)^2+tlreac(kk)^2;
    YR(ii,ii) = YR(ii,ii) + tlresis(kk)/denom + 0.5*tlcond(kk);
    YI(ii,ii) = YI(ii,ii) - tlreac(kk)/denom + 0.5*tlsuscep(kk);
% Diagonal principal - soma os elementos shunt nas admitâncias conectadas na barra - elementos conectados à barra ii
    YR(ii,jj) = YR(ii,jj) - tlresis(kk)/denom;
    YI(ii,jj) = YI(ii,jj) + tlreac(kk)/denom;
% Outros elementos - não leva em consideração os elementos shunt, e são o negativo das admitâncias entre as barras
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
flag = 0; %é utilizada para parar o algoritmo quando esse valor passa a ser 1 na comparação em outra função
it = 1;
% CALCULATE NET POWERS
[PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD);
while ( it < itmax && flag==0 )
    % CALCULATED POWERS (Pesp de cada barra)
    [PCAL,QCAL] = CalculatedPowers(nbb,VM,VA,YR,YI);
    % CHECK FOR POSSIBLE GENERATOR’S REACTIVE POWERS LIMITS VIOLATIONS
    [QNET,bustype] = GeneratorsLimits(ngn,genbus,bustype,QGEN,QMAX,QMIN,QCAL,QNET, QLOAD, it, VM, nld, loadbus);
    % POWER MISMATCHES (calular as diferenças - Deltas)
    [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL);
    % JACOBIAN FORMATION
    [JAC] = NewtonRaphsonJacobian(nmax,nbb,bustype,PCAL,QCAL,VM,VA,YR,YI);
    % SOLVE FOR THE STATE VARIABLES VECTOR
    D = JAC\DPQ';
    % UPDATE STATE VARIABLES
    [VA,VM] = StateVariablesUpdates(nbb,D,VA,VM);
    if flag == 1
    break
    end
    it=it+1;
end
end
% End function Newton-Raphson

% Function to calculate the net scheduled powers
function [PNET,QNET] = NetPowers(nbb,ngn,nld,genbus,loadbus,PGEN,QGEN,PLOAD,QLOAD);
% CALCULATE NET POWERS
PNET = zeros(1,nbb);
QNET = zeros(1,nbb);
for ii = 1:ngn
    PNET(genbus(ii)) = PNET(genbus(ii)) + PGEN(ii);
    QNET(genbus(ii)) = QNET(genbus(ii)) + QGEN(ii);
end
for ii = 1:nld
    PNET(loadbus(ii)) = PNET(loadbus(ii)) - PLOAD(ii);
    QNET(loadbus(ii)) = QNET(loadbus(ii)) - QLOAD(ii);
end
end
%Observa-se aqui que o que foi feito é a fórmula estudada: Pesp=Pger-Pcons;
%End function NetPowers

%Function to calculate injected bus powers
function [PCAL,QCAL] = CalculatedPowers(nbb,VM,VA,YR,YI)
% Include all entries
PCAL = zeros(1,nbb);
QCAL = zeros(1,nbb);
%Nessa função, será feito o cálculo do "somatório" utilizado na diferença de potências
for ii = 1: nbb
%ii é equivalente ao k na teoria estudada
    PSUM = 0;
    QSUM = 0;
    for jj = 1: nbb
%jj é equivalente ao m na teoria estudada
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
function [DPQ,DP,DQ,flag] = PowerMismatches(nmax,nbb,tol,bustype,flag,PNET,QNET,PCAL,QCAL)
% POWER MISMATCHES
DPQ = zeros(1,nmax);
DP = zeros(1,nbb);
DQ = zeros(1,nbb);
DP = PNET - PCAL;
DQ = QNET - QCAL;
% PNET=Pesp e QNET=Qesp
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
function [JAC] = NewtonRaphsonJacobian(nmax,nbb,bustype,PCAL,QCAL,VM,VA,YR,YI)
% JACOBIAN FORMATION
% Include all entries
JAC = zeros(nmax,nmax);
iii = 1;
for ii = 1: nbb
    jjj = 1;
    for jj = 1: nbb
        if ii == jj %Diagonal principal
            JAC(iii,jjj) = -QCAL(ii) - VM(ii)^2*YI(ii,ii); %Derivada de P em relação a theta
            JAC(iii,jjj+1) = PCAL(ii) + VM(ii)^2*YR(ii,ii); %Derivada de P em relação a V
            JAC(iii+1,jjj) = PCAL(ii) - VM(ii)^2*YR(ii,ii); %Derivada de Q em relação a theta
            JAC(iii+1,jjj+1) = QCAL(ii) - VM(ii)^2*YI(ii,ii); %Derivada de Q em reação a V
%Fórmulas contidas no livro, e destrinchadas nas anotações e ChatGPT
        else
            JAC(iii,jjj) = VM(ii)*VM(jj)*(YR(ii,jj)*sin(VA(ii)-VA(jj))-YI(ii,jj)*cos(VA(ii)-VA(jj))); %Derivada de P em relação a theta
            JAC(iii+1,jjj) = -VM(ii)*VM(jj)*(YI(ii,jj)*sin(VA(ii)-VA(jj))+YR(ii,jj)*cos(VA(ii)-VA(jj))); %Derivada de Q em relação a theta
            JAC(iii,jjj+1) = -JAC(iii+1,jjj); %Derivada de P em relação a V
            JAC(iii+1,jjj+1) = JAC(iii,jjj); %Derivada de Q em reação a V
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
% Calculate active and reactive powers at the sending and receiving ends of tranmsission lines
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

