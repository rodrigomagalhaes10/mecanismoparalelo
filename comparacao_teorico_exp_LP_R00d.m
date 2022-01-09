% Universidade de São Paulo
% Programa de Pós Graduação em Engenharia Mecânica
% Escola Politécnica
% Nome: Rodrigo Magalhães Cruz Alves Silva
% NUSP: 9853634


% Cinemática direta e inversa de um manipulador assimétrico de cinemática paralela
  
% deleta todas as variáveis na memória
clear all

% fecha todas as janelas abertas do sistema
close all

% limpa a tela
clc

% Fator de diferença entre as molas
fator=1.0;

% dados da geometria
L0=260e-3;
L1=345e-3;
L2=440e-3;
L3=L1;
L4=L2;

% L5=0.00010;
L5=160e-3/2;

Lc=200e-3;

% Dados cadeia passiva
L6=100e-3;

L7=385e-3;
L8=L7;
A6=pi/4*(14.3e-3)^2;
A7=pi/4*(14.3e-3)^2;
A8=A7;
I6=(30.7e-3)^3*(6.1e-3)/12;

% Parâmetros de inércia Mk mi Izi para os 8 elos
m1=275.166e-3;
Iz1=4831958.5e-9+(15+142)*(343/8)^2*10^(-7);
 
m2=2700*(pi/4*(14.3e-3)^2)*422e-3;
% m2=0;
Iz2=m2*L2^2/12;

m3=275.166e-3;
Iz3=4831958.5e-9+(15+142)*(343/8)^2*10^(-7);

m4=2700*(pi/4*(14.3e-3)^2)*422e-3;
% m4=0;
Iz4=m4*L4^2/12;

m5=0.65;
Iz5=m5*L5^2/12;
% Iz5=0;

m6=143e-3;
Iz6=m6*L6^2/12;

m7=134e-3;
Iz7=m7*L7^2/12;

m8=134e-3;
Iz8=m8*L8^2/12;

% Distâncias dos centros de massa dos elos (1), (2), (4) e (5) (usualmente
% metade do comprimento do elo
r1=L1/2;
r2=L2/2;
r3=L3/2;
r4=L4/2;
% r5=L5/2;
% r6=L6/2;
r7=L7/2;
r8=L8/2;

% Ângulo de inclinação do efetuador thetaB
thetaB=0;

% Matriz Kq
Amat=pi/4*(14.3e-3)^2;
Emat=69e9;
Imat=(30.7e-3)^3*(6.1e-3)/12;

% rigidez conjunto motoredutor
kred=859.43;
k1=(kred*3*Emat*Imat/L1*fator)/(kred+3*Emat*Imat/L1*fator);
k2=(kred*3*Emat*Imat/L1*fator)/(kred+3*Emat*Imat/L1*fator);

% k1=3*Emat*Imat/L1*fator;
% k2=3*Emat*Imat/L3;
% k1=((3*Emat*Imat/L1)^(-1)+(4.89)^(-1)+(2.80)^(-1))^(-1)*fator;
% k2=((3*Emat*Imat/L3)^(-1)+(4.89)^(-1)+(2.80)^(-1))^(-1);
k3=Emat*Amat/L1;
k4=Emat*Amat/L3;
k5=Emat*Amat/L2;
k6=Emat*Amat/L4;

% % Mola k7 provisória
% k7=100000000000;

% % Matriz de massa dos elos ativos 1 e 1'
% m1=275.166e-3;
% m1_1=275.166e-3;
% IG1=4831958.5e-9;
% IG1_1=4831958.5e-9;
% l1=222.7e-3;
% l1_1=222.7e-3;
% IO1=IG1+m1*l1^2;
% IO1_1=IG1_1+m1_1*l1_1^2;
% 
% % Parâmetros de inércia das barras 2 e 2'
% m2=2700*(pi/4*(14.3e-3)^2)*422e-3;
% m2_1=2700*(pi/4*(14.3e-3)^2)*422e-3;
% IG2=m2*L2^2/12;
% IG2_1=m2_1*L2^2/12;

% Criação da malha do espaço de trabalho a ser mapeado
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[XB,YB]=meshgrid(x,y);

% Inicialização das matrizes com o tamanho correto das características a 
% serem mapeadas 
deslx=zeros(size(XB));
desly=zeros(size(XB));
omega=zeros(size(XB));
autovecx=zeros(size(XB));
autovecy=zeros(size(XB));
numcond=zeros(size(XB));

% Determinação das características em cada ponto do espaço de trabalho a
% ser mapeado  
for k=1:size(XB,1)
for p=1:size(YB,2)
% Variáveis locais xB e yB apenas para a posição dentro do loop for
xB=XB(k,p);
yB=YB(k,p);
    
% coeficientes do sistema de equações não lineares da cinemática inversa
% G*sin(q1)+H*cos(q1)+I=0
% J*sin(q2)+K*cos(q2)+L=0

% % Antiga
% G=-2*L1*yB;
% H=-2*L1*xB+2*L1*L0-2*L1*L3;
% I=xB^2-2*L0*xB+2*L3*xB+L1^2-2*L3*L0+L3^2-L2^2+L0^2+yB^2;
% J=-2*L1*yB;
% K=-2*L1*xB-2*L1*L0+2*L1*L3;
% L=xB^2+2*L0*xB-2*L3*xB+L1^2+L0^2-2*L0*L3+L3^2-L2^2+yB^2;

% % Nova nomenclatura 2gdl Equações originais
% % L0^2+L1^2-L2^2-2*L0*L5+L5^2-2*L0*xB+2*L5*xB+xB^2+yB^2+2*L1*(L0-L5-xB)*cosq1-2*L1*yB*sinq1
% % L0^2+L3^2-L4^2-2*L0*L5+L5^2+2*L0*xB-2*L5*xB+xB^2+yB^2+2*L3*(L0-L5+xB)*cosq2-2*L3*yB*sinq2
% G=-2*L1*yB;
% H=2*L1*(L0-L5-xB);
% I=L0^2+L1^2-L2^2-2*L0*L5+L5^2-2*L0*xB+2*L5*xB+xB^2+yB^2;
% J=-2*L3*yB;
% K=2*L3*(L0-L5+xB);
% L=L0^2+L3^2-L4^2-2*L0*L5+L5^2+2*L0*xB-2*L5*xB+xB^2+yB^2;

% Nova nomenclatura 3gdl Equações originais
% L0^2-L2^2-2*L0*xB+xB^2+yB^2+L1^2-2*L0*L5*cos(thetaB)+2*L5*xB*cos(thetaB)+L5^2+cosq1*(2*L0*L1-2*L1*xB-2*L1*L5*cos(thetaB))+2*L5*yB*sin(thetaB)+sinq1*(-2*L1*yB-2*L1*L5*sin(thetaB))
% L0^2-L4^2+2*L0*xB+xB^2+yB^2+L3^2-2*L0*L5*cos(thetaB)-2*L5*xB*cos(thetaB)+L5^2+cosq2*(2*L0*L3+2*L3*xB-2*L3*L5*cos(thetaB))-2*L5*yB*sin(thetaB)+sinq2*(-2*L3*yB+2*L3*L5*sin(thetaB))
G=(-2*L1*yB-2*L1*L5*sin(thetaB));
H=(2*L0*L1-2*L1*xB-2*L1*L5*cos(thetaB));
I=L0^2-L2^2-2*L0*xB+xB^2+yB^2+L1^2-2*L0*L5*cos(thetaB)+2*L5*xB*cos(thetaB)+L5^2+2*L5*yB*sin(thetaB);
J=(-2*L3*yB+2*L3*L5*sin(thetaB));
K=(2*L0*L3+2*L3*xB-2*L3*L5*cos(thetaB));
L=L0^2-L4^2+2*L0*xB+xB^2+yB^2+L3^2-2*L0*L5*cos(thetaB)-2*L5*xB*cos(thetaB)+L5^2-2*L5*yB*sin(thetaB);


% usando as identidades trigonométricas
% u=tan(teta/2)
% cos(teta)=(1-u^2)/(1+u^2)
% sin(teta)=(2*u)/(1+u^2)
% pode-se encontrar q1 e q2
% fazendo esta substituição no membro do lado esquerdo
% (I-H)*u^2+(2*G)*u+(H+I)=0
% para o membro do lado direito
% (L-K)*v^2+(2*J)*v+(K+L)=0
% resolvendo as duas equações tem-se para o membro esquerdo
u1=(-G+sqrt(G^2-I^2+H^2))/(I-H);
u2=(-G-sqrt(G^2-I^2+H^2))/(I-H);

% para o membro do lado direito
v1=(-J+sqrt(J^2-L^2+K^2))/(L-K);
v2=(-J-sqrt(J^2-L^2+K^2))/(L-K);

% retomando a definição das identidades trigonométricas, tem-se para o
% membro do lado esquerdo
cosq11=(1-u1^2)/(1+u1^2);
sinq11=(2*u1)/(1+u1^2);
% q11=atand((2*u1)/(1-u1^2));
q11=atand((sinq11)/(cosq11));

cosq12=(1-u2^2)/(1+u2^2);
sinq12=(2*u2)/(1+u2^2);
% q12=atand((2*u2)/(1-u2^2));
q12=atand((sinq12)/(cosq12));

% para o lado direito
cosq21=(1-v1^2)/(1+v1^2);
sinq21=(2*v1)/(1+v1^2);
% q21=atand((2*v1)/(1-v1^2));
q21=atand((sinq21)/(cosq21));

cosq22=(1-v2^2)/(1+v2^2);
sinq22=(2*v2)/(1+v2^2);
% q22=atand((2*v2)/(1-v2^2));
q22=atand((sinq22)/(cosq22));

% Solução escolhida para a cinemática inversa
sinq1=sinq12;
cosq1=cosq12;
sinq2=sinq22;
cosq2=cosq22;
q1=q12;
q2=q22;

% pontos para plotagem do grafico

% pontos O1 e O2
O1=[L0 0];
O2=[-L0 0];

% % pontos B1 e B2 2gdl
% B1=[xB+L5 yB];
% B2=[xB-L5 yB];

% pontos B1 e B2 3gdl
B1=[xB+L5*cos(thetaB) yB+L5*sin(thetaB)];
B2=[xB-L5*cos(thetaB) yB-L5*sin(thetaB)];


% % ponto escolhido (configurações das cadeia periféricas para fora)2gdl
% A1=[L1*cosq1+L0 L1*sinq1];
% A2=[-L3*cosq2-L0 L3*sinq2];

% ponto escolhido (configurações das cadeia periféricas para fora)3gdl
A1=[L1*cosq1+L0 L1*sinq1];
A2=[-L3*cosq2-L0 L3*sinq2];

% Matriz jacobiana
% Jx*dx=Jq*dq
% Matriz Jx

% % Antiga
% % Jx=[Q3(i,j)-L3+L0+L1*cosq2 Q4(i,j)-L1*sinq2 ; Q3(i,j)+L3-L0-L1*cosq1 Q4(i,j)-L1*sinq1];
% A=-2*L1*cosq1-2*L0+2*L3;
% B=-2*L1*sinq1;
% D=2*L0-2*L3-2*L1*cosq2;
% E=-2*L1*sinq2;
% Jx=[2*xB+A 2*yB+B;2*xB+D 2*yB+E];

% % Jx Nova 2gdl
% Jx=[2*(-L0+L5+xB-L1*cosq1) 2*(yB-L1*sinq1);
%     2*(L0-L5+xB+L3*cosq2)  2*(yB-L3*sinq2)];

% Jx Nova 3gdl
Jx=zeros(3,3);
Jx(1,1)=2*(-L0+xB-L1*cosq1+L5*cos(thetaB)); 
Jx(1,2)=2*(yB-L1*sinq1+L5*sin(thetaB));
Jx(1,3)=2*L5*cos(thetaB)*(yB-L1*sinq1+L5*sin(thetaB))-2*L5*(-L0+xB-L1*cosq1+L5*cos(thetaB))*sin(thetaB);
Jx(2,1)=2*(L0+xB+L3*cosq2-L5*cos(thetaB));
Jx(2,2)=2*(yB-L3*sinq2-L5*sin(thetaB));
Jx(2,3)=2*L5*(L0+xB+L3*cosq2-L5*cos(thetaB))*sin(thetaB)-2*L5*cos(thetaB)*(yB-L3*sinq2-L5*sin(thetaB));
Jx(3,1)=0;
Jx(3,2)=0;
Jx(3,3)=1;

% % Matriz Jq
% % J1=(Q3(i,j)+L3-L0-L1*cosq1)*(L1*sinq1)+(Q4(i,j)-L1*sinq1)*(L1*cosq1);
% % J2=(Q3(i,j)-L3+L0+L1*cosq2)*(L1*sinq21)+(Q4(i,j)-L1*sinq2)*(L1*cosq2);
% % Antiga
% J1=G*cosq1-H*sinq1;
% J2=J*cosq2-K*sinq2;
% J3=2*(L1*cosq1+L0-xB-L3)*cosq1+2*(L1*sinq1-yB)*sinq1;
% J4=2*(L1*cosq2-L0-xB+L3)*cosq2+2*(L1*sinq2-yB)*sinq2;
% Jq=-1*[J1 0 J3 0 -2*L2 0;0 J2 0 J4 0 -2*L2];


% Jq nova 2gdl
% Jq=[-2*L1*(-L0+L5+xB-L1*cosq1)*sinq1+2*L1*cosq1*(yB-L1*sinq1) 0;
%     0                                                         2*L3*(L0-L5+xB+L3*cosq2)*sinq2+2*L3*cosq2*(yB-L3*sinq2)];

% Jq nova 3gdl
Jq=zeros(3,3);
Jq(1,1)=2*L1*cosq1*(yB-L1*sinq1+L5*sin(thetaB))-2*L1*(-L0+xB-L1*cosq1+L5*cos(thetaB))*sinq1; 
Jq(1,2)=0;
Jq(1,3)=0;
Jq(2,1)=0;
Jq(2,2)=2*L3*(L0+xB+L3*cosq2-L5*cos(thetaB))*sinq2+2*L3*cosq2*(yB-L3*sinq2-L5*sin(thetaB));
Jq(2,3)=0;
Jq(3,1)=0;
Jq(3,2)=0;
Jq(3,3)=1;

% % Jqext nova 2 gdl
% Jqext=[-2*L1*(-L0+L5+xB-L1*cosq1)*sinq1+2*L1*cosq1*(yB-L1*sinq1) 0                                                       2*cosq1*(-L0+L5+xB-L1*cosq1)+2*sinq1*(yB-L1*sinq1)  0                                                 2*L2 0;
%         0                                                        2*L3*(L0-L5+xB+L3*cosq2)*sinq2+2*L3*cosq2*(yB-L3*sinq2) 0                                                  -2*cosq2*(L0-L5+xB+L3*cosq2)+2*sinq2*(yB-L3*sinq2) 0    2*L4];


% Jqext nova 3 gdl  
Jqext=zeros(3,7);
Jqext(1,1)=2*L1*cosq1*(yB-L1*sinq1+L5*sin(thetaB))-2*L1*(-L0+xB-L1*cosq1+L5*cos(thetaB))*sinq1; 
Jqext(1,2)=0;
Jqext(1,3)=2*cosq1*(-L0+xB-L1*cosq1+L5*cos(thetaB))+2*sinq1*(yB-L1*sinq1+L5*sin(thetaB));
Jqext(1,4)=0;
Jqext(1,5)=2*L2;
Jqext(1,6)=0;
Jqext(1,7)=0;
Jqext(2,1)=0; 
Jqext(2,2)=2*L3*(L0+xB+L3*cosq2-L5*cos(thetaB))*sinq2+2*L3*cosq2*(yB-L3*sinq2-L5*sin(thetaB));
Jqext(2,3)=0;
Jqext(2,4)=2*sinq2*(yB-L3*sinq2-L5*sin(thetaB))-2*cosq2*(L0+xB+L3*cosq2-L5*cos(thetaB));
Jqext(2,5)=0;
Jqext(2,6)=2*L4;
Jqext(2,7)=0;
Jqext(3,1)=0; 
Jqext(3,2)=0;
Jqext(3,3)=0;
Jqext(3,4)=0;
Jqext(3,5)=0;
Jqext(3,6)=0;
Jqext(3,7)=1;

% Rigidez da cadeia passiva k7
% ângulo fi do paralelogramo
cosfi=xB/L7;
sinfi=sqrt(L7^2-xB^2)/L7;
tanfi=sinfi/cosfi;
% % Rigidez original
% k7=(1/(2*L6*tanfi^2*A6*Emat)+L7/(4*L6^2*sinfi^2*A7*Emat)+L8/(4*L6^2*cosfi^2*A8*Emat)+L6^5/(6*Emat*I6))^(-1);

% Rigidez simplificada
% k7=(L7/(4*L6^2*sinfi^2*A7*Emat)+L8/(4*L6^2*sinfi^2*A8*Emat)+L6^5/(6*Emat*I6))^(-1);
k7=(L7/(4*L6^2*sinfi^2*A7*Emat)+L8/(4*L6^2*sinfi^2*A8*Emat)+L6/(6*Emat*I6))^(-1);

% % Matriz de rigidez do espaço dos atuadores 2 gdl
% Kq=[k1 0 0 0 0 0;
%     0 k2 0 0 0 0;
%     0 0 k3 0 0 0;
%     0 0 0 k4 0 0;
%     0 0 0 0 k5 0;
%     0 0 0 0 0 k6];

% Matriz de rigidez do espaço dos atuadores 3 gdl
Kq=[k1 0  0  0  0  0  0;
    0  k2 0  0  0  0  0;
    0  0  k3 0  0  0  0;
    0  0  0  k4 0  0  0;
    0  0  0  0  k5 0  0;
    0  0  0  0  0  k6 0;
    0  0  0  0  0  0  k7];
% Matriz de compliância para evitar a inversão de uma matriz retangular
Comp=inv(Jx)*Jqext*inv(Kq)*transpose(Jqext)*transpose(inv(Jx));

% % Deslocamento
% delta=Comp*[0;1;0];
% 
% % dados para plotagem da superfície
% deslx(k,p)=delta(1,1);
% desly(k,p)=delta(2,1);

% dados para plotagem da superfície 

% kinverso=inv(Comp);
% deslx(k,p)=kinverso(1,1);
% desly(k,p)=kinverso(2,2);

% Rigidez pela inversão dos elementos da matriz compliância
deslx(k,p)=1/Comp(1,1);
desly(k,p)=1/Comp(2,2);

% % Para a plotagem da flexibilidade (compliância)
% deslx(k,p)=Comp(1,1);
% desly(k,p)=Comp(2,2);



% % nova formulação de rigidez para compatibilidade com os dados
% % experimentais na tese
% deslx(k,p)=1/sqrt(Comp(1,1)^2+Comp(2,1)^2);
% desly(k,p)=1/sqrt(Comp(1,2)^2+Comp(2,2)^2);






% % ANTIGA FORMULAÇÃO
% % Cálculo da matriz de rigidez do sistema via conceito de massa concentrada
% % do sistema
% 
% % Matriz de massa dos elos ativos
% 
% % Matriz jacobiana para o cálculo da matriz de massa das barras em rotação
% Jm=(-1*[J1 0;0 J2])\Jx;
% 
% % Matriz de massa dos elos ativos
% M_elos_at=[IO1*Jm(1,1)^2+IO1_1*Jm(2,1)^2              IO1*Jm(1,1)*Jm(1,2)+IO1_1*Jm(2,1)*Jm(2,2);
%            IO1*Jm(1,2)*Jm(1,1)+IO1_1*Jm(2,2)*Jm(2,1)  IO1*Jm(1,2)^2+IO1_1*Jm(2,2)^2];
% % meq1=IG1/L1+m1/4;
% % meq1_1=IG1_1/L1+m1_1/4;
% 
% 
% % Matriz de massa das duas barras em rototranslação 2 e 2', dividindo sua
% % massa entre as duas extremidades do elo
% 
% % Matriz de massa dos elos passivos, primeiro derivando os termos da
% % equação de Lagrange
% dT2ddotq3(1)=(0.25*m2*((1-Jm(1,1)*L1*sinq12)^2+(Jm(1,1)*L1*cosq12)^2)+IG2*(Jm(1,1)*L1*cosq12)^2/(L0+L1*cosq12-XB(k,p)-L3)^2);
% dT2ddotq3(2)=(0.25*m2*((-Jm(1,2)*L1*sinq12)*(1-Jm(1,1)*L1*sinq12)+(1+Jm(1,2)*L1*cosq12)*(Jm(1,1)*L1*cosq12))+IG2*(Jm(1,2)*L1*cosq12-1)*(Jm(1,1)*L1*cosq12)/(L0+L1*cosq12-XB(k,p)-L3)^2);
% 
% dT2ddotq4(1)=(0.25*m2*((-Jm(1,2)*L1*sinq12)*(1-Jm(1,1)*L1*sinq12)+(1+Jm(1,2)*L1*cosq12)*(Jm(1,1)*L1*cosq12))+IG2*(Jm(1,2)*L1*cosq12-1)*(Jm(1,1)*L1*cosq12)/(L0+L1*cosq12-XB(k,p)-L3)^2);
% dT2ddotq4(2)=(0.25*m2*((-Jm(1,2)*L1*sinq12)^2+(1+Jm(1,2)*L1*cosq12)^2)+IG2*(Jm(1,2)*L1*cosq12-1)^2/(L0+L1*cosq12-XB(k,p)-L3)^2);
% 
% dT2_1ddotq3(1)=(0.25*m2_1*((1-Jm(2,1)*L1*sinq21)^2+(Jm(2,1)*L1*cosq21)^2)+IG2_1*(Jm(2,1)*L1*cosq21)^2/(-L0+L1*cosq21-XB(k,p)+L3)^2);
% dT2_1ddotq3(2)=(0.25*m2_1*((-Jm(2,2)*L1*sinq21)*(1-Jm(2,1)*L1*sinq21)+(1+Jm(2,2)*L1*cosq21)*(Jm(2,1)*L1*cosq21))+IG2_1*(Jm(2,2)*L1*cosq21-1)*(Jm(2,1)*L1*cosq21)/(-L0+L1*cosq21-XB(k,p)+L3)^2);
% 
% dT2_1ddotq4(1)=(0.25*m2_1*((-Jm(2,2)*L1*sinq21)*(1-Jm(2,1)*L1*sinq21)+(1+Jm(2,2)*L1*cosq21)*(Jm(2,1)*L1*cosq21))+IG2_1*(Jm(2,2)*L1*cosq21-1)*(Jm(2,1)*L1*cosq21)/(-L0+L1*cosq21-XB(k,p)+L3)^2);
% dT2_1ddotq4(2)=(0.25*m2_1*((-Jm(2,2)*L1*sinq21)^2+(1+Jm(2,2)*L1*cosq21)^2)+IG2_1*(Jm(2,2)*L1*cosq21-1)^2/(-L0+L1*cosq21-XB(k,p)+L3)^2);
% 
% % montagem da matriz de massa dos elos passivos
% M_elos_pas=[dT2ddotq3(1)+dT2_1ddotq3(1) dT2ddotq3(2)+dT2_1ddotq3(2);
%             dT2ddotq4(1)+dT2_1ddotq4(1) dT2ddotq4(2)+dT2_1ddotq4(2)];
% % M_elos_pas=0.5*[m2*L1^2*Jm(1,1)^2+m2_1*L1^2*Jm(2,1)^2+(m2+m2_1)    m2*L1^2*Jm(1,1)*Jm(1,2)+m2_1*L1^2*Jm(2,1)*Jm(2,2);
% %                 m2*L1^2*Jm(1,2)*Jm(1,1)+m2_1*L1^2*Jm(2,2)*Jm(2,1)  m2*L1^2*Jm(1,2)^2+m2_1*L1^2*Jm(2,2)^2+(m2+m2_1)];
% 
% 
% % Matriz de massa
% Mx=M_elos_at+M_elos_pas;
% % matmassa=[meqx 0;
% %           0 meqy];
%   


% NOVA FORMULAÇÃO
J=inv(Jq)*Jx;
% J=Jq\Jx;

% % Matriz Cx1 gdl 2
% Cx1=[-r1*sinq1 0;
%       r1*cosq1 0;
%       1        0]*J;
% % Matriz Cx1 gdl 2   
% Cx3=[0 r3*sinq2;
%      0 r3*cosq2;
%      0 -1]*J;

% Matriz Cx1 gdl 3
Cx1=[-r1*sinq1 0 0;
      r1*cosq1 0 0;
      1        0 0]*J;
% Matriz Cx1 gdl 3   
Cx3=[0 r3*sinq2 0;
     0 r3*cosq2 0;
     0 -1       0]*J; 
 
% % Matriz Cx2 gdl 2
% Cx2=[r2/L2-((L2-r2)*(-L0+L5+xB-L1*cosq1)*sinq1)/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))  (-(((L2-r2)*sinq1*(yB-L1*sinq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))));
%      ((L2-r2)*cosq1*(-L0+L5+xB-L1*cosq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))       r2/L2+((L2-r2)*cosq1*(yB-L1*sinq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1));
%      ((-1*(L0-L5-xB+L1*cosq1)*sinq1)/(yB*cosq1+(L0-L5-xB)*sinq1)+1)/(L1*sinq1-yB)   (((sinq1*(yB-L1*sinq1))/((L1*sinq1-yB)*(yB*cosq1+(L0-L5-xB)*sinq1))))];
% 
% % Matriz Cx4 gdl 2
% Cx4=[r4/L4+((L4-r4)*(L0-L5+xB+L3*cosq2)*sinq2)/(L4*(yB*cosq2+(L0-L5+xB)*sinq2))  ((L4-r4)*sinq2*(yB-L3*sinq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2));
%      ((L4-r4)*cosq2*(L0-L5+xB+L3*cosq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2))      r4/L4+((L4-r4)*cosq2*(yB-L3*sinq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2));
%      (((-L0+L5-xB-L3*cosq2)*sinq2)/(yB*cosq2+(L0-L5+xB)*sinq2)+1)/(L3*sinq2-yB)  (-((sinq2*(yB-L3*sinq2))/((L3*sinq2-yB)*(yB*cosq2+(L0-L5+xB)*sinq2))))];

% Matriz Cx2 gdl 3
Cx2=zeros(3,3);
Cx2(1,1)=r2/L2-((L2-r2)*(-L0+xB-L1*cosq1+L5*cos(thetaB))*sinq1)/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))); 
Cx2(1,2)=-(((L2-r2)*sinq1*(yB-L1*sinq1+L5*sin(thetaB)))/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))));
Cx2(1,3)=(-((L5*r2*sin(thetaB))/L2)-(L5*(L2-r2)*sinq1*(yB*cos(thetaB)-L1*(sinq1*cos(thetaB)-sin(thetaB)*cosq1)+(L0-xB)*sin(thetaB)))/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))));
Cx2(2,1)=((L2-r2)*cosq1*(-L0+xB-L1*cosq1+L5*cos(thetaB)))/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1)));
Cx2(2,2)=(r2/L2+((L2-r2)*cosq1*(yB-L1*sinq1+L5*sin(thetaB)))/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))));
Cx2(2,3)=((L5*r2*cos(thetaB))/L2+(L5*(L2-r2)*cosq1*(yB*cos(thetaB)-L1*(sinq1*cos(thetaB)-sin(thetaB)*cosq1)+(L0-xB)*sin(thetaB)))/(L2*(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))));
Cx2(3,1)=(((-L0+xB-L1*cosq1+L5*cos(thetaB))*sinq1)/(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))+1)/(-yB+L1*sinq1-L5*sin(thetaB));
Cx2(3,2)=(sinq1*(yB-L1*sinq1+L5*sin(thetaB)))/((yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))*(-yB+L1*sinq1-L5*sin(thetaB)));
Cx2(3,3)=((L5*sinq1*(yB*cos(thetaB)-L1*(sinq1*cos(thetaB)-sin(thetaB)*cosq1)+(L0-xB)*sin(thetaB)))/(yB*cosq1+(L0-xB)*sinq1-L5*(sinq1*cos(thetaB)-sin(thetaB)*cosq1))-L5*sin(thetaB))/(-yB+L1*sinq1-L5*sin(thetaB));

% Matriz Cx4 gdl 3
Cx4=zeros(3,3);
Cx4(1,1)=(r4/L4+((L4-r4)*(L0+xB+L3*cosq2-L5*cos(thetaB))*sinq2)/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))); 
Cx4(1,2)=((L4-r4)*sinq2*(yB-L3*sinq2-L5*sin(thetaB)))/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)));
Cx4(1,3)=((L5*r4*sin(thetaB))/L4+(L5*(L4-r4)*sinq2*(-yB*cos(thetaB)+(L0+xB)*sin(thetaB)+L3*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2))));
Cx4(2,1)=((L4-r4)*cosq2*(L0+xB+L3*cosq2-L5*cos(thetaB)))/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)));
Cx4(2,2)=(r4/L4+((L4-r4)*cosq2*(yB-L3*sinq2-L5*sin(thetaB)))/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2))));
Cx4(2,3)=((L5*(L4-r4)*cosq2*(-yB*cos(thetaB)+(L0+xB)*sin(thetaB)+L3*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))/(L4*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))-(L5*r4*cos(thetaB))/L4);
Cx4(3,1)=(((-L0-xB-L3*cosq2+L5*cos(thetaB))*sinq2)/(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2))+1)/(-yB+L3*sinq2+L5*sin(thetaB));
Cx4(3,2)=-((sinq2*(yB-L3*sinq2-L5*sin(thetaB)))/((-yB+L3*sinq2+L5*sin(thetaB))*(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2))));
Cx4(3,3)=(L5*sin(thetaB)-(L5*sinq2*(-yB*cos(thetaB)+(L0+xB)*sin(thetaB)+L3*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))/(yB*cosq2+(L0+xB)*sinq2-L5*(sinq2*cos(thetaB)+sin(thetaB)*cosq2)))/(-yB+L3*sinq2+L5*sin(thetaB));

Cx5=[1 0 0;
     0 1 0;
     0 0 1];
  
Cx6=[1          0 0;
     -1/(tanfi) 0 0;
     0          0 0];

Cx7=[r7/L7          0 0;
     -r7/(L7*tanfi) 0 0;
     -1/(L7*sinfi)  0 0];
 
Cx8=[r8/L8          0 0;
     -r8/(L8*tanfi) 0 0;
     -1/(L8*sinfi)  0 0];

M1=[m1 0   0;
    0   m1 0;
    0   0  Iz1];
M2=[m2 0   0;
    0   m2 0;
    0   0  Iz2];
M3=[m3 0   0;
    0   m3 0;
    0   0  Iz3];
M4=[m4 0   0;
    0   m4 0;
    0   0  Iz4];

M5=[m5 0   0;
     0   m5 0;
     0   0   Iz5]; 
 
M6=[m6 0   0;
     0   m6 0;
     0   0   Iz6]; 

M7=[m7 0   0;
     0   m7 0;
     0   0   Iz7]; 

M8=[m8 0   0;
     0   m8 0;
     0   0   Iz8];

Mx=Cx1'*M1*Cx1+Cx2'*M2*Cx2+Cx3'*M3*Cx3+Cx4'*M4*Cx4+Cx5'*M5*Cx5+Cx6'*M6*Cx6+Cx7'*M7*Cx7+Cx8'*M8*Cx8;
% Mx=Cx1'*M1*Cx1+Cx3'*M3*Cx3;

 
 
% Matriz de rigidez a partir da matriz compliância
Kx=inv(Comp);

% Número de condicionamento da matriz de rigidez 2x2
numcond(k,p)=cond([Kx(1,1) Kx(1,2);
                   Kx(2,1) Kx(2,2)]);

% Cálculo da primeira frequência natural a partir das matrizes de massa e 
% de rigidez e ordenando em ordem crescente as frequências
[V,W] = eig(Mx\Kx);
[d,ind] = sort(diag(W));
Wordenado=W(ind,ind);
Vordenado=V(:,ind);
omega1(k,p)=(sqrt(Wordenado(1,1)))/(2*pi);
omega2(k,p)=(sqrt(Wordenado(2,2)))/(2*pi);
% Autovetores 
% autovecx(k,p)=Vordenado(1,1);
% autovecy(k,p)=Vordenado(2,1);
autovecx1(k,p)=Vordenado(1,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2+Vordenado(3,1)^2);
autovecy1(k,p)=Vordenado(2,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2+Vordenado(3,1)^2);

autovecx2(k,p)=Vordenado(1,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2);
autovecy2(k,p)=Vordenado(2,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2);
% [W] = eig(matmassa\Ktotal);
% omega(k,p)=(sqrt(min(W)))/(2*pi);
end
end


% Dados experimentais mapeamentos

stiffxexptransp = [ ...
              135         410         431         615; ...
              183         537         767        1029; ...
              217        1587        2106        1654; ...
              247        1950       31585       12148; ...
%               247        1950       0       12148; ...
              341        7896       26321        5264; ...
%               341        7896       0        5264; ...
              423        3948        3899        8312; ...
              363        1564        1449        2871; ...
              292         659         845         921; ...
              423         682         726        1654; ...
              560         468         565        1136; ...
             1174         600         651        1533; ...
             4327         882         890        1187; ...
            10892       31585        1962        1429; ...
%             10892       0        1962        1429; ...
             1587        1026        2078        4387; ...
             1385         996        1469        1174; ...
             1093         733         629         377; ...
              810         445         371         267; ...
              536         321         266         265; ...
              341         231         311         572; ...
              378         368         776        1815; ...
              389         539        1416        2771; ...
              675        1774        1815        2272; ...
             1322        1880        2632        3008; ...
             1556        2149        3037        3509; ...
             1564        2038        6074        5353; ...
      ...
     ];
stiffxexp=stiffxexptransp';

stiffyexptransp = [ ...
              276         900         267         382; ...
              331        1579         251         519; ...
              439         612         282         522; ...
              534         492         356         345; ...
              863         496         442         363; ...
            21057         393         606         398; ...
%             0         393         606         398; ...
            45122         347         752         404; ...
%             0         347         752         404; ...
              926         494         671         420; ...
             2106         629        8774         369; ...
             2289        3673         788         484; ...
             1436        4859        1836         658; ...
             1938        3549        3998        1436; ...
             2149       78964        1926        1373; ...
%              2149       0       1926        1373; ...
             5640        1224        1350        1398; ...
            10892        1201        1192         685; ...
            10189         825         584         307; ...
             1504         512         802         207; ...
              649         390        1533         179; ...
              297         228         887         217; ...
              241         291        1815         363; ...
              281         360         699         480; ...
              305         532         391         442; ...
              905        2012         362         500; ...
             1196        2820         341         508; ...
             1268        3223         314         539; ...
      ...
     ];
stiffyexp=stiffyexptransp';

freqxexptranp = [ ...
        15.8000   14.5000   13.5000   10.7000; ...
        15.5000   15.5000   14.5000   13.0000; ...
        16.5000   15.8000   15.5000   14.3000; ...
        17.5000   15.8000   15.7000   15.7000; ...
        17.3000   16.2000   17.0000   14.5000; ...
        17.0000   16.0000   16.5000   15.7000; ...
        18.3000   16.7000   15.5000   16.2000; ...
        18.2000   17.0000   16.5000   14.8000; ...
        19.3000   17.0000   16.2000   15.8000; ...
        19.5000   17.0000   15.7000   17.0000; ...
        19.8000   17.0000   17.2000   15.3000; ...
        19.5000   16.2000   14.2000   16.0000; ...
        19.3000   14.5000   12.8000   10.5000; ...
        16.7000   14.5000   13.0000   11.3000; ...
        15.7000   17.0000   16.0000   14.8000; ...
        17.5000   16.0000   12.7000   14.3000; ...
        17.2000   15.0000   15.7000   16.2000; ...
        16.2000   17.5000   17.0000   16.2000; ...
        18.7000   17.5000   17.2000   17.0000; ...
        19.2000   17.7000   17.5000   16.8000; ...
        19.0000   17.8000   17.3000   16.3000; ...
        19.0000   17.8000   17.0000   15.5000; ...
        18.0000   17.5000   16.7000   15.5000; ...
        16.8000   17.0000   15.3000   15.2000; ...
        17.3000   16.8000   15.7000   14.5000; ...
      ...
     ];
freqxexp=freqxexptranp';

freqyexptranp = [ ...
        13.5000   12.5000   13.0000   13.0000; ...
        16.5000   15.8000   13.2000   13.5000; ...
        16.5000   16.3000   15.5000   13.7000; ...
        17.5000   16.7000   16.0000   15.5000; ...
        17.5000   17.3000    5.5000   15.7000; ...
         4.8000    4.8000    5.2000   16.5000; ...
         5.7000    4.7000    5.2000    5.2000; ...
         4.5000    4.7000    5.2000    5.2000; ...
         4.5000    4.5000    5.5000    5.2000; ...
         6.2000    5.0000    6.0000    5.0000; ...
         5.8000    4.8000    5.5000    4.8000; ...
         6.5000    7.2000    9.8000    5.5000; ...
         7.7000    5.3000    9.3000    6.2000; ...
         5.7000    5.7000    5.3000    4.5000; ...
         4.8000    5.3000    4.8000    4.5000; ...
         4.7000    5.8000    5.0000    6.3000; ...
         5.5000    5.0000    5.3000    5.8000; ...
         5.2000    5.3000    4.8000    5.8000; ...
         4.7000    5.3000    4.7000    5.0000; ...
         5.2000    5.3000    5.0000    5.5000; ...
         6.0000    5.0000   17.0000   16.5000; ...
         5.7000   16.0000   16.5000   15.0000; ...
        16.5000   15.5000   15.8000   14.5000; ...
        16.0000   14.2000   14.3000   13.7000; ...
        14.3000   13.5000   14.0000   13.5000; ...
      ...
     ];
freqyexp=freqyexptranp';

% Novos dados experimentais
    exatidaoYcarga0transp = [ ...
        10.7300   17.3600   14.4500   16.7400; ...
        10.3200   15.6500   14.0100   15.0400; ...
         9.9300   15.3400   13.8600   14.7100; ...
         9.5100   14.8300   13.3800   14.5100; ...
         7.6600   12.4600   11.2200   11.6900; ...
         5.0900    9.5500    9.0600    9.1900; ...
         3.5300    7.7100    7.2500    7.9900; ...
         1.9000    6.4600    6.1800    6.7800; ...
         0.8500    4.1400    3.5200    5.3800; ...
         0.2400    1.2400    0.7800    3.6600; ...
        -0.1000    0.5400   -0.4800    2.6700; ...
        -0.2700    0.4900   -0.8200    1.5900; ...
        -1.8100   -1.8600   -1.7400         0; ...
        -3.2100   -4.2900   -2.5500   -1.7700; ...
        -3.4700   -4.4600   -2.9200   -2.6400; ...
        -4.0400   -5.0800   -4.1900   -5.4500; ...
        -4.4100   -7.2300   -6.8900   -8.8900; ...
        -5.1100   -9.6800   -9.9600  -11.3800; ...
        -8.3000  -13.4100  -11.8300  -11.9200; ...
       -11.9500  -15.0900  -12.5200  -12.6300; ...
       -12.4800  -15.3500  -13.5100  -14.1300; ...
       -12.9700  -15.5900  -14.5400  -15.4600; ...
       -13.4400  -15.9700  -14.8700  -15.8000; ...
       -13.9600  -16.3800  -15.2700  -15.8900; ...
       -13.8000  -16.1500  -14.8200  -15.8000; ...
      ...
     ];
exatidaoYcarga0=exatidaoYcarga0transp';


    exatidaoZcarga0transp = [ ...
       -10.1000  -11.3000  -10.7000  -10.6000; ...
       -10.6700   -9.7200  -10.0600   -9.3900; ...
       -10.7200  -10.0600  -10.2200   -9.5000; ...
       -11.0600  -10.4400  -10.2800   -9.8900; ...
       -11.0500  -10.5000  -10.7500  -10.4000; ...
        -9.5000  -10.0000  -10.2200  -10.7200; ...
        -8.3300   -9.5600   -9.5000   -9.1700; ...
        -4.4400   -6.0000   -8.0000   -7.1100; ...
        -2.9400   -4.7500   -6.4000   -8.1500; ...
        -1.0000   -1.8300   -4.5000   -6.6100; ...
        -0.6100   -1.0000   -2.3300   -5.3300; ...
        -1.0000   -1.0600   -1.6100   -3.5000; ...
        -0.4500   -1.2000   -2.1000   -3.5000; ...
        -0.3300   -1.8300   -2.7800   -3.2200; ...
        -0.7200   -2.0000   -3.3300   -5.1700; ...
        -1.3300   -3.5000   -6.2200  -10.3300; ...
        -3.4500   -6.3500   -9.2500  -13.7000; ...
        -6.3300   -8.9400  -13.6700  -15.2200; ...
       -10.6100  -14.1100  -15.2200  -15.0600; ...
       -15.2200  -15.5000  -14.7200  -14.6100; ...
       -14.9000  -14.8000  -14.0000  -14.3000; ...
       -14.3900  -14.1700  -13.8300  -13.9400; ...
       -13.5600  -13.4400  -13.0000  -12.9400; ...
       -12.8300  -12.6100  -12.0000  -12.0000; ...
       -11.8000  -11.6000  -11.0000  -11.0000; ...
      ...
     ];
exatidaoZcarga0=exatidaoZcarga0transp';


    exatidaoYcargaVtransp = [ ...
        24.7500   14.7800   28.6400    9.4100; ...
        23.3200   14.2200   28.6900    9.6300; ...
        21.1000   11.7300   27.5000    9.4000; ...
        19.9700   10.5700   25.2300    7.1700; ...
        16.3800    8.5400   22.0500    5.7500; ...
        10.5400    5.1100   18.3200    5.2000; ...
         8.9400    3.1800   16.0500    3.8700; ...
         8.1600    2.3500   14.8200    2.3000; ...
         4.2500    0.8600    7.6500    0.7100; ...
         1.2700    0.2300    0.2600   -0.2600; ...
         0.4700   -0.4600   -0.2400   -0.4900; ...
         0.0300   -0.8300    0.1900   -0.7700; ...
        -0.5400   -0.5200   -0.4400         0; ...
        -1.8000   -0.3800   -1.0300    0.2800; ...
        -2.0800   -0.8000   -1.4700   -0.4200; ...
        -2.4900   -1.3000   -1.8300   -0.7700; ...
        -2.2200   -3.1600   -7.5600   -0.6600; ...
        -1.8000   -4.9500  -14.8500   -1.3000; ...
        -2.3700   -5.5900  -15.6800   -3.2800; ...
        -5.6100   -9.5800  -20.1000   -7.4100; ...
        -6.9500  -11.0600  -22.6600   -9.2100; ...
        -7.6600  -13.1800  -26.2400   -9.4500; ...
       -12.3700  -16.4000  -26.7100  -10.2100; ...
       -13.4200  -17.0500  -27.1800  -10.4100; ...
       -13.0400  -16.7400  -27.1200  -10.1800; ...
      ...
     ];
 exatidaoYcargaV= exatidaoYcargaVtransp';


    exatidaoZcargaVtransp = [ ...
        -8.5000   -8.8000   -8.2000   -6.7000; ...
        -7.1100   -8.2200   -8.1100   -6.5600; ...
        -5.5000   -6.0000   -7.2800   -6.5600; ...
        -4.5600   -5.0600   -5.1700   -4.3900; ...
        -4.9000   -5.0500   -5.3000   -3.9000; ...
        -1.3900   -2.7200   -4.5000   -3.3300; ...
         0.8900   -0.1700   -1.7800   -2.0000; ...
         1.1100    0.6100   -0.1700    0.2800; ...
         1.7000    0.9500    0.3500    1.0000; ...
         2.0600    1.3300    0.5000    1.0000; ...
         2.7800    1.7200    0.6100    1.0600; ...
         2.6700    1.8900    1.0000    1.4400; ...
         3.3000    2.2000    1.0000    1.2000; ...
         3.3300    2.0600    1.0000    1.3900; ...
         3.2200    2.1100    1.0000    1.1100; ...
         3.0600    1.9400    1.0000    1.1700; ...
         2.7000    1.4000    0.7000    0.8500; ...
         2.7200    1.1700         0   -0.4400; ...
         1.5600    0.6100   -1.2200   -3.2800; ...
        -2.7200   -4.9400   -6.3300   -7.6100; ...
        -4.3000   -5.9000   -7.8000   -9.9000; ...
        -4.7800   -7.4400  -10.1700   -9.8900; ...
        -9.5000  -10.1700   -9.8900   -9.7200; ...
        -9.2800   -9.5600   -9.1700   -8.8900; ...
        -8.7000   -8.8000   -8.8000   -8.7000; ...
      ...
     ];
exatidaoZcargaV=exatidaoZcargaVtransp';


    exatidaoYcargaHtransp = [ ...
        33.8700   23.2700   20.1500   20.0300; ...
        29.1300   19.9800   16.8300   16.7500; ...
        26.6900   16.7200   14.6300   15.4400; ...
        25.1500   16.0300   13.1600   13.5300; ...
        19.8500   12.1400   11.3300   10.9800; ...
        14.6000    8.5300    9.5300    8.9400; ...
        13.9600    8.0400    9.0000    8.3100; ...
        12.4800    7.2900    8.4200    7.8900; ...
         7.7300    5.2300    5.6300    5.8600; ...
         3.6000    3.8100    3.2800    4.0900; ...
         2.1800    2.4900    1.8000    2.8300; ...
         1.4600    1.8100    0.9100    1.8900; ...
        -1.5800   -1.5800   -1.4200    0.3800; ...
        -1.2300   -0.6700   -1.5200   -1.2400; ...
        -1.2500   -0.8400   -1.5000   -2.3800; ...
        -1.3300   -1.0200   -1.6700   -2.6400; ...
        -1.6800   -1.7700   -2.5200   -3.7300; ...
        -2.2500   -2.4400   -3.9400   -5.6300; ...
        -3.5600   -4.0700   -6.4600   -9.7900; ...
        -7.6500   -8.3900  -10.4600  -12.6200; ...
        -8.0800  -11.0100  -12.4500  -14.3900; ...
       -10.4400  -14.2100  -14.1900  -15.4100; ...
       -12.3200  -14.8000  -14.7800  -15.9400; ...
       -13.0100  -15.3300  -15.1200  -16.1600; ...
       -12.6700  -15.0400  -14.8900  -16.0600; ...
      ...
     ];
exatidaoYcargaH=exatidaoYcargaHtransp';


    exatidaoZcargaHtransp = [ ...
       -17.4000  -16.3000  -15.3000  -14.9000; ...
       -13.5600  -13.8300  -13.1700  -12.3300; ...
       -11.8300  -11.5600  -11.6700  -11.7800; ...
       -10.8300  -11.5600  -10.7200  -10.7800; ...
       -11.0000  -11.4500  -10.8000  -10.3000; ...
       -10.9400  -11.8900  -10.8900  -10.4400; ...
       -10.9400  -11.8300  -10.8900  -10.3300; ...
        -9.4400  -11.5000  -11.0000  -10.6700; ...
        -7.1500   -9.6000  -10.2000  -10.1000; ...
        -5.6100   -8.1100   -9.6100   -9.5000; ...
        -2.5000   -5.8900   -7.0000   -7.5000; ...
        -1.0000   -4.3900   -5.2800   -6.2200; ...
         0.1000   -1.7000   -4.1000   -5.7000; ...
         0.1700   -1.4400   -1.6700   -4.2200; ...
         0.2200   -1.5000   -1.7200   -2.0000; ...
         0.2800   -1.5600   -1.8300   -2.0000; ...
         0.2500   -1.8000   -1.6500   -2.5500; ...
         0.0600   -2.2800   -3.1100   -4.3300; ...
        -2.2200   -4.1700   -6.4400   -9.5600; ...
        -7.8900  -10.0000  -11.1100  -12.2200; ...
        -7.9000  -10.8500  -11.9000  -12.3500; ...
       -10.3300  -13.0600  -11.6100  -11.8300; ...
       -11.2800  -12.2200  -11.2200  -11.0000; ...
       -10.8300  -11.5600  -10.4400  -10.0000; ...
       -10.0000  -10.5000  -10.0000   -9.5000; ...
      ...
     ];
exatidaoZcargaH=exatidaoZcargaHtransp';


    precisaoYcarga0transp = [ ...
         6.0200    1.8100    3.9000   11.0200; ...
         3.1900    6.0500    4.6400    5.9700; ...
         3.6800    5.9300    4.2900    6.3400; ...
         3.1400    4.9500    4.2700    6.7500; ...
         5.4700    8.3400    7.7900   10.2600; ...
         1.3600    1.6100    5.4700    2.9700; ...
         3.7700    3.8400    5.2100    1.7200; ...
         2.8000    4.8100    6.5800    3.0800; ...
         1.4700    3.3800    4.4100    1.7800; ...
         2.2200    3.9000    3.6400    2.4500; ...
         2.7200    4.0000    3.1700    1.1000; ...
         2.8800    3.9000    3.2700    2.1000; ...
         6.5500    9.6100    4.2300    7.1900; ...
         3.4900    4.8200    2.0100    6.1000; ...
         3.2800    4.3400    1.8900    6.9000; ...
         2.7100    3.7200    2.2000    6.8700; ...
         4.0400    3.5000    1.8100    4.5200; ...
         3.8300    3.3500    1.0300    2.0700; ...
         3.9900    3.1700    2.5300    2.3400; ...
         2.7700    2.5900    2.6900    2.3300; ...
         1.6800    1.4800    3.3500    4.7800; ...
         1.1900    2.9700    2.6800    2.6900; ...
         1.5700    3.0500    2.7500    2.8400; ...
         1.8400    2.9000    2.6000    2.2900; ...
         1.0600    1.1600    0.6500    1.5300; ...
      ...
     ];
precisaoYcarga0=precisaoYcarga0transp';


    precisaoZcarga0transp = [ ...
         0.6700    6.2400    2.0100    5.8500; ...
         1.0600    1.0900    0.5000    1.6400; ...
         1.0900    0.9000    1.0900    1.3000; ...
         0.5000    1.3900    1.0900    1.2500; ...
         0.4700    1.2200    1.0600    1.1800; ...
         5.4600    3.1800    1.7000    0.7900; ...
         9.8400    5.8300    5.0900    5.0900; ...
         8.8900   11.0500    7.7200    5.4500; ...
         5.2700    9.0600    8.2400    9.4400; ...
              0    3.2700    8.0100    9.0200; ...
         0.6600    1.5000    2.8100    5.9500; ...
              0    1.9000    1.9500    1.3000; ...
         1.4900    1.2600    2.7200    1.4100; ...
         0.7500    0.7500    0.7900    1.0900; ...
         0.7900    1.3000    1.3000    1.9800; ...
         1.5000    2.8100    5.0900    5.6600; ...
         6.3400    4.3000    1.6200    4.4300; ...
         6.1400    4.7100    3.5200    0.7900; ...
         4.9600    4.8500    0.7900    0.5000; ...
         1.3200    1.0600    0.7900    0.6600; ...
         1.7000    1.0500         0    1.4500; ...
         1.4600    0.7500    0.7500    0.9000; ...
         0.5000    0.9000         0    0.9000; ...
         0.7500    0.6600         0         0; ...
         1.3400    1.2500         0         0; ...
      ...
     ];
precisaoZcarga0=precisaoZcarga0transp';


    precisaoYcargaVtransp = [ ...
         7.1600    3.7800    5.4400    0.8300; ...
         7.2000    2.6900    3.8100    1.6000; ...
         8.7700    6.1200    3.3600    2.0700; ...
         8.9600    7.5900    6.7100    4.1300; ...
         7.5300    6.3800    6.8700    5.3400; ...
         7.5700    8.6400    8.7900    5.4300; ...
         8.2900    8.1600    7.6600    4.1800; ...
         8.4700    8.0500    7.5500    3.9400; ...
        13.5000    7.1700   25.5000    5.5100; ...
         3.6900    4.2300   12.6200    3.4700; ...
         3.3600    4.4800   12.7500    4.2200; ...
         3.0000    4.5100   15.5600    4.2700; ...
         4.5200    3.0400   11.0600    2.5100; ...
         3.7900    4.0000    3.0700    3.5800; ...
         3.9200    3.6700    3.3500    3.3100; ...
         3.6400    3.6600    3.2100    3.1600; ...
         1.4100    7.2800   24.1100    2.1600; ...
         2.9600    3.0700   13.2300    2.4800; ...
         2.7400    2.7700   10.3700    4.9200; ...
         3.6700    2.5900   10.3000    4.6500; ...
        10.2500    9.9600    8.2300    1.8200; ...
         5.2200    3.7600    2.5700    2.3600; ...
         3.0700    2.2700    2.9000    2.7300; ...
         2.9100    2.3300    2.9600    2.7000; ...
         2.3500    1.5800    1.7000    1.6600; ...
      ...
     ];
precisaoYcargaV=precisaoYcargaVtransp';


    precisaoZcargaVtransp = [ ...
         2.3700    1.7100    2.2700    0.8200; ...
         5.8500    2.8200    1.4600    1.5800; ...
         9.8900    8.2500    4.5700    1.1700; ...
        11.6300   10.9800    8.8400    7.0700; ...
        11.7900   11.8800   10.5200    9.8700; ...
         5.0700    7.8000   11.1200   10.4700; ...
         1.0000    1.6800    4.3200    7.1900; ...
         1.0000    0.6600    0.7500    1.8500; ...
         0.7700    0.4700    0.7200         0; ...
         0.5000    1.0600         0         0; ...
         0.7900    1.0900    0.6600    0.5000; ...
         1.0600    0.6600         0    0.5000; ...
         0.7700    0.7700         0    0.7700; ...
         0.7500    0.5000         0    0.6600; ...
         0.7900    0.6600         0    0.6600; ...
         1.1700    0.5000         0    0.7500; ...
         0.7700    0.6300    0.7700    0.7200; ...
         1.3200    1.0600         0    4.4700; ...
         1.1700    1.2500    4.2500   10.4500; ...
         6.8800    0.9000    3.9000    9.9400; ...
        15.3500   16.2000   10.9100    2.2100; ...
         7.9100    4.9400    1.0600    1.0000; ...
         2.9000    1.0600    0.6600    1.0900; ...
         1.7000    1.3900    1.0600    0.6600; ...
         0.8200    0.8200    0.8200    0.8200; ...
      ...
     ];
precisaoZcargaV=precisaoZcargaVtransp';


    precisaoYcargaHtransp = [ ...
         5.4500    2.7300    3.2400    1.6100; ...
         9.2900    8.6200    7.9600    8.8100; ...
         6.3200    2.9000    3.8300    5.0700; ...
         5.1100    2.7800    2.9400    1.8500; ...
        16.7300   11.7100    6.3100    7.4300; ...
         6.6900    2.4000    4.3700    5.0400; ...
         6.5900    2.1300    4.2300    4.7500; ...
         5.5400    1.5500    4.1100    4.3700; ...
         9.1300    1.4800    6.8700    6.4500; ...
         5.5900    1.6200    1.0200    2.2500; ...
         5.3600    2.4000    2.0300    1.8400; ...
         5.4900    1.5500    2.5600    2.7300; ...
         6.7600    4.2900    0.7500    1.6500; ...
         2.6500    2.4600    2.5000    1.2600; ...
         2.6800    2.5900    2.6800    2.7100; ...
         2.7000    2.3600    2.5600    2.3900; ...
         2.8000    3.5000    3.6800    3.2600; ...
         1.9500    2.7800    1.1800    1.1700; ...
         1.4800    2.5200    3.0100    3.5500; ...
         1.0400    3.3900    2.9000    0.6100; ...
        11.0200    2.5500    2.6100    4.5700; ...
         3.7100    1.9700    2.5100    3.7200; ...
         2.7100    2.3000    2.8400    3.6200; ...
         3.0600    2.7000    2.5900    3.4700; ...
         2.0900    1.2200    0.9300    3.2100; ...
      ...
     ];
precisaoYcargaH=precisaoYcargaHtransp';


    precisaoZcargaHtransp = [ ...
         2.4600    1.7100    1.3400    1.6400; ...
        10.9100   11.3500   10.0100   10.8400; ...
         5.4600    2.1800    6.1800    7.6900; ...
         0.7500    1.9000    1.5200    3.5300; ...
         1.0000    1.1100    0.7700    1.4500; ...
         0.5000    0.6600    0.6600    0.9000; ...
         1.7500    0.7500    0.6600    1.0600; ...
         5.6800    1.8400         0    1.0600; ...
        11.0000    5.6000    2.2600    0.6300; ...
         9.3600    8.2800    3.4200    1.9800; ...
         7.4200   11.0400    9.6300    7.2700; ...
         2.4900    8.6100   10.9700   10.1200; ...
         0.6300    0.7700    6.7000    9.8500; ...
         0.7500    0.5000    1.3000    7.9100; ...
         0.7900         0    0.7900         0; ...
         0.7900    0.5000    0.7500         0; ...
         0.7900    0.7700    1.2300    1.4900; ...
         1.3900    2.0000    4.1600    5.3000; ...
         6.0500    6.1800    8.8600    8.7300; ...
         4.4200    6.5000    7.1500    2.8200; ...
        16.4900    8.7500    3.6600    0.7200; ...
         6.2300    0.9000    0.6600    0.7500; ...
         1.3200    0.7900    0.7900         0; ...
         1.0600    0.5000    0.9000         0; ...
              0         0         0         0; ...
      ...
     ];
precisaoZcargaH=precisaoZcargaHtransp';




% % Plotagem dos dados experimentais
% % EXATIDÃO
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoYcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em X sem carga (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoXcarga0.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoZcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em Y sem carga (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoYcarga0.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoYcargaV,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em X carga vertical (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoXcargaV.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoZcargaV,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em Y carga vertical (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoYcargaV.svg'])
% 
% 
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoYcargaH,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em X carga horizontal (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoXcargaH.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoZcargaH,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de exatidão em Y carga horizontal (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['exatidaoYcargaH.svg'])
% 
% % PRECISÃO
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoYcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em X sem carga (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoXcarga0.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoZcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em Y sem carga (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoYcarga0.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoYcargaV,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em X carga vertical (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoXcargaV.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoZcargaV,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em Y carga vertical (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoYcargaV.svg'])
% 
% 
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoYcargaH,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em X carga horizontal (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoXcargaH.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,precisaoZcargaH,'b.','MarkerSize',10) 
% view(-64,32)
% title(['erro de precisão em Y carga horizontal (mm)'])
% zlabel('erro (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['precisaoYcargaH.svg'])
% 
% 





% % Deslocamentos via diferença de exatidão
% 
% % figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoYcargaH-exatidaoYcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['deslocamento y (exatidaoYcargaH-exatidaoYcarga0) (mm)'])
% zlabel('deslocamento (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['deslocamentoyexp.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,exatidaoZcargaV-exatidaoZcarga0,'b.','MarkerSize',10) 
% view(-64,32)
% title(['deslocamento z (exatidaoZcargaV-exatidaoZcarga0) (mm)'])
% zlabel('deslocamento (mm)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% saveas(gcf,['deslocamentozexp.svg'])










% Elementos da matriz compliância experimental
cexp11=(exatidaoYcargaH-exatidaoYcarga0)/3.2/1000;
cexp12=(exatidaoYcargaV-exatidaoYcarga0)/3.2/1000;
cexp21=(exatidaoZcargaH-exatidaoZcarga0)/3.2/1000;
cexp22=(exatidaoZcargaV-exatidaoZcarga0)/3.2/1000;


% Matriz de rigidez experimental
kexp11=zeros(4,25);
kexp12=zeros(4,25);
kexp21=zeros(4,25);
kexp22=zeros(4,25);
numcondexp=zeros(4,25);
for i=1:4
   for j=1:25 
       kexp=inv([cexp11(i,j) cexp12(i,j);
                 cexp21(i,j) cexp22(i,j)]);
        kexp11(i,j)=kexp(1,1);
        kexp12(i,j)=kexp(1,2);
        kexp21(i,j)=kexp(2,1);
        kexp22(i,j)=kexp(2,2);
        numcondexp(i,j)=cond(kexp);
end
end


% Rigidez em x via inversão da matriz compliância
figure
surf(XB*1000,YB*1000,deslx)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,kexp21,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in x(N/m) compliancia efetiva fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
zlim([0,1.35e4])
% axis equal
saveas(gcf,['dx diferenca exatidao compliancia efetiva fator',num2str((fator-1)*100),'%.svg'])
%  
% 
% Rigidez em y via inversão da matriz compliância
figure
surf(XB*1000,YB*1000,desly)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,kexp22,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in y(N/m) compliancia efetiva fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['dy diferenca exatidao compliancia efetiva fator',num2str((fator-1)*100),'%.svg'])












% n=5;
% levels=round(linspace(-0.0070,0.0070,25).*10^n)./10^n;
% subplot(1,2,1)
% surf(XB*1000,YB*1000,deslx)
% alpha 0.3
% hold on
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,stiffxexp,'b.','MarkerSize',10) 

% set(gca, 'YDir','reverse')
% title('Displacement in x (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dx(mm)')
% subplot(1,2,2)
% [~,h] = contour(XB*1000,YB*1000,deslx*1000,levels);
% [~,h] = contour(XB*1000,YB*1000,deslx);
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(stiffxexp,[],1))))
% 




% % set(gca, 'YDir','reverse')
% title(['Stiffness in x(N/m) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% axis equal
% saveas(gcf,['dx fator',num2str((fator-1)*100),'%.svg'])
%  


% % n=5;
% % levels=round(linspace(0.016,0.028,15).*10^n)./10^n;
% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,desly)
% % alpha 0.3
% % hold on
% % x=-240e-3:20e-3:240e-3;
% % y=410e-3:20e-3:470e-3;
% % [x1,y1]=meshgrid(x,y);
% % stem3(x1*1000,y1*1000,stiffyexp,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
% % title('Displacement in y (Concentrated Parameters)')
% % xlabel('x(mm)')
% % ylabel('y(mm)')
% % zlabel('dy(mm)')
% % subplot(1,2,2)
% % [~,h] = contour(XB*1000,YB*1000,desly*1000,levels);
% [~,h] = contour(XB*1000,YB*1000,desly);
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(stiffyexp,[],1))))
% 
% 
% title(['Stiffness in y(N/m) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% axis equal 
% saveas(gcf,['dy fator',num2str((fator-1)*100),'%.svg'])
% 
% 
% 
% % n=2;
% % omegaarred=ceil(omega.*10^n)./10^n;
% 
% 
% 
% % n=1;
% % levels=round(linspace(30,40,20).*10^n)./10^n;
% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,omega1)
% % alpha 0.3
% % hold on
% % x=-240e-3:20e-3:240e-3;
% % y=410e-3:20e-3:470e-3;
% % [x1,y1]=meshgrid(x,y);
% % stem3(x1*1000,y1*1000,freqyexp,'b.','MarkerSize',10) 
% 
% % set(gca, 'YDir','reverse')
% % title('First natural frequency (Concentrated Parameters)')
% % xlabel('x(m)')
% % ylabel('y(m)')
% % zlabel('\omega (Hz)')
% % subplot(1,2,2)
% % [~,h] = contour(XB*1000,YB*1000,omega1,levels);
% [~,h] = contour(XB*1000,YB*1000,omega1);
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(freqyexp,[],1))))
% 
% title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% axis equal 
% saveas(gcf,['freq 1st fator',num2str((fator-1)*100),'%.svg'])
% 
% hold on
% % plot([0 0 -180 180 -180 180 -110 110 -50 50 -185 185],[-500 -425 -500 -500 -405 -405 -445 -445 -475 -475 -515 -515],'xr')
% 
% 
% 
% 
% 
% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,autovecy)
% % set(gca, 'YDir','reverse')
% % title('Displacement in y (Concentrated Parameters)')
% % xlabel('x(mm)')
% % ylabel('y(mm)')
% % zlabel('dy(mm)')
% % subplot(1,2,2)
% [~,h] = contour(XB*1000,YB*1000,autovecx1,15);
% % set(gca, 'YDir','reverse')
% title(['Component x eigenvector without the passive chain 1st frequency fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 


% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,autovecy)
% % set(gca, 'YDir','reverse')
% % title('Displacement in y (Concentrated Parameters)')
% % xlabel('x(mm)')
% % ylabel('y(mm)')
% % zlabel('dy(mm)')
% % subplot(1,2,2)
% [~,h] = contour(XB*1000,YB*1000,autovecy1,15);
% % set(gca, 'YDir','reverse')
% title(['Component y eigenvector without the passive chain 1st frequency fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 


% figure
% quiver(XB*1000,YB*1000,autovecx1,autovecy1);
% % set(gca, 'YDir','reverse')
% title(['Eigenvector mapping of the first natural frequency(2gdl) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% saveas(gcf,['eigen1st(2gdl) mob2 fator',num2str((fator-1)*100),'%.bmp'])
% % axis equal 



% n=1;
% 
% % levels=round(linspace(30,40,20).*10^n)./10^n;
% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,omega2)
% % alpha 0.3
% % hold on
% % x=-240e-3:20e-3:240e-3;
% % y=410e-3:20e-3:470e-3;
% % [x1,y1]=meshgrid(x,y);
% % stem3(x1*1000,y1*1000,freqxexp,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
% % title('First natural frequency (Concentrated Parameters)')
% % xlabel('x(m)')
% % ylabel('y(m)')
% % zlabel('\omega (Hz)')
% % subplot(1,2,2)
% % [~,h] = contour(XB*1000,YB*1000,omega2,levels);
% [~,h] = contour(XB*1000,YB*1000,omega2);
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(freqxexp,[],1))))
% 
% title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% axis equal 
% saveas(gcf,['freq 2nd fator',num2str((fator-1)*100),'%.svg'])
% 
% hold on
% % plot([0 0 -180 180 -180 180 -110 110 -50 50 -185 185],[-500 -425 -500 -500 -405 -405 -445 -445 -475 -475 -515 -515],'xr')


% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,autovecy)
% % set(gca, 'YDir','reverse')
% % title('Displacement in y (Concentrated Parameters)')
% % xlabel('x(mm)')
% % ylabel('y(mm)')
% % zlabel('dy(mm)')
% % subplot(1,2,2)
% [~,h] = contour(XB*1000,YB*1000,autovecx2,15);
% % set(gca, 'YDir','reverse')
% title(['Component x eigenvector without the passive chain 2nd frequency fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 


% figure
% % subplot(1,2,1)
% % surf(XB*1000,YB*1000,autovecy)
% % set(gca, 'YDir','reverse')
% % title('Displacement in y (Concentrated Parameters)')
% % xlabel('x(mm)')
% % ylabel('y(mm)')
% % zlabel('dy(mm)')
% % subplot(1,2,2)
% [~,h] = contour(XB*1000,YB*1000,autovecy2,15);
% % set(gca, 'YDir','reverse')
% title(['Component y eigenvector without the passive chain 2nd frequency fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 


% figure
% quiver(XB*1000,YB*1000,autovecx2,autovecy2);
% % set(gca, 'YDir','reverse')
% title(['Eigenvector mapping of the second natural frequency(2gdl) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% saveas(gcf,['eigen2nd(2gdl) mob2 fator',num2str((fator-1)*100),'%.bmp'])
% % axis equal 



























% % Mapeamento com density plot
% 
% contourf(XB*1000,YB*1000,deslx);
% colorbar
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(stiffxexp,[],1))))
% 
% % set(gca, 'YDir','reverse')
% title(['Stiffness in x(N/m) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal
% saveas(gcf,['dx fator',num2str((fator-1)*100),'%.svg'])
%  
% 
% figure
% 
% contourf(XB*1000,YB*1000,desly);
% colorbar
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(stiffyexp,[],1))))
% 
% 
% title(['Stiffness in y(N/m) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 
% saveas(gcf,['dy fator',num2str((fator-1)*100),'%.svg'])
% 
% 
% figure
% 
% contourf(XB*1000,YB*1000,omega1);
% colorbar
% 
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(freqyexp,[],1))))
% 
% title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 
% saveas(gcf,['freq 1st fator',num2str((fator-1)*100),'%.svg'])
% 
% hold on
% 
% 
% figure
% 
% contourf(XB*1000,YB*1000,omega2);
% colorbar
% % set(gca, 'YDir','reverse')
% hold on
% x1=-240e-3:20e-3:240e-3;
% y1=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% plot(x1*1000,y1*1000,'rx')
% text(reshape(x1*1000,[],1),reshape(y1*1000,[],1),cellstr(num2str(reshape(freqxexp,[],1))))
% 
% title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
% grid on
% % axis equal 
% saveas(gcf,['freq 2nd fator',num2str((fator-1)*100),'%.svg'])
% 
% hold on















% Mapeamento em superfície
figure
surf(XB*1000,YB*1000,deslx)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,stiffxexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in x(N/m) fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal
saveas(gcf,['dx fator',num2str((fator-1)*100),'%.svg'])
 


% gráfico rigidez em x com zoom

surf(XB*1000,YB*1000,deslx)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,stiffxexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in x(N/m) fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
zlim([0,1.35e4])
% axis equal
saveas(gcf,['dx zoom fator',num2str((fator-1)*100),'%.svg'])



figure
surf(XB*1000,YB*1000,desly)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,stiffyexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in y(N/m) fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['dy fator',num2str((fator-1)*100),'%.svg'])


% gráfico rigidez em y com zoom

figure
surf(XB*1000,YB*1000,desly)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,stiffyexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in y(N/m) fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
zlim([0,6000])
% axis equal 
saveas(gcf,['dy zoom fator',num2str((fator-1)*100),'%.svg'])



figure
surf(XB*1000,YB*1000,omega1)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,freqyexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
zlabel('frequência (Hz)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['freq 1st fator',num2str((fator-1)*100),'%.svg'])


figure
surf(XB*1000,YB*1000,omega2)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
stem3(x1*1000,y1*1000,freqxexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
zlabel('frequência (Hz)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['freq 2nd fator',num2str((fator-1)*100),'%.svg'])











% Mapeamento em superfície considerando a precisão como fonte de folgas a
% serem eliminadas no cálculo da rigidez
figure
surf(XB*1000,YB*1000,deslx)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
% rigidez com a precisão já descontada como folga
clearancex=0.001*abs(((precisaoYcargaH).^2+(precisaoZcargaH).^2).^(1/2)-((precisaoYcarga0).^2+(precisaoZcarga0).^2).^(1/2));
% clearancex=0.0000001;
stiffxexpnoclearances=3.2./abs(3.2./(stiffxexp)-clearancex);

stem3(x1*1000,y1*1000,stiffxexpnoclearances,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in x(N/m) prec folga fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% zlim([0,1.35e4])
% axis equal
saveas(gcf,['dx prec folga fator',num2str((fator-1)*100),'%.svg'])
 


figure
surf(XB*1000,YB*1000,desly)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
% rigidez com a precisão já descontada como folga
clearancey=0.001*abs(((precisaoYcargaV).^2+(precisaoZcargaV).^2).^(1/2)-((precisaoYcarga0).^2+(precisaoZcarga0).^2).^(1/2));
% clearancey=0.00001;
stiffyexpnoclearances=3.2./abs(3.2./(stiffyexp)-clearancey);

stem3(x1*1000,y1*1000,stiffyexpnoclearances,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in y(N/m) prec folga fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% zlim([0,6000])
% axis equal 
saveas(gcf,['dy prec folga fator',num2str((fator-1)*100),'%.svg'])



% figure
% surf(XB*1000,YB*1000,omega1)
% alpha 0.3
% hold on
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqyexp,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
% view(-64,32)
% title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% % axis equal 
% saveas(gcf,['freq 1st fator',num2str((fator-1)*100),'%.svg'])
% 
% 
% figure
% surf(XB*1000,YB*1000,omega2)
% alpha 0.3
% hold on
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqxexp,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
% view(-64,32)
% title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% % axis equal 
% saveas(gcf,['freq 2nd fator',num2str((fator-1)*100),'%.svg'])
 







% % Mapeamento em superfície considerando exatidão como ponto de partida para
% % o cálculo da rigidez, sendo possível determinar os quatro elementos da
% % matriz de rigidez do efetuador
figure
surf(XB*1000,YB*1000,deslx)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % % Eliminando os pontos de transiçao do gráfico
% % x(:,[1 5 9 13 17 21 25])=[];
% 
y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);
% 
% % stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcarga0));
% % % Eliminando os pontos de transiçao do gráfico
% % stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];
% 
% % stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH));
% % % Eliminando os pontos de transiçao do gráfico
% % stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];
% 
stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH-exatidaoYcarga0));
% % % Eliminando os pontos de transiçao do gráfico
% % stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];
% 
% stiffxexpbaseexatidao=3.2./(0+abs(0.001*(exatidaoYcargaH-exatidaoYcarga0)));
% % % Eliminando os pontos de transiçao do gráfico
% % stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];
% 
% 
% % stem3(x1*1000,y1*1000,exatidaoYcarga0,'b.','MarkerSize',10) 
% % stem3(x1*1000,y1*1000,exatidaoYcargaH,'b.','MarkerSize',10) 
stem3(x1*1000,y1*1000,stiffxexpbaseexatidao,'b.','MarkerSize',10) 
% % set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in x(N/m) compliancia fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
zlim([0,1.35e4])
% axis equal
saveas(gcf,['dx diferenca exatidao compliancia fator',num2str((fator-1)*100),'%.svg'])
%  
% 
% 
figure
surf(XB*1000,YB*1000,desly)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % Eliminando os pontos de transiçao do gráfico
% x(:,[1 5 9 13 17 21 25])=[];

y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV-exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=3.2./(-0.0015+abs(0.001*(exatidaoZcargaV-exatidaoZcarga0)));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stem3(x1*1000,y1*1000,exatidaoZcarga0,'b.','MarkerSize',10)
% stem3(x1*1000,y1*1000,exatidaoZcargaV,'b.','MarkerSize',10)

stem3(x1*1000,y1*1000,stiffyexpbaseexatidao,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['Stiffness in y(N/m) compliancia fator',num2str((fator-1)*100),'%'])
zlabel('rigidez (N/m)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['dy diferenca exatidao compliancia fator',num2str((fator-1)*100),'%.svg'])
% 
% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqyexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
% view(-64,32)
% title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% zlabel('frequência (Hz)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% axis equal 
% saveas(gcf,['freq 1st fator EXP',num2str((fator-1)*100),'%.svg'])


% figure
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqxexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
% view(-64,32)
% title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% zlabel('frequência (Hz)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% axis equal 
% saveas(gcf,['freq 2nd fator EXP',num2str((fator-1)*100),'%.svg'])




% figure
% surf(XB*1000,YB*1000,omega1)
% alpha 0.3
% hold on
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqyexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
% view(-64,32)
% title(['First natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% zlabel('frequência (Hz)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% axis equal 
% saveas(gcf,['freq 1st fator',num2str((fator-1)*100),'%.svg'])


% figure
% surf(XB*1000,YB*1000,omega2)
% alpha 0.3
% hold on
% x=-240e-3:20e-3:240e-3;
% y=410e-3:20e-3:470e-3;
% [x1,y1]=meshgrid(x,y);
% stem3(x1*1000,y1*1000,freqxexp,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
% view(-64,32)
% title(['Second natural frequency(Hz) fator',num2str((fator-1)*100),'%'])
% zlabel('frequência (Hz)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% axis equal 
% saveas(gcf,['freq 2nd fator',num2str((fator-1)*100),'%.svg'])

% 
% 
% 
% 
% 
% 
% % Número de condicionamento da matriz de rigidez 2x2
% figure
% surf(XB*1000,YB*1000,numcond)
% alpha 0.3
% % hold on
% % x=-240e-3:20e-3:240e-3;
% % y=410e-3:20e-3:470e-3;
% % [x1,y1]=meshgrid(x,y);
% % 
% % % número de condicionamento da matriz de rigidez experimental obtida a
% % % partir da matriz compliância
% % stem3(x1*1000,y1*1000,numcondexp,'b.','MarkerSize',10)
% 
% view(-64,32)
% title(['Conditioning number 2x2 stiffness matrix'])
% xlabel('x(mm)')
% ylabel('y(mm)')
% grid on
% xticks(-240:80:240)
% yticks(410:20:470)
% % axis equal 
% saveas(gcf,['numcond','%.svg'])





% Mapeamento em superfície considerando exatidão como ponto de partida para
% o cálculo da compliância, sendo possível determinar os quatro elementos da
% matriz de compliância do efetuador
figure
surf(XB*1000,YB*1000,deslx*3.2*1000)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % Eliminando os pontos de transiçao do gráfico
% x(:,[1 5 9 13 17 21 25])=[];

y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH-exatidaoYcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

compxexpbaseexatidao=abs(exatidaoYcargaH-exatidaoYcarga0);
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];


% stem3(x1*1000,y1*1000,exatidaoYcarga0,'b.','MarkerSize',10) 
% stem3(x1*1000,y1*1000,exatidaoYcargaH,'b.','MarkerSize',10) 
stem3(x1*1000,y1*1000,compxexpbaseexatidao,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['flexibilidade em y(mm/3,2N) abs(exatidaoYcargaH-exatidaoYcarga0)'])
xlabel('x(mm)')
ylabel('y(mm)')
zlabel('deslocamento (mm)')

grid on
xticks(-240:80:240)
yticks(410:20:470)
% zlim([0,1.35e4])
% axis equal
saveas(gcf,['compy.svg'])
 


figure
surf(XB*1000,YB*1000,desly*3.2*1000)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % Eliminando os pontos de transiçao do gráfico
% x(:,[1 5 9 13 17 21 25])=[];

y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV-exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

compyexpbaseexatidao=abs(exatidaoZcargaV-exatidaoZcarga0);
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stem3(x1*1000,y1*1000,exatidaoZcarga0,'b.','MarkerSize',10)
% stem3(x1*1000,y1*1000,exatidaoZcargaV,'b.','MarkerSize',10)

stem3(x1*1000,y1*1000,compyexpbaseexatidao,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['flexibilidade em z(mm/3,2N) abs(exatidaoZcargaV-exatidaoZcarga0)'])
xlabel('x(mm)')
ylabel('y(mm)')
zlabel('deslocamento (mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['compz.svg'])





% Mapeamento em superfície apenas da compliância teórica para a a carga de 
% 3,2N tanto na vertical como na horizontal
figure
surf(XB*1000,YB*1000,deslx*3.2*1000)
alpha 0.3
hold on
view(-64,32)
title(['flexibilidade em y(mm/3,2N) apenas teórico'])
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% zlim([0,1.35e4])
% axis equal
saveas(gcf,['compyteo.svg'])
 


figure
surf(XB*1000,YB*1000,desly*3.2*1000)
alpha 0.3
hold on
view(-64,32)
title(['flexibilidade em z(mm/3,2N) apenas teorico'])
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['compzteo.svg'])








% Mapeamento em superfície considerando exatidão como ponto de partida para
% o cálculo da compliância, supondo os erro de precisão como fonte de
% folgas
figure
surf(XB*1000,YB*1000,deslx*3.2*1000)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % Eliminando os pontos de transiçao do gráfico
% x(:,[1 5 9 13 17 21 25])=[];

y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffxexpbaseexatidao=1./(0.001/3.2*abs(exatidaoYcargaH-exatidaoYcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

compxexpbaseexatidao=abs(abs(exatidaoYcargaH-exatidaoYcarga0)-3/3*abs(precisaoYcargaH-precisaoYcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffxexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];


% stem3(x1*1000,y1*1000,exatidaoYcarga0,'b.','MarkerSize',10) 
% stem3(x1*1000,y1*1000,exatidaoYcargaH,'b.','MarkerSize',10) 
stem3(x1*1000,y1*1000,compxexpbaseexatidao,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['flexibilidade em y(mm/3,2N) com precisão como folga'])
zlabel('deslocamento (mm)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% zlim([0,1.35e4])
% axis equal
saveas(gcf,['compyfolgasprec.svg'])
 


figure
surf(XB*1000,YB*1000,desly*3.2*1000)
alpha 0.3
hold on
x=-240e-3:20e-3:240e-3;
% % Eliminando os pontos de transiçao do gráfico
% x(:,[1 5 9 13 17 21 25])=[];

y=410e-3:20e-3:470e-3;
[x1,y1]=meshgrid(x,y);

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stiffyexpbaseexatidao=1./(0.001/3.2*abs(exatidaoZcargaV-exatidaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

compyexpbaseexatidao=abs(abs(exatidaoZcargaV-exatidaoZcarga0)-3/3*abs(precisaoZcargaV-precisaoZcarga0));
% % Eliminando os pontos de transiçao do gráfico
% stiffyexpbaseexatidao(:,[1 5 9 13 17 21 25])=[];

% stem3(x1*1000,y1*1000,exatidaoZcarga0,'b.','MarkerSize',10)
% stem3(x1*1000,y1*1000,exatidaoZcargaV,'b.','MarkerSize',10)

stem3(x1*1000,y1*1000,compyexpbaseexatidao,'b.','MarkerSize',10) 
% set(gca, 'YDir','reverse')
view(-64,32)
title(['flexibilidade em z(mm/3,2N) com precisão como folga'])
zlabel('deslocamento (mm)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
xticks(-240:80:240)
yticks(410:20:470)
% axis equal 
saveas(gcf,['compzfolgaprec.svg'])


