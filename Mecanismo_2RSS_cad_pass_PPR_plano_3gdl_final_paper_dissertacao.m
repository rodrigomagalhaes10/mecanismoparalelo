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

% dados da geometria
L0=(260e-3-160e-3/2);
L1=350e-3;
L2=422e-3;
L3=L1;
L4=L2;

% L5=0.000000000010;
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
Iz1=4831958.5e-9;
 
m2=2700*(pi/4*(14.3e-3)^2)*422e-3;
% m2=0;
Iz2=m2*L2^2/12;

m3=275.166e-3;
Iz3=4831958.5e-9;

m4=2700*(pi/4*(14.3e-3)^2)*422e-3;
% m4=0;
Iz4=m4*L4^2/12;

m5=0.65;
Iz5=m5*L5^2/12;
% Iz5=0;

m6=143e-3;
% m6=0;
Iz6=m6*L6^2/12;

m7=134e-3;
% m7=0;
Iz7=m7*L7^2/12;

m8=134e-3;
% m8=0;
Iz8=m8*L8^2/12;

% Distâncias dos centros de massa dos elos (1), (2), (4) e (5) (usualmente
% metade do comprimento do elo
r1=222.7e-3;
r2=L2/2;
r3=222.7e-3;
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
k1=3*Emat*Imat/L1;
k2=3*Emat*Imat/L3;
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
x=-200e-3:20e-3:200e-3;
y=400e-3:20e-3:520e-3;
[XB,YB]=meshgrid(x,y);

% Inicialização das matrizes com o tamanho correto das características a 
% serem mapeadas 
deslx=zeros(size(XB));
desly=zeros(size(XB));
omega1=zeros(size(XB));
autovecx1=zeros(size(XB));
autovecy1=zeros(size(XB));

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
k7=(L7/(4*L6^2*sinfi^2*A7*Emat)+L8/(4*L6^2*sinfi^2*A8*Emat)+L6^5/(6*Emat*I6))^(-1);

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

% Deslocamento
delta=Comp*[0;1;0];

% dados para plotagem da superfície
deslx(k,p)=delta(1,1);
desly(k,p)=delta(2,1);

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

% Cálculo da primeira frequência natural a partir das matrizes de massa e 
% de rigidez e ordenando em ordem crescente as frequências
[V,W] = eig(Mx\Kx);
[d,ind] = sort(diag(W));
Wordenado=W(ind,ind);
Vordenado=V(:,ind);
omega1(k,p)=(sqrt(Wordenado(1,1)))/(2*pi);

% Segunda frequencia
omega2(k,p)=(sqrt(Wordenado(2,2)))/(2*pi);



% Autovetores 
% autovecx(k,p)=Vordenado(1,1);
% autovecy(k,p)=Vordenado(2,1);
autovecx1(k,p)=Vordenado(1,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2+Vordenado(3,1)^2);
autovecy1(k,p)=Vordenado(2,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2+Vordenado(3,1)^2);



% autovetores da segunda frequencia
autovecx2(k,p)=Vordenado(1,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2+Vordenado(3,2)^2);
autovecy2(k,p)=Vordenado(2,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2+Vordenado(3,2)^2);



% [W] = eig(matmassa\Ktotal);
% omega(k,p)=(sqrt(min(W)))/(2*pi);
end
end

n=5;
levels=round(linspace(-0.0070,0.0070,25).*10^n)./10^n;

% subplot(1,2,1)
% surf(XB*1000,YB*1000,dx)
% set(gca, 'YDir','reverse')
% title('Displacement in x (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dx(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,deslx*1000,levels);
set(gca, 'YDir','reverse')
title('Displacement in x(mm) with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'dx(3gdl) mob3.bmp')
% axis equal 


n=5;
levels=round(linspace(0.016,0.028,15).*10^n)./10^n;
figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,dy)
% set(gca, 'YDir','reverse')
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,desly*1000,levels);
set(gca, 'YDir','reverse')
title('Displacement in y(mm) with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'dy(3gdl) mob3.bmp')
% axis equal 



n=1;
levels=round(linspace(30,40,20).*10^n)./10^n;
figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,omega)
% set(gca, 'YDir','reverse')
% title('First natural frequency (Concentrated Parameters)')
% xlabel('x(m)')
% ylabel('y(m)')
% zlabel('\omega (Hz)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,omega1,levels);
set(gca, 'YDir','reverse')
title('First natural frequency(Hz) with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'freq(3gdl) mob3.bmp')
% axis equal 
hold on
% plot([0 0 -180 180 -180 180 -110 110 -50 50 -185 185],[-500 -425 -500 -500 -405 -405 -445 -445 -475 -475 -515 -515],'xr')


figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,autovecy)
% set(gca, 'YDir','reverse')
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,autovecx1,15);
set(gca, 'YDir','reverse')
title('Component x eigenvector with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
% axis equal 


figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,autovecy)
% set(gca, 'YDir','reverse')
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,autovecy1,15);
set(gca, 'YDir','reverse')
title('Component y eigenvector with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
% axis equal 


figure
% quiver(XB*1000,YB*1000,autovecx1,autovecy1);
quiver(XB*1000,YB*1000,abs(autovecx1).*sign(XB),+abs(autovecy1));
set(gca, 'YDir','reverse')
% title('Eigenvector mapping of the first natural frequency with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
axis equal 



% Segunda frequencia natural

n=1;
levels=round(linspace(40,50,20).*10^n)./10^n;
figure

% subplot(1,2,1)
% surf(XB*1000,YB*1000,omega)
% set(gca, 'YDir','reverse')
% title('First natural frequency (Concentrated Parameters)')
% xlabel('x(m)')
% ylabel('y(m)')
% zlabel('\omega (Hz)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,omega2,levels);
set(gca, 'YDir','reverse')
% title('Second natural frequency(Hz) with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'freq2(3gdl) mob3.bmp')
% axis equal 
hold on
% plot([0 0 -180 180 -180 180 -110 110 -50 50 -185 185],[-500 -425 -500 -500 -405 -405 -445 -445 -475 -475 -515 -515],'xr')


figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,autovecy)
% set(gca, 'YDir','reverse')
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,autovecx2,15);
set(gca, 'YDir','reverse')
title('Component x eigenvector 2 with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
% axis equal 


figure
% subplot(1,2,1)
% surf(XB*1000,YB*1000,autovecy)
% set(gca, 'YDir','reverse')
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(XB*1000,YB*1000,autovecy2,15);
set(gca, 'YDir','reverse')
title('Component y eigenvector 2 with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
% axis equal 


figure
quiver(XB*1000,YB*1000,abs(autovecx2),-abs(autovecy2).*sign(XB));
% quiver(XB*1000,YB*1000,autovecx2,autovecy2);

set(gca, 'YDir','reverse')
% title('Eigenvector mapping of the second natural frequency with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
axis equal 



