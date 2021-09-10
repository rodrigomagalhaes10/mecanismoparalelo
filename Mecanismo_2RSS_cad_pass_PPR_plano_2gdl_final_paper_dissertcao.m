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
autovecx=zeros(size(XB));
autovecy=zeros(size(XB));

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

% Nova nomenclatura Equações originais
% L0^2+L1^2-L2^2-2*L0*L5+L5^2-2*L0*xB+2*L5*xB+xB^2+yB^2+2*L1*(L0-L5-xB)*cosq1-2*L1*yB*sinq1
% L0^2+L3^2-L4^2-2*L0*L5+L5^2+2*L0*xB-2*L5*xB+xB^2+yB^2+2*L3*(L0-L5+xB)*cosq2-2*L3*yB*sinq2
G=-2*L1*yB;
H=2*L1*(L0-L5-xB);
I=L0^2+L1^2-L2^2-2*L0*L5+L5^2-2*L0*xB+2*L5*xB+xB^2+yB^2;
J=-2*L3*yB;
K=2*L3*(L0-L5+xB);
L=L0^2+L3^2-L4^2-2*L0*L5+L5^2+2*L0*xB-2*L5*xB+xB^2+yB^2;


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
% q11=atand((2*u1)/(1-u1^2))
q11=atand((sinq11)/(cosq11));

cosq12=(1-u2^2)/(1+u2^2);
sinq12=(2*u2)/(1+u2^2);
% q12=atand((2*u2)/(1-u2^2))
q12=atand((sinq12)/(cosq12));

% para o lado direito
cosq21=(1-v1^2)/(1+v1^2);
sinq21=(2*v1)/(1+v1^2);
% q21=atand((2*v1)/(1-v1^2))
q21=atand((sinq21)/(cosq21));

cosq22=(1-v2^2)/(1+v2^2);
sinq22=(2*v2)/(1+v2^2);
% q22=atand((2*v2)/(1-v2^2))
q22=atand((sinq22)/(cosq22));

% Solução escolhida para a cinemática inversa
sinq1=sinq12;
cosq1=cosq12;
sinq2=sinq22;
cosq2=cosq22;

% pontos para plotagem do grafico

% pontos O1 e O2
O1=[L0 0];
O2=[-L0 0];

% pontos B1 e B2
B1=[xB+L5 yB];
B2=[xB-L5 yB];

% ponto escolhido (configurações das cadeia periféricas para fora)
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

Jx=[2*(-L0+L5+xB-L1*cosq1) 2*(yB-L1*sinq1);
    2*(L0-L5+xB+L3*cosq2)  2*(yB-L3*sinq2)];

% % Matriz Jq
% % J1=(Q3(i,j)+L3-L0-L1*cosq1)*(L1*sinq1)+(Q4(i,j)-L1*sinq1)*(L1*cosq1);
% % J2=(Q3(i,j)-L3+L0+L1*cosq2)*(L1*sinq21)+(Q4(i,j)-L1*sinq2)*(L1*cosq2);
% % Antiga
% J1=G*cosq1-H*sinq1;
% J2=J*cosq2-K*sinq2;
% J3=2*(L1*cosq1+L0-xB-L3)*cosq1+2*(L1*sinq1-yB)*sinq1;
% J4=2*(L1*cosq2-L0-xB+L3)*cosq2+2*(L1*sinq2-yB)*sinq2;
% Jq=-1*[J1 0 J3 0 -2*L2 0;0 J2 0 J4 0 -2*L2];

Jq=[-2*L1*(-L0+L5+xB-L1*cosq1)*sinq1+2*L1*cosq1*(yB-L1*sinq1) 0;
    0                                                         2*L3*(L0-L5+xB+L3*cosq2)*sinq2+2*L3*cosq2*(yB-L3*sinq2)];

Jqext=[-2*L1*(-L0+L5+xB-L1*cosq1)*sinq1+2*L1*cosq1*(yB-L1*sinq1) 0                                                       2*cosq1*(-L0+L5+xB-L1*cosq1)+2*sinq1*(yB-L1*sinq1)  0                                                 2*L2 0;
        0                                                        2*L3*(L0-L5+xB+L3*cosq2)*sinq2+2*L3*cosq2*(yB-L3*sinq2) 0                                                  -2*cosq2*(L0-L5+xB+L3*cosq2)+2*sinq2*(yB-L3*sinq2) 0    2*L4];
% Matriz de rigidez do espaço dos atuadores
Kq=[k1 0 0 0 0 0;
    0 k2 0 0 0 0;
    0 0 k3 0 0 0;
    0 0 0 k4 0 0;
    0 0 0 0 k5 0;
    0 0 0 0 0 k6];


% Matriz de compliância para evitar a inversão de uma matriz retangular
Comp=inv(Jx)*Jqext*inv(Kq)*transpose(Jqext)*transpose(inv(Jx));

% Deslocamento
delta=Comp*[0;1];

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

Cx1=[-r1*sinq1 0;
      r1*cosq1 0;
      1        0]*J;
 
Cx3=[0 r3*sinq2;
     0 r3*cosq2;
     0 -1]*J;

% Cx2=[11 12;
%      21 22;
%      31 32];
 
Cx2=[r2/L2-((L2-r2)*(-L0+L5+xB-L1*cosq1)*sinq1)/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))  (-(((L2-r2)*sinq1*(yB-L1*sinq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))));
     ((L2-r2)*cosq1*(-L0+L5+xB-L1*cosq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1))       r2/L2+((L2-r2)*cosq1*(yB-L1*sinq1))/(L2*(yB*cosq1+(L0-L5-xB)*sinq1));
     ((-1*(L0-L5-xB+L1*cosq1)*sinq1)/(yB*cosq1+(L0-L5-xB)*sinq1)+1)/(L1*sinq1-yB)   (((sinq1*(yB-L1*sinq1))/((L1*sinq1-yB)*(yB*cosq1+(L0-L5-xB)*sinq1))))];

Cx4=[r4/L4+((L4-r4)*(L0-L5+xB+L3*cosq2)*sinq2)/(L4*(yB*cosq2+(L0-L5+xB)*sinq2))  ((L4-r4)*sinq2*(yB-L3*sinq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2));
     ((L4-r4)*cosq2*(L0-L5+xB+L3*cosq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2))      r4/L4+((L4-r4)*cosq2*(yB-L3*sinq2))/(L4*(yB*cosq2+(L0-L5+xB)*sinq2));
     (((-L0+L5-xB-L3*cosq2)*sinq2)/(yB*cosq2+(L0-L5+xB)*sinq2)+1)/(L3*sinq2-yB)  (-((sinq2*(yB-L3*sinq2))/((L3*sinq2-yB)*(yB*cosq2+(L0-L5+xB)*sinq2))))];
 
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

Mx=Cx1'*M1*Cx1+Cx2'*M2*Cx2+Cx3'*M3*Cx3+Cx4'*M4*Cx4+[1 0;0 1;0 0]'*[0.65 0 0;0 0.65 0;0 0 1]*[1 0;0 1;0 0];
% Mx=Cx1'*M1*Cx1+Cx2'*M2*Cx2+Cx3'*M3*Cx3+Cx4'*M4*Cx4;
 
 
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
autovecx(k,p)=Vordenado(1,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2);
autovecy(k,p)=Vordenado(2,1)/sqrt(Vordenado(1,1)^2+Vordenado(2,1)^2);

% autovetores da segunda frequencia
autovecx2(k,p)=Vordenado(1,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2);
autovecy2(k,p)=Vordenado(2,2)/sqrt(Vordenado(1,2)^2+Vordenado(2,2)^2);



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
title('Displacement in x(mm) without the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'dx(2gdl) mob2.bmp')
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
title('Displacement in y(mm) without the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'dy(2gdl) mob2.bmp')
% axis equal 


% n=2;
% omegaarred=ceil(omega.*10^n)./10^n;



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
title('First natural frequency(Hz) without the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
saveas(gcf,'freq(2gdl) mob2.bmp')
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
[~,h] = contour(XB*1000,YB*1000,autovecx,15);
set(gca, 'YDir','reverse')
title('Component x eigenvector without the passive chain')
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
[~,h] = contour(XB*1000,YB*1000,autovecy,15);
set(gca, 'YDir','reverse')
title('Component y eigenvector without the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
% axis equal 


figure
quiver(XB*1000,YB*1000,autovecx,autovecy);
set(gca, 'YDir','reverse')
% title('Eigenvector mapping of the first natural frequency(2gdl)')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
% axis equal 

% Segunda frequencia natural

n=1;
levels=round(linspace(40,55,20).*10^n)./10^n;
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
% title('Second natural frequency(Hz) without the passive chain')
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
% quiver(XB*1000,YB*1000,abs(autovecx2).*sign(XB),+abs(autovecy2).*sign(YB));
quiver(XB*1000,YB*1000,autovecx2,autovecy2);

set(gca, 'YDir','reverse')
% title('Eigenvector mapping of the second natural frequency with the passive chain')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
% axis equal 
