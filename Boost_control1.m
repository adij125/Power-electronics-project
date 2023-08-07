% State-space average model of a Boost Converter
% u=[delta];  x=[i_L; v_C]; 
% Component Parameters 
L=3.9e-3;  r_L=31e-3;  C=100e-6;  r_C=0.091;  R=251;

% Circuit Conditions 
V_I=100;Delta=0.75; 

% First form models of the two states 
A_on = [-r_L/L,  0;  0,  -1/C/(R+r_C)]; 
B_on = [1/L;  0];
C_on = [0,  R/(R+r_C)]; 

A_off = [-r_L/L-r_C*R/L/(R+r_C),  -R/L/(R+r_C);  R/C/(R+r_C),  -1/C/(R+r_C)]; 
B_off = [1/L;  0];
C_off = [r_C*R/(R+r_C),  R/(R+r_C)]; 

% Average the two models 
% Operating point model 
A = Delta*A_on + (1-Delta)*A_off; 
B = Delta*B_on + (1-Delta)*B_off; 
C = Delta*C_on + (1-Delta)*C_off; 

U=V_I;
X = -inv(A)*B*U; 
Y = C*X; 

%Small-signal model 
E = (A_on-A_off)*X + (B_on-B_off)*U;
F = (C_on-C_off)*X; 
boost_ss = ss(A,E,C,F); 

help=tf(boost_ss)
boost_ssCL = feedback(0.004* boost_ss, 1)
figure
rlocus(help)

% Find unit step-response 
T = 0:0.001:0.5;
[y,T,x] = step(boost_ssCL,T); 

% Scale step-response and add to operating point
iL = 0.05*x(:,1)+X(1)*ones(length(T),1);
y = 0.05*y+Y*ones(length(T),1);

T = [0; 0.15; T+0.15];
iL = [X(1); X(1); iL]; 
y = [Y; Y; y]; 

figure;

subplot(3,1,1)
plot([0, 0.15, 0.1501, 0.2], [Delta, Delta, Delta+0.5, Delta+0.5])
axis([0.145, 0.195, 0.45, 0.6])

subplot(3,1,2)
plot(T,iL)


subplot(3,1,3)
plot(T,y)

z = -1000;
p_plus_i = tf([1,-z],[1,0]);
piboost = series(p_plus_i, boost_ss)
k=1.5e-6;
piboost_closed1 = feedback(k*piboost, 1);
piboost_closed2 = feedback(0.0010*piboost, 1)
figure
rlocus(piboost);
axis([-1500, 1000, -400, 400])

% Find unit step-response 
T = 0:0.001:0.5;
[y,T,x] = step(piboost_closed2,T); 

% Scale step-response and add to operating point
iL = 0.05*x(:,1)+X(1)*ones(length(T),1);
y = 0.05*y+Y*ones(length(T),1);

T = [0; 0.15; T+0.15];
iL = [X(1); X(1); iL]; 
y = [Y; Y; y]; 

figure;

subplot(3,1,1)
plot([0, 0.15, 0.1501, 0.2], [0.75, 0.75, 0.80, 0.80])
axis([0.145, 0.195, 0.45, 0.6])

subplot(3,1,2)
plot(T,iL)


subplot(3,1,3)
plot(T,y)







