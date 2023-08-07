% State-space average model of a Flyback Converter
% U=[V_in];  x=[big phi; v_C]; 
% u=[delta] 
% Component Parameters 
L1=0.01255; 
L2=8.5418438e-7; 
Lratio=sqrt(L1/L2);%(value for L2)
r_L=31e-3;  C=100e-6;  r_C=0.091;  R=0.09076;

% L = (big phi)/ I1

% Circuit Conditions 
V_I=400;Delta=0.5; 

% First form models of the two states 
A_on = [0,  -r_L/L1;  0,  -1/(C*(R + r_C))]; 
B_on = [1/L1;  0];
C_on = [0,  R/(R+r_C)]; 

A_off = [(r_L/L2-r_C*R/(L2*(R+r_C))),  -R/(Lratio*L2*((R+r_C)));  R*Lratio/(C*(R+r_C)),  (1)/(C*(r_C+R))]; 
B_off = 0;
C_off = [Lratio*r_C*R/(R+r_C),  R/(R+r_C)]; 

% Average the two models 
% Operating point model 
A = Delta*A_on + (1-Delta)*A_off; 
B = Delta*B_on + (1-Delta)*B_off; 
C = Delta*C_on + (1-Delta)*C_off; 

U=V_I;
X = -inv(A)*B*U;
Y = C*X

%Small-signal model 
% E = (A_on-A_off)*X + (B_on-B_off)*U;
% F = (C_on-C_off)*X; 
% control=tf([0,0.01444,514],[0,1,0])
% boost_ss = ss(A,E,C,F);
% rlocus(boost_ss)
% close all

% help=tf((1/2000)*boost_ss)
% boost_ssCL = feedback(boost_ss, 1)
% figure
% rlocus(help)
% figure
% bode(help)

% Find unit step-response 
% T = 0:0.001:0.05;
% [y,T,x] = step(boost_ss,T); 

% Scale step-response and add to operating point
% iL = 0.05*x(:,1)+X(1)*ones(length(T),1);
% y = 0.05*y+Y*ones(length(T),1);
% 
% T = [0; 0.15; T+0.15];
% iL = [X(1); X(1); iL]; 
% y = [Y; Y; y]; 
% 
% figure;
% 
% subplot(3,1,1)
% plot([0, 0.15, 0.1501, 0.2], [Delta, Delta, Delta+0.5, Delta+0.5])
% axis([0.145, 0.195, 0.45, 0.6])
% 
% subplot(3,1,2)
% plot(T,iL)
% 
% 
% subplot(3,1,3)
% plot(T,y)

% z = -1000;
% %p_plus_i = tf([1,200, 5e04],[1,1000,0]);
% p1 = 51.9e3;
% z1 = 481.4;
% p0 = 600.49;
% g = (2*pi*p0*p1*p1)/(2000*(z1*z1))
% %p_plus_i = tf([1, -z],[1,0])
% p_plus_i = tf([g, g*4*pi*z1, g*4*(pi^2)*(z1^2)],[1,4*pi*p1,4*(pi^2)*(p1^2),0])
% piboost = series(p_plus_i, boost_ss)
% k=0.003;
% piboost_closed1 = feedback(k*piboost, 1);
% piboost_closed2 = feedback(0.0010*piboost, 1)
% figure
% rlocus((1/2000) * piboost);
% figure
% bode((1/2000) * piboost);
% axis([-1500, 1000, -400, 400])
% 
% % Find unit step-response 
% T = 0:0.001:0.5;
% [y,T,x] = step(piboost_closed1,T); 
% 
% % Scale step-response and add to operating point
% iL = 0.05*x(:,1)+X(1)*ones(length(T),1);
% y = 0.05*y+Y*ones(length(T),1);
% 
% T = [0; 0.15; T+0.15];
% iL = [X(1); X(1); iL]; 
% y = [Y; Y; y]; 
% 
% figure;
% 
% subplot(3,1,1)
% plot([0, 0.15, 0.1501, 0.2], [0.75, 0.75, 0.80, 0.80])
% axis([0.145, 0.195, 0.45, 0.6])
% 
% subplot(3,1,2)
% plot(T,iL)
% 
% 
% subplot(3,1,3)
% plot(T,y)







