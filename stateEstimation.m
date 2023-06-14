clear all 
close all
clc
format compact
%%system definition 
A = [-0.4 -1;1 0]
B = [1;0]
C = [0 1]
D = 0 
x0 = [-0;0] % inital condition 

%Compute gain L of the observer 
%Checking the observability and reachability
%%
%bot reachability and observability matrix rank shoukd be full 

M_r = ctrb(A, B)
rho_M = rank(M_r)
Mo = obsv(A,C)
rho_Mo= rank(Mo)



%%
%unitary DC gain with matrix N 

s_hat = 0.08
t_s2 = 4
zeta = abs(log(s_hat)) / sqrt(pi^2 +(log(s_hat))^2)
wn = log((1/100)^(-1)/(zeta*t_s2))
lambda1 = -zeta *wn + j*wn*sqrt(1-zeta^2)
lambda2 = -zeta *wn - j*wn*sqrt(1-zeta^2)
lamda_des = [lambda1 lambda2]
K = place(A,B, lamda_des)


A_c = A-B*K
B_c =B
C_c = C
D_c = 0

sys_c = ss(A_c, B_c, C_c, D_c)
N = 1/dcgain(sys_c)
%%
%OBSERVER DESIGH ->define the eig of the observer, how?
% obs eig must be choosen reala and coinsedent, damping=1 no oscillation 

%placing the eigenvalues in non oscillative fashion
lambda_obs_des = [-zeta*wn -zeta*wn]*5; %-zeta*wn was the real part of eigenvalues computed before
%why we multiplied by the constant 5, time cont of the observer is 5 time
%...faster ?????
L = acker(A', C', lambda_obs_des)' %L has to be column don't forget to transpose!!


%define state equation of the observer with the relevant matrices
                % A ,     B,     C,    D
sys_obsv = ss(A - L*C, [B L], eye(2), 0)
sys_x = ss(A,B,eye(2),0)
t_sim=50

%observer is a dyn system x_hat so we should have an initital estimated
%...state. Since we dont know it we sat x_hat to 0
x0_hat = [0;0];% initial estimated state
%%
%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


out= sim("lesson26ex2sim.slx");
figure
stepinfo(out.y.data, out.y.time, 1, 0,'RiseTimeLimits', [0 1], 'SettlingTimeTreshold', 0.05) 
plot(out.y.time, out.y.data);title('Transient requirements')
grid on
%%


figure
plot(out.x.time, out.x.data(:,1), 'b', 'linewidth',  1.5) ; hold on;
grid on;
plot(out.x_hat.time, out.x_hat.data(:,1), 'r', 'linewidth',  1.5) ; grid on;

figure
plot(out.x.time, out.x.data(:,2), 'b', 'linewidth',  1.5) ; hold on;
grid on;
plot(out.x_hat.time, out.x_hat.data(:,2), 'r', 'linewidth',  1.5) ; grid on;
return
%there are two inputs of the observer u and y


