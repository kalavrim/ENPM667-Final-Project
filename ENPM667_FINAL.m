%ENPM667 Final Project - Matthew Kalavritinos
%% Part C - Controllability
% State space form has been derived for the linearized system. This was
% found by hand as is shown in the report
clear;
close all;
clc;

%Comment if trying to find numerical ans
syms M m1 m2 l1 l2 g

%Leave this commented for symbolic form
% M = 10;
% m1 = 0;
% m2 = 0;
% l1 = 1;
% l2 = 1;
% g = 9.81;

%A matrix
A = [
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    0, (m1*g)/M, (m2*g)/M, 0, 0, 0;
    0, -(M+m1)*g/(M*l1), -m2*g/(M*l1), 0, 0, 0;
    0, -m1*g/(M*l2), -(M+m2)*g/(M*l2), 0, 0, 0,
    ];

%B matrix
B = [
    0;
    0;
    0;
    1/M;
    -1/(M*l1);
    -1/(M*l2);
];

C_sym = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];
disp(C_sym)
rank_smy = rank(C_sym)
if(rank_smy == 6)
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end
%% Part D - LQR Controller & Lyapunov's Indirect Method
clc;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

%A matrix
A = [
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    0, (m1*g)/M, (m2*g)/M, 0, 0, 0;
    0, -(M+m1)*g/(M*l1), -m2*g/(M*l1), 0, 0, 0;
    0, -m1*g/(M*l2), -(M+m2)*g/(M*l2), 0, 0, 0,
    ];

%B matrix
B = [
    0;
    0;
    0;
    1/M;
    -1/(M*l1);
    -1/(M*l2);
];

C_sym = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];
rank_smy = rank(C_sym)
if(rank_smy == 6)
    disp('The system is controllable.');
else
    disp('The system is not controllable');
end


%setting up LQR
Q = diag([100, 100, 100, 100, 100, 100]);
R = 1;

[K, ~, E] = lqr(A, B, Q, R);

fprintf('\n LQR Gain: \n');
disp(K);

fprintf('Closed-loop eigenvalues: \n');
disp(E);

X_0 = [0; 0.2; 0.15; 0; 0; 0;];
C = eye(6);
D = 0;

LQR_TEST = ss(A-(B*K), B, C, D);

figure();
initial(LQR_TEST, X_0)
%% Part D.2
% Tuned LQR controller
clc;

%setting up LQR
Q = diag([500, 10000, 1000, 10000, 10000, 10000]);
R = .01;

[K_lin, ~, E] = lqr(A, B, Q, R);

fprintf('\n LQR Gain: \n');
disp(K_lin);

fprintf('Closed-loop eigenvalues: \n');
disp(E);

X_0 = [0; 0.2; 0.15; 0; 0; 0;];
C = eye(6);
D = 0;

LQR_TEST = ss(A-(B*K_lin), B, C, D);

figure();
initial(LQR_TEST, X_0)
%% Part D.3
% Nonlinear System Response with LQR

function dstate = NL_system(state, K)
    % Parameters
    M = 1000;
    m1 = 100;
    m2 = 100;
    l1 = 20;
    l2 = 10;
    g = 9.81;

    % States
    x = state(1);
    theta1 = state(2);
    theta2 = state(3);
    x_dot = state(4);
    theta1_dot = state(5);
    theta2_dot = state(6);

    % State feedback control
    u = -K*state;

    % Mass matrix
    Mass = [M + m1 + m2, m1*l1*cos(theta1), m2*l2*cos(theta2);
            m1*l1*cos(theta1), m1*l1^2, 0;
            m2*l2*cos(theta2), 0, m2*l2^2];

    % Right-hand side
    RHS = [u - m1*l1*theta1_dot^2*sin(theta1) - m2*l2*theta2_dot^2*sin(theta2);
           -m1*g*l1*sin(theta1);
           -m2*g*l2*sin(theta2)];

    % Solve for accelerations
    acc = Mass\RHS;

    % Derivative of state
    dstate = zeros(6,1);
    dstate(1) = x_dot;
    dstate(2) = theta1_dot;
    dstate(3) = theta2_dot;
    dstate(4:6) = acc;
end


% 
Q = diag([500, 10000, 1000, 10000, 10000, 10000]);
R = 0.01;
[K_nonlin, ~, eigen] = lqr(A, B, Q, R);

% Initial condition
X_0 = [0; 0.3; 0.2; 0; 0; 0];

% Time vector
tspan = 0:0.1:150;

% Solve nonlinear system
[t, state_history] = ode45(@(t,state) NL_system(state,K_nonlin), tspan, X_0);

% Plotting
figure();

% x position
subplot(6,1,1);
plot(t, state_history(:,1),'k'); grid on;
ylabel('x (m)');
title('Nonlinear system response with LQR');

% theta1
subplot(6,1,2);
plot(t, state_history(:,2),'r'); grid on;
ylabel('\theta_1 (rad)');

% theta2
subplot(6,1,3);
plot(t, state_history(:,3),'b'); grid on;
ylabel('\theta_2 (rad)');

% x velocity
subplot(6,1,4);
plot(t, state_history(:,4),'k'); grid on;
ylabel('x-dot (m/s)');

% theta1 velocity
subplot(6,1,5);
plot(t, state_history(:,5),'r'); grid on;
ylabel('\theta_1-dot (rad/s)');

% theta2 velocity
subplot(6,1,6);
plot(t, state_history(:,6),'b'); grid on;
ylabel('\theta_2-dot (rad/s)');
xlabel('Time (s)');

fprintf('Closed-loop eigenvalues: \n');
disp(eigen);



%% Part E - Observability
clc;
close all; 

% C matricies:
% C1: [x(t)], C2: [theta1(t), theta2(t)], C3: [x(t), theta2(t)], 
% C4:[x(t), theta1(t), theta2(t)

C1 = [1, 0, 0, 0, 0, 0];
C2 = [0, 1, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0];
C3 = [1, 0, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0];
C4 = [1, 0, 0, 0, 0, 0;
      0, 1, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0];

% build each of the observability matricies
function O = observability_manual(A, C)
    n = size(A,1);   % number of states
    O = C;
    Ak = A;
    for k = 1:n-1
        O = [O; C*Ak];
        Ak = Ak*A;
    end
end

C_list = {C1, C2, C3, C4};
names  = {'C1','C2','C3','C4'};

for k = 1:4
    Ck = C_list{k};

    O = observability_manual(A, Ck);
    r = rank(O);

    fprintf('%s rank = %d\n', names{k}, r);
end

%% Part F - Lunenberger Observer
clc;
close all;

D = 0;

% desired poles
obs_poles = [-5, -4, -2, -8, -9, -3];

%lunenberger observers
L1 = place(A', C1', obs_poles)';
L3 = place(A', C3', obs_poles)';
L4 = place(A', C4', obs_poles)';

%closed loop matricies
A1_cl = [((A-B*K_lin)) (B*K_lin); zeros(size(A)) (A-(L1*C1))];
A3_cl = [((A-B*K_lin)) (B*K_lin); zeros(size(A)) (A-(L3*C3))];
A4_cl = [((A-B*K_lin)) (B*K_lin); zeros(size(A)) (A-(L4*C4))];

B1_cl = [B; zeros(size(B))];
B3_cl = [B; zeros(size(B))];
B4_cl = [B; zeros(size(B))];

C1_cl = [C1 zeros(size(C1))];
C3_cl = [C3 zeros(size(C3))];
C4_cl = [C4 zeros(size(C4))];

D1_cl = D;
D3_cl = D;
D4_cl = D;

%initial conditions
X_0 = [0, .3, .2, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%%
%simulating and plotting the observers to the linearized systems
obs1_sys = ss(A1_cl, B1_cl, C1_cl, D1_cl);
figure()
step(obs1_sys);
figure()
initial(obs1_sys, X_0);
grid on

obs3_sys = ss(A3_cl, B3_cl, C3_cl, D3_cl);
figure()
step(obs3_sys);
figure()
initial(obs3_sys, X_0);
grid on

obs4_sys = ss(A4_cl, B4_cl, C4_cl, D4_cl);
figure()
step(obs4_sys);
figure()
initial(obs4_sys, X_0);
grid on

%% nonlinear

clc;
close all;

% Time span
tspan = 0:0.01:200;

% Initial conditions
z0_true = [0; 0.3; 0.2; 0; 0; 0];  % nonlinear plant
zhat0 = zeros(6,1);                % observer initial state
zaug0 = [z0_true; zhat0];          % augmented initial condition

% Observers and gains
C_list = {C1, C3, C4};
L_list = {L1, L3, L4};
names = {'C1','C3','C4'};
state_names = {'x','theta1','theta2','xdot','theta1dot','theta2dot'};

% Inputs
u_IC = 0;       % IC response
u_step = 500;   % step response

for i = 1:length(C_list)
    C = C_list{i};
    L = L_list{i};
    n_outputs = size(C,1);

    % Initial Condition Response
    [t_ic, zaug_ic] = ode45(@(t,zaug) NL_plant_aug_step(t, zaug, K_lin, A, B, C, L, u_IC), tspan, zaug0);
    z_true_ic = zaug_ic(:,1:6);   
    zhat_ic   = zaug_ic(:,7:12);  

    % Step Response
    [t_step, zaug_step] = ode45(@(t,zaug) NL_plant_aug_step(t, zaug, K_lin, A, B, C, L, u_step), tspan, zaug0);
    z_true_step = zaug_step(:,1:6);
    zhat_step   = zaug_step(:,7:12);

    % Plot IC Response
    figure('Name',['IC Response - Observer ', names{i}]);
    for j = 1:n_outputs
        subplot(n_outputs,1,j);
        plot(t_ic, z_true_ic(:,find(C(j,:))), '--','LineWidth',1.5); hold on;
        plot(t_ic, zhat_ic(:,find(C(j,:))), '-','LineWidth',1.5);
        ylabel(state_names{find(C(j,:))});
        if j==n_outputs, xlabel('Time (s)'); end
        legend('True','Observer Estimate');
        grid on;
    end
    sgtitle(['Initial Condition Response - Nonlinear Plant + Observer ', names{i}]);

    % Plot Step Response
    figure('Name',['Step Response - Observer ', names{i}]);
    for j = 1:n_outputs
        subplot(n_outputs,1,j);
        plot(t_step, z_true_step(:,find(C(j,:))), '--','LineWidth',1.5); hold on;
        plot(t_step, zhat_step(:,find(C(j,:))), '-','LineWidth',1.5);
        ylabel(state_names{find(C(j,:))});
        if j==n_outputs, xlabel('Time (s)'); end
        legend('True','Observer Estimate');
        grid on;
    end
    sgtitle(['Step Response - Nonlinear Plant + Observer ', names{i}]);
end

% functions
function dzaug = NL_plant_aug_step(~, zaug, K, A, B, C, L, u_ext)
    n = 6;
    z = zaug(1:n);        % nonlinear plant
    zhat = zaug(n+1:end); % observer

    % Plant uses true state feedback
    u = -K*z + u_ext;

    % Nonlinear plant dynamics
    dz = NL_system_aug(z, K, u);

    % Observer dynamics
    y = C*z;
    dzhat = A*zhat + B*u + L*(y - C*zhat);

    dzaug = [dz; dzhat];
end

function dstate = NL_system_aug(state, K, u)
    % Parameters
    M = 1000; m1 = 100; m2 = 100; l1 = 20; l2 = 10; g = 9.81;

    % States
    x = state(1); theta1 = state(2); theta2 = state(3);
    x_dot = state(4); theta1_dot = state(5); theta2_dot = state(6);

    % Mass matrix
    Mass = [M+m1+m2, m1*l1*cos(theta1), m2*l2*cos(theta2);
            m1*l1*cos(theta1), m1*l1^2, 0;
            m2*l2*cos(theta2), 0, m2*l2^2];

    % Right-hand side
    RHS = [u - m1*l1*theta1_dot^2*sin(theta1) - m2*l2*theta2_dot^2*sin(theta2);
           -m1*g*l1*sin(theta1);
           -m2*g*l2*sin(theta2)];

    % Accelerations
    acc = Mass\RHS;

    % Derivative
    dstate = [x_dot; theta1_dot; theta2_dot; acc];
end

%% Part G LQG Controller

close all;
clc;

C=C1;
%LQR
Q = diag([500, 10000, 1000, 10000, 10000, 10000]);
R = 0.01;
[K, P, eigen] = lqr(A, B, Q, R);

%Kalman
W = 0.1*eye(6);
V = 1;
G = eye(6);

Lk = lqe(A, G, C, W, V);

fprintf('LQR K (1x6):\n'); disp(K);
fprintf('Kalman gain Lk (6x1):\n'); disp(Lk);

u_IC = 0;        % no external input for IC response
u_step = 500;    % step response

% initial condition reposnse
[t_ic, zaug_ic] = ode45(@(t,zaug) NL_plant_LQG_aug(t, zaug, K, A, B, C, Lk, u_IC), tspan, zaug0);
z_true_ic = zaug_ic(:,1:6);
zhat_ic   = zaug_ic(:,7:12);

% step response
[t_step, zaug_step] = ode45(@(t,zaug) NL_plant_LQG_aug(t, zaug, K, A, B, C, Lk, u_step), tspan, zaug0);
z_true_step = zaug_step(:,1:6);
zhat_step   = zaug_step(:,7:12);

% plots
figure('Name','LQG: IC - True vs Estimate (x)');
subplot(2,1,1);
plot(t_ic, z_true_ic(:,1), '--','LineWidth',1.5); hold on;
plot(t_ic, zhat_ic(:,1), '-','LineWidth',1.5);
ylabel('x (m)'); legend('True x','xhat'); grid on;
title('LQG - Initial Condition Response (no external force)');

subplot(2,1,2);
% Plot control input applied (u = -K*xhat) for IC simulation
u_ic_time = - (zhat_ic * K');   % zhat_ic is N x 6, K' is 6x1 -> N x 1
plot(t_ic, u_ic_time, 'LineWidth',1.2);
ylabel('u (N)'); xlabel('Time (s)'); grid on;

figure('Name','LQG: Step - True vs Estimate (x)');
subplot(2,1,1);
plot(t_step, z_true_step(:,1), '--','LineWidth',1.5); hold on;
plot(t_step, zhat_step(:,1), '-','LineWidth',1.5);
ylabel('x (m)'); legend('True x','xhat'); grid on;
title('LQG - Step Response (external force applied)');

subplot(2,1,2);
u_step_time = - (zhat_step * K');
plot(t_step, u_step_time, 'LineWidth',1.2);
ylabel('u (N)'); xlabel('Time (s)'); grid on;

% quantity checks
err_ic = max(abs(z_true_ic - zhat_ic), [], 1);    % max abs error per state (IC)
err_step = max(abs(z_true_step - zhat_step), [], 1); % max abs error per state (step)
fprintf('Max abs estimation error (IC) per state:\n'); disp(err_ic);
fprintf('Max abs estimation error (Step) per state:\n'); disp(err_step);

% LQG Estimator
function dzaug = NL_plant_LQG_aug(~, zaug, K, A, B, C, Lk, u_ext)
    n = 6;
    z = zaug(1:n);        % true nonlinear plant state
    zhat = zaug(n+1:end); % estimator state

    % Controller uses estimate
    u = -K * zhat + u_ext;

    % Nonlinear plant dynamics
    dz = NL_system_for_control(z, u);

    y = C * z;   % measurement from the nonlinear plant
    dzhat = A * zhat + B * u + Lk * (y - C * zhat);

    dzaug = [dz; dzhat];
end

function dstate = NL_system_for_control(state, u)
    
    M = 1000; m1 = 100; m2 = 100; l1 = 20; l2 = 10; g = 9.81;

    x = state(1); theta1 = state(2); theta2 = state(3);
    x_dot = state(4); theta1_dot = state(5); theta2_dot = state(6);

    % Mass matrix
    Mass = [M + m1 + m2, m1*l1*cos(theta1), m2*l2*cos(theta2);
            m1*l1*cos(theta1), m1*l1^2, 0;
            m2*l2*cos(theta2), 0, m2*l2^2];

    % RHS
    RHS = [u - m1*l1*theta1_dot^2*sin(theta1) - m2*l2*theta2_dot^2*sin(theta2);
           -m1*g*l1*sin(theta1);
           -m2*g*l2*sin(theta2)];

    acc = Mass \ RHS;   % accelerations [xddot; theta1ddot; theta2ddot]

    dstate = [x_dot; theta1_dot; theta2_dot; acc];
end