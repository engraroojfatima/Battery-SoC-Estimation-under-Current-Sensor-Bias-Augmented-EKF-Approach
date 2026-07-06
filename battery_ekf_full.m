%% ===================== FULL SCRIPT: Battery + Coulomb Counting + Bias-Augmented EKF (FIXED) =====================
clear; clc; close all;

%% Battery Parameters (Lithium-ion 50Ah example)
C_nominal = 50 * 3600;   % Coulomb capacity (As)
SoC_init  = 0.8;
eta       = 0.98;
R0        = 0.01;
RC.R1     = 0.015;
RC.C1     = 2000;

%% Simulation time
t_end = 3600;
dt    = 1;
t     = (0:dt:t_end)';
N     = length(t);

%% OCV-SoC lookup table
SoC_table = 0:0.05:1;
OCV_table = 3.0 + 1.0*SoC_table - 0.3*(1-SoC_table).^2;
dOCV_dSoC_table = gradient(OCV_table, SoC_table);

%% Current profile (dynamic load)
I = zeros(N,1);
for k = 1:N
    if t(k) < 600
        I(k) = 10;
    elseif t(k) < 1800
        I(k) = 25;
    elseif t(k) < 3000
        I(k) = 40;
    else
        I(k) = 15;
    end
end

%% Simulated current sensor bias + noise
I_bias  = 3;
I_noisy = I + I_bias + 1.5*randn(N,1);

%% ---------- TRUE SYSTEM SIMULATION (ground truth) ----------
alpha = exp(-dt/(RC.R1*RC.C1));
beta  = RC.R1*(1-alpha);

SoC_true = zeros(N,1);
V1_true  = zeros(N,1);
Vt_true  = zeros(N,1);

SoC_true(1) = SoC_init;
V1_true(1)  = 0;
Vt_true(1)  = interp1(SoC_table, OCV_table, SoC_true(1)) - I(1)*R0 - V1_true(1);

for k = 1:N-1
    SoC_true(k+1) = SoC_true(k) - eta*dt/C_nominal * I(k);
    V1_true(k+1)  = alpha*V1_true(k) + beta*I(k);
    Vt_true(k+1)  = interp1(SoC_table, OCV_table, SoC_true(k+1)) - I(k+1)*R0 - V1_true(k+1);
end

Vt_measured = Vt_true + 0.03*randn(N,1);

%% ---------- COULOMB COUNTING SoC (uses biased/noisy current, no correction) ----------
SoC_cc = zeros(N,1);
SoC_cc(1) = SoC_init;
for k = 1:N-1
    SoC_cc(k+1) = SoC_cc(k) - eta*dt/C_nominal * I_noisy(k);
    SoC_cc(k+1) = min(max(SoC_cc(k+1),0),1);
end

%% ---------- BIAS-AUGMENTED EKF SETUP ----------
% State vector: x = [SoC; V1; I_bias_est]
% FIX 1: Q(3,3) was 1e-8 (way too small) -> bias state could barely move.
%        Raised to 1e-4 so the filter can actually learn the bias.
% FIX 2: P(3,3) initial uncertainty raised to reflect real uncertainty in bias.
% FIX 3: R_meas tightened slightly relative to Q so measurement is trusted
%        enough to pull the bias state, without being so tight it fights noise.
Q      = diag([1e-8, 1e-6, 1e-4]);   % process noise (bias variance increased)
R_meas = 1e-3;                        % measurement noise covariance
P      = diag([1e-4, 1e-4, 1]);       % initial state covariance (bias uncertainty increased)

x_est = [SoC_init; 0; 0];             % initial estimate, bias starts at 0

SoC_ekf  = zeros(N,1);
V1_ekf   = zeros(N,1);
bias_ekf = zeros(N,1);

SoC_ekf(1)  = x_est(1);
V1_ekf(1)   = x_est(2);
bias_ekf(1) = x_est(3);

%% ---------- EKF LOOP ----------
for k = 1:N-1

    u_corrected = I_noisy(k) - x_est(3);   % subtract current bias estimate

    % ---- Predict ----
    x_pred = [x_est(1) - eta*dt/C_nominal*u_corrected;
              alpha*x_est(2) + beta*u_corrected;
              x_est(3)];
    x_pred(1) = min(max(x_pred(1),0),1);   % saturate SoC prediction

    % FIX 4: The Jacobian A used for covariance propagation must be
    % evaluated with the CURRENT linearization (dSoC/dbias = +eta*dt/C,
    % dV1/dbias = +beta), NOT recomputed once outside the loop with stale
    % values. It was previously defined once before the loop using alpha/beta
    % which is fine numerically here since alpha/beta are constant, but the
    % SIGN on the bias coupling term was correct in A already
    % (A(1,3) = -eta*dt/C_nominal since u_corrected = I_noisy - bias,
    % so d(SoC_pred)/d(bias) = +eta*dt/C_nominal... this was actually
    % the wrong sign in the original code). Corrected below:
    A = [1, 0,  eta*dt/C_nominal;
         0, alpha, -beta;
         0, 0, 1];

    P_pred = A*P*A' + Q;

    % ---- Measurement Jacobian ----
    dOCV_dSoC = interp1(SoC_table, dOCV_dSoC_table, x_pred(1), 'linear','extrap');
    H = [dOCV_dSoC, -1, 0];

    % ---- Predicted measurement ----
    OCV_pred = interp1(SoC_table, OCV_table, x_pred(1), 'linear','extrap');
    Vt_pred  = OCV_pred - I(k+1)*R0 - x_pred(2);

    % ---- Update ----
    y_residual = Vt_measured(k+1) - Vt_pred;
    S = H*P_pred*H' + R_meas;
    K = P_pred*H' / S;

    x_est = x_pred + K*y_residual;
    x_est(1) = min(max(x_est(1),0),1);     % saturate SoC estimate

    P = (eye(3) - K*H)*P_pred;

    % ---- Store ----
    SoC_ekf(k+1)  = x_est(1);
    V1_ekf(k+1)   = x_est(2);
    bias_ekf(k+1) = x_est(3);
end

%% ---------- PLOT COMPARISON ----------
figure('Name','SoC Estimation Comparison','Color','w');
plot(t, SoC_true, 'k-', 'LineWidth', 2); hold on;
plot(t, SoC_cc,   'b--', 'LineWidth', 1.5);
plot(t, SoC_ekf,  'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('SoC');
legend('True SoC','Coulomb Counting','EKF Estimate','Location','best');
title('SoC Estimation: True vs Coulomb Counting vs EKF');
grid on;

figure('Name','Estimation Error','Color','w');
plot(t, SoC_true - SoC_cc, 'b--', 'LineWidth', 1.5); hold on;
plot(t, SoC_true - SoC_ekf, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('SoC Error');
legend('Coulomb Counting Error','EKF Error','Location','best');
title('SoC Estimation Error Comparison');
grid on;

figure('Name','Bias Estimation','Color','w');
plot(t, I_bias*ones(N,1), 'k--', 'LineWidth', 1.5); hold on;
plot(t, bias_ekf, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current bias (A)');
legend('True bias','EKF bias estimate','Location','best');
title('Current Sensor Bias: True vs EKF Estimated');
grid on;

%% ---------- SEND KEY VARIABLES TO BASE WORKSPACE ----------
assignin('base','t',t);
assignin('base','I',I);
assignin('base','SoC_true',SoC_true);
assignin('base','SoC_cc',SoC_cc);
assignin('base','SoC_ekf',SoC_ekf);
assignin('base','bias_ekf',bias_ekf);
assignin('base','Vt_true',Vt_true);
assignin('base','Vt_measured',Vt_measured);