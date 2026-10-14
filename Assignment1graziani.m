% Assignment 1 (Graziani Matteo Marwan 307177)
%% INITIALIZATION (MATLAB R2024b)
clear; clc; close all;

% Structural and aerodynamic parameters (Data provided in the text) 
EJ = 1e7; % Bending stiffness [Nm^2]
GJ = 1e6; % Torsional stiffness [Nm^2]
EA = 1e8; % Strut axial stiffness [N]
Cl_alpha = 2 * pi; % Wing lift slope
Cmac_beta = -0.1; % Aileron moment coefficient
b = 14; % Wingspan [m]
c = 1; % Chord [m]
e = c/4; % Aerodynamic arm (distance AC-EA) 
ba = 2.5; % Aileron span [m]
gamma_nom = 15 * pi / 180; % Nominal strut angle
Cl_beta = pi; % Aileron lift slope

% Strut parameters
% Calculation of the wall attachment height (za) to maintain constant gamma
za = b/2 * tan(gamma_nom);

%% Point 1.a (q_D VS position of B along the chord)
% Divergence pressure analysis varying the strut position x_B 
k_eff_nom = 2 * EA * cos(gamma_nom) * sin(gamma_nom)^2 / b; 
x_range = linspace(-c/2,c/2,1800); % Range from leading edge to trailing edge
N = length(x_range);
q_D = zeros(N,1);

for i=1:N
    x_pos = x_range(i);
    % Analytical expression for divergence pressure (Ritz 1-DOF bending + 1-DOF torsion)
    q_D(i) = ((4 * EJ * GJ)/(b^4) + (k_eff_nom * x_pos^2 * EJ)/(b^3) + (k_eff_nom * GJ)/(16 * b))/ ...
        ((4 * c * e * Cl_alpha * EJ)/ (3 * b^2) + (k_eff_nom * c * e * Cl_alpha * b)/(48) + (c * x_pos * k_eff_nom * Cl_alpha * b)/(32));
end

q_D_plot = q_D; 
figure('Name', '1.a Divergence Analysis');
plot(x_range, q_D_plot, 'LineWidth', 2);
hold on;
ylim([-2e5, 6e5]);
xlim([-0.5, 0.5]);
xline(0, '--r', 'Elastic Axis', 'LineWidth', 1);
xline(-c/4, '--b', 'AC', 'LineWidth', 1);
xlabel('x_B [m]');
ylabel('q_D [Pa]');
title('Divergence Pressure vs Strut Position');
grid on;

%% Point 1.b (Control Effectiveness) 
% Study of control effectiveness varying dynamic pressure q and angle gamma 
x_fixed = -0.1;

% Simplified structural and aerodynamic stiffness matrices
K_1 = [4*EJ/b^3, 0; 0, GJ/b];
K_a = [0, c*Cl_alpha*b/4; 0, c*e*Cl_alpha*b/3];

% Load vector due to aileron deflection
f_beta = [c*Cl_beta*(b^3-(b-ba)^3)/(3*b^2); ...
          c^2*Cmac_beta*(b^2-(b-ba)^2)/(2*b) + c*e*Cl_beta*(b^2-(b-ba)^2)/(2*b)];

gamma_range = 10:2.5:20;
q_vec = linspace(0,3e5,1000);

figure('Color', 'w', 'Name', '1.b Control Effectiveness'); 
hold on; grid on;

for g = gamma_range
    g_rad = g * pi / 180;
    y_B = za / tan(g_rad); % Position y_B varies to keep the wall attachment point A fixed
    k_strut = EA * cos(g_rad) * sin(g_rad)^2 / y_B;
    
    % Additional strut stiffness matrix (bending-torsional coupling)
    K_2 = k_strut * [(y_B/b)^4, -x_fixed*(y_B/b)^3; -x_fixed*(y_B/b)^3, x_fixed^2*(y_B/b)^2];
    
    Ec = zeros(size(q_vec));
    for j = 1:length(q_vec)
        q = q_vec(j);
        % Solving the static aeroelastic system: (K_tot - q*Ka) * u = q * f_beta * beta
        u = ( (K_1 + K_2) - q*K_a ) \ (q * f_beta);
        L_R = q * c * Cl_beta * ba; % Rigid wing lift
        L_theta = q * c * Cl_alpha * u(2) * b / 2; % Additional lift due to elastic deformation
        Ec(j) = (L_R + L_theta) / L_R; % Effectiveness ratio 
    end
    plot(q_vec, Ec, 'LineWidth', 2, 'DisplayName', sprintf('\\gamma = %g°', g));
end

ylim([-80, 80]);
xlim([0, 6e4]);
yline(0, 'k--', 'Reversal');
xlabel('q [Pa]'); ylabel('E_C');
legend('Location', 'southwest'); title('1.b) Control Effectiveness');

%% Point 1.c (Deformed Shapes Comparison) 
% Comparison between clean wing and braced wing at q = 0.4 * q_D_unbraced 
q_D_nostrut = (3 * GJ) / (c * e * Cl_alpha * b^2);
q_target = 0.4 * q_D_nostrut; 
beta_0 = 10 * pi / 180; 

% Solution for the unbraced case
u_nostrut = (K_1 - q_target * K_a) \ (q_target * beta_0 * f_beta);

% Solution for the braced case (nominal configuration)
y_B_nom = za / tan(gamma_nom);
k_nom = EA * cos(gamma_nom) * sin(gamma_nom)^2 / y_B_nom;
K_2_nom = k_nom * [(y_B_nom/b)^4, -x_fixed*(y_B_nom/b)^3; -x_fixed*(y_B_nom/b)^3, x_fixed^2*(y_B_nom/b)^2];
u_strut = ( (K_1 + K_2_nom) - q_target * K_a ) \ (q_target * beta_0 * f_beta);

% Reconstruction of bending and torsional shapes along the span y
y_span = 0:0.1:b;
w_n = (y_span/b).^2 * u_nostrut(1);
th_n = (y_span/b) * u_nostrut(2);
w_s = (y_span/b).^2 * u_strut(1);
th_s = (y_span/b) * u_strut(2);

figure('Color', 'w', 'Name', '1.c Comparison');
subplot(2,1,1);
plot(y_span, w_s, 'b-', 'LineWidth', 2); hold on;
plot(y_span, w_n, 'r--', 'LineWidth', 2); grid on;
ylabel('Bending w [m]'); title('1.c) Bending Shape Comparison');
legend('With Strut (Blue)', 'Without Strut (Red)');

subplot(2,1,2);
plot(y_span, th_s*180/pi, 'b-', 'LineWidth', 2); hold on;
plot(y_span, th_n*180/pi, 'r--', 'LineWidth', 2); grid on;
xlabel('y [m]'); ylabel('\theta [deg]'); title('1.c) Torsional Shape Comparison');

%% Point 1.d (Reversal Analysis) 
% Control reversal pressure analysis 
x_rev_range = linspace(-c/2,c/2,1000);
N_rev = length(x_rev_range);
gamma = 15 * pi / 180;
k_eff = 2 * EA * cos(gamma) * sin(gamma)^2 / b; 
q_R_braced = zeros(N_rev,1);

% Theoretical reversal pressure for unbraced wing
q_R_unbraced = (GJ / b) / (c * Cl_alpha * e * b/3 - 5.6 * 1.5605) * ones(N_rev,1); 

% Inclusion of aileron effects in the aerodynamic matrix for eigenvalue calculation
K_beta = [0, -5.6 * 6.535; 0, -5.6 * 1.5605];
K_a_rev = K_a + K_beta;

for i=1:N_rev
    x = x_rev_range(i);
    K_2 = [ k_eff / 16, - k_eff * x / 8;  - k_eff * x / 8, k_eff * x^2 / 4 ];
    K_tot = K_1 + K_2;
    % Reversal pressure is the eigenvalue of the system (K_structural) * u = q * (K_total_aerodynamic) * u
    eigvals_nom = eig(K_tot, K_a_rev);
    eigvals_nom = real(eigvals_nom(abs(imag(eigvals_nom)) < 1e-6));
    q_R_braced(i) = min(eigvals_nom);
end

% Calculation of target point x_B = -0.1
x_target = -0.1;
K_2_tgt = [ k_eff / 16, - k_eff * x_target / 8;  - k_eff * x_target / 8, k_eff * x_target^2 / 4 ];
K_tot_tgt = K_1 + K_2_tgt;
eigvals_tgt = eig(K_tot_tgt, K_a_rev);
q_R_target = min(real(eigvals_tgt(abs(imag(eigvals_tgt)) < 1e-6)));

% Identification of the asymptote where effectiveness tends to infinity/changes sign
den_fun = @(x) (32 * 5.6 * 1.5605 * EJ) ./ (k_eff .* x .* b^3) + (5.6 * 1.5605) ./ (2 .* x) + (5.6 * 6.535);
x_asymp = fzero(den_fun, -0.05);

figure('Color', 'w', 'Name', '1.d Control Reversal vs Strut Position');
hold on; grid on;
q_R_plot = q_R_braced;
q_R_plot(abs(q_R_plot) > 1e6) = NaN; 
plot(x_rev_range, q_R_unbraced, '--r', 'LineWidth', 1.5, 'DisplayName', 'q_R unbraced');
plot(x_rev_range, q_R_plot, '-b', 'LineWidth', 1.5, 'DisplayName', 'q_R braced');
plot(x_target, q_R_target, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Target (x_B = -0.1)');
xline(x_asymp, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
yline(0, '-k', 'HandleVisibility', 'off');
ylim([-5e5 5e5]); xlim([-0.5 0.5]);
xlabel('Strut chordwise position x_B [m]'); ylabel('q_R [Pa]');
title('1.d: Control Reversal vs Strut Position');
legend('Location', 'northeast');

%% Point 1.e (Convergence Analysis) 
% Verification of q_D convergence by increasing the number of shape functions 
N_max = 10; 
qD_convergence = zeros(N_max, 1);
N_vector = 1:N_max;

for N_shape = 1:N_max
    % Definition of polynomial shape functions (Ritz)
    phi = cell(N_shape, 1);
    phi_dotdot = cell(N_shape, 1);
    psi = cell(N_shape, 1); 
    psi_dot = cell(N_shape, 1);
    for i = 1:N_shape
        phi{i}        = @(y) (y/b).^(i+1); % Bending (C0 and C1 zero at y=0)
        phi_dotdot{i} = @(y) ((i+1)*i/b^2) .* (y/b).^(i-1);
        psi{i}        = @(y) (y/b).^i;     % Torsion (C0 zero at y=0)
        psi_dot{i}    = @(y) (i/b) .* (y/b).^(i-1);
    end
    
    % Matrix assembly via numerical integration along the span
    K_1_loop = zeros(2*N_shape, 2*N_shape);
    K_a_loop = zeros(2*N_shape, 2*N_shape);
    for i = 1:N_shape
        for j = 1:N_shape
            K_1_loop(i,j) = integral(@(y) EJ .* phi_dotdot{i}(y) .* phi_dotdot{j}(y), 0, b);
            K_1_loop(N_shape+i, N_shape+j) = integral(@(y) GJ .* psi_dot{i}(y) .* psi_dot{j}(y), 0, b);
            K_a_loop(i, N_shape+j) = c * Cl_alpha * integral(@(y) phi{i}(y) .* psi{j}(y), 0, b);
            K_a_loop(N_shape+i, N_shape+j) = c * e * Cl_alpha * integral(@(y) psi{i}(y) .* psi{j}(y), 0, b);
        end
    end
    
    % Evaluation of shape functions at attachment point B (y=b/2)
    v_phi = zeros(N_shape,1); v_psi = zeros(N_shape,1);
    for n = 1:N_shape
        v_phi(n) = phi{n}(b/2); v_psi(n) = psi{n}(b/2); 
    end
    K_2_loop = k_eff * [ (v_phi * v_phi'), (0.1 * v_phi * v_psi'); (0.1 * v_psi * v_phi'), (0.01 * v_psi * v_psi') ];
    
    % Calculation of q_D as the minimum positive eigenvalue
    eigvals = eig(K_1_loop + K_2_loop, K_a_loop);
    qD_convergence(N_shape) = min(real(eigvals(real(eigvals)>0)));
end

figure('Color', 'w', 'Name', '1.e Convergence');
plot(N_vector, qD_convergence, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'y');
grid on; xlabel('Number of Shape Functions (N)'); ylabel('q_D [Pa]');
title('1.e) Convergence of the Divergence Solution');

%% Point 1.f (Design Optimization) 
% Optimization of angle gamma to increase q_D by 40% 
N_f = 5;
x_f = -0.1;

% Calculation of base matrices for N=5
phi = cell(N_f, 1); 
phi_dotdot = cell(N_f, 1);
psi = cell(N_f, 1); 
psi_dot = cell(N_f, 1);
for i = 1:N_f
    phi{i}        = @(y) (y/b).^(i+1);
    phi_dotdot{i} = @(y) ((i+1)*i/b^2) .* (y/b).^(i-1);
    psi{i}        = @(y) (y/b).^i;
    psi_dot{i}    = @(y) (i/b) .* (y/b).^(i-1);
end

K_s1 = zeros(2*N_f, 2*N_f);
K_a = zeros(2*N_f, 2*N_f);
for i = 1:N_f
    for j = 1:N_f
        K_s1(i,j) = integral(@(y) EJ .* phi_dotdot{i}(y) .* phi_dotdot{j}(y), 0, b);
        K_s1(N_f+i, N_f+j) = integral(@(y) GJ .* psi_dot{i}(y) .* psi_dot{j}(y), 0, b);
        K_a(i, N_f+j)     = c * Cl_alpha * integral(@(y) phi{i}(y) .* psi{j}(y), 0, b);
        K_a(N_f+i, N_f+j)   = c * e * Cl_alpha * integral(@(y) psi{i}(y) .* psi{j}(y), 0, b);
    end
end

angle_range = linspace(15, 40, 25);
M = length(angle_range);
q_D_new = zeros(M, 1);
L_strut_new = zeros(M, 1);

for i = 1:M
    gamma = angle_range(i) * pi / 180;
    L_strut_new(i) = za / sin(gamma); % Variable strut length
    k_eff_new = EA * sin(gamma)^2 / L_strut_new(i);
    y_B = za / tan(gamma);
    v_phi = zeros(N_f, 1); 
    v_psi = zeros(N_f, 1);
    for j = 1:N_f
        v_phi(j) = phi{j}(y_B);
        v_psi(j) = psi{j}(y_B);
    end
    K_s2_new = k_eff_new * [ (v_phi * v_phi'), ( (-x_f) .* (v_psi * v_phi') ); ( (-x_f) .* (v_phi * v_psi') ), ( (x_f)^2 .* (v_psi * v_psi') ) ];
    K_tot_new = K_s1 + K_s2_new;
    eigvals_nom = eig(K_tot_new, K_a);
    q_D_new(i) = min(real(eigvals_nom(real(eigvals_nom) > 0)));
end

% Calculation of Baseline (gamma = 15°) for reference
v_phi = zeros(N_f, 1); 
v_psi = zeros(N_f, 1);
for j = 1:N_f
    v_phi(j) = phi{j}(b/2);
    v_psi(j) = psi{j}(b/2);
end
gamma_b = 15 * pi / 180;
k_eff_base = 2 * EA * cos(gamma_b) * sin(gamma_b)^2 /b;
K_s2_base = k_eff_base * [ (v_phi * v_phi'), ( (-x_f) .* (v_psi * v_phi') ); ( (-x_f) .* (v_phi * v_psi') ), ( (x_f)^2 .* (v_psi * v_psi') ) ];
K_tot_base = K_s1 + K_s2_base;
eigvals_nom = eig(K_tot_base, K_a);
q_D_baseline = min(real(eigvals_nom(real(eigvals_nom) > 0)));

% Find the angle that guarantees q_D_target = 1.4 * q_D_baseline
q_target_val = 1.4 * q_D_baseline; 
gamma_target = interp1(q_D_new, angle_range, q_target_val);

figure('Color', 'w', 'Name', '1.f Design Optimization Plot');
grid on; hold on;
yyaxis left
plot(angle_range, q_D_new, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'q_D');
ylabel('Divergence Pressure q_D [Pa]');
ylim([1.5e4, 2.4e4]);
set(gca, 'YColor', 'b');

yyaxis right
plot(angle_range, L_strut_new, '--or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'L_s');
ylabel('Strut Length L_s [m]');
ylim([2.5, 7.5]);
set(gca, 'YColor', 'r');

yyaxis left
plot(gamma_target, q_target_val, 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
text(gamma_target + 0.5, q_target_val, 'Optimal Design', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Strut Angle \gamma [deg]');
title('5: Divergence and Geometry vs \gamma (x_B = -0.1, N = 5)');
legend('q_D', 'L_s', 'Location', 'northeast');