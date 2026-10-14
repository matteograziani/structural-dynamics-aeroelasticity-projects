%  STRUCTURAL DYNAMICS AND AEROELASTICITY — ASSIGNMENT 2
%  Politecnico di Milano — 2025/2026
%  Prof. Giuseppe Quaranta
%  Student : Matteo Marwan Graziani
%  Matricola: 307177
%  MATLAB   : R2024b
%  Date     : May, 2026

%  DESCRIPTION:
%  Ritz-Galerkin aeroelastic model of a strut-braced wing. The analyses
%  cover:
%    1.a  Modal convergence study (structural, no aerodynamics)
%    1.b  Aeroelastic FRF — aileron deflection beta -> tip displacement w_P
%    1.c  Flutter speed as a function of strut chordwise position
%    1.d  Root bending moment response to a step vertical gust (direct recovery)
%    1.e  RMS root bending moment under unit-intensity white-noise gust
%         (Lyapunov equation)

clear; clc; close all;

set(groot, 'DefaultAxesFontName',  'Times New Roman', ...
           'DefaultAxesFontSize',  11, ...
           'DefaultTextFontName',  'Times New Roman', ...
           'DefaultLineLineWidth', 1.6, ...
           'DefaultAxesLineWidth', 0.8, ...
           'DefaultAxesBox',       'on', ...
           'DefaultAxesGridAlpha', 0.35, ...
           'DefaultFigureColor',   'w');

C1 = [0.122, 0.467, 0.706];   % blue
C2 = [0.839, 0.153, 0.157];   % red
C3 = [0.173, 0.627, 0.173];   % green
C4 = [0.580, 0.404, 0.741];   % purple


%% PROBLEM PARAMETERS

% Structural parameters 
EJ = 1e7;             % [N m^2]  Bending stiffness
GJ = 1e6;             % [N m^2]  Torsional stiffness
EA = 1e8;             % [N] Strut axial stiffness
b = 14;               % [m] Wing semi-span
c = 1;                % [m] Chord
ba = 2.5;             % [m] Aileron spanwise extent
e = 0.25;             % [m] EA-AC distance (= c/4)
d = 0.15;             % [m] CG-EA offset (= 0.15 c), CG aft of EA
m = 27;               % [kg/m] Mass per unit span
J_cg = 7.56;          % [kg m^2/m] Polar inertia per unit span about CG
J = J_cg + m*d^2;     % [kg m^2/m] Polar inertia about EA (parallel-axis)
gamma = 15 * pi / 180;% [rad]   Strut angle w.r.t. spanwise axis

% Tip tank
MT = 100;               % [kg]    Tank mass
L1 = 18;                % position of 2nd surname letter 'R' in alphabet
L2 =  1;                % position of 3rd surname letter 'A' in alphabet
yT = (L1*b)/26;         % [m]     Tank spanwise station
zT = -1 - 0.5*(L2/26);  % [m]     Tank vertical offset (below EA)

% Aerodynamic parameters
rho = 1.225;            % [kg/m^3] Air density (ISA, sea level)
CL_alpha = 2*pi;        % [-] Lift-curve slope
CMAC_beta = -0.1;       % [-] Pitching-moment coefficient vs aileron
CL_beta = pi;           % [-] Lift coefficient vs aileron deflection

% Strut equivalent stiffness
L_strut = (b/2)/cos(gamma);             % [m]   Strut length
k_eq = (EA/L_strut)*(sin(gamma))^2;     % [N/m] Equivalent vertical spring
xB_init = -0.5*c;                       % [m]   Default strut: leading edge


% SHAPE FUNCTIONS (polynomial Ritz basis)

% Bending modes : 
% Phi(y,i) — satisfy clamped-free BCs at y=0 and y=b 
% Phi_yy(y,i) — second spatial derivative of Phi
% Torsion modes : 
% Psi(y,j) — satisfy zero-twist BC at root (y=0)
% Psi_y — first  spatial derivative of Psi 

Phi = @(y,i) (y/b).^(i+1) - (i*(i+1))/((i+2)*(i+3)) .* (y/b).^(i+2) + (i*(i-1))/((i+3)*(i+4)) .* (y/b).^(i+3);

Phi_yy = @(y,i) (1/b^2) .* (i*(i+1) .*(y/b).^(i-1) - (i*(i+1)^2)/(i+3).*(y/b).^i + (i*(i-1)*(i+2))/(i+4) .*(y/b).^(i+1) );

Psi = @(y,j) (y/b).^j - (j/(j+1)).*(y/b).^(j+1);

Psi_y  = @(y,j) (j/b) .* ( (y/b).^(j-1) - (y/b).^j );


%  POINT 1.a — MODAL CONVERGENCE STUDY

fprintf('\nPOINT 1.a — Modal Convergence \n');

N_max = 10;
freq_bend = zeros(N_max, 1);   % 1st bending natural frequency [rad/s]
freq_tors = zeros(N_max, 1);   % 1st torsion natural frequency [rad/s]

for n = 1:N_max
    Nw = n;            % number of bending shape functions
    Nt = n;            % number of torsion shape functions
    N_tot = Nw + Nt;
    
    M_n = zeros(N_tot);
    K_n = zeros(N_tot);
    
    % Mass matrix
    % M_ww : bending-bending inertia (distributed mass + point tank mass)
    for i = 1:Nw
        for k = 1:Nw
            fun = @(y) Phi(y,i).*Phi(y,k);
            M_n(i,k) = m*integral(fun,0,b) + MT*Phi(yT,i)*Phi(yT,k);
        end
    end
    % M_wt & M_tw : bending-torsion inertia coupling (CG offset d)
    for i = 1:Nw
        for j = 1:Nt
            fun = @(y) Phi(y,i).*Psi(y,j);
            M_n(i,Nw+j) = -m*d*integral(fun,0,b);
            M_n(Nw+j,i) = M_n(i,Nw+j);
        end
    end
    % M_tt : torsion-torsion inertia (polar inertia about EA + tank contribution)
    for j = 1:Nt
        for k = 1:Nt
            fun = @(y) Psi(y,j).*Psi(y,k);
            M_n(Nw+j,Nw+k) = J*integral(fun,0,b) + MT*zT^2*Psi(yT,j)*Psi(yT,k);
        end
    end
    
    % Stiffness matrix
    % K_ww : bending stiffness (EJ) + strut equivalent spring at y=b/2
    for i = 1:Nw
        for k = 1:Nw
            fun = @(y) Phi_yy(y,i).*Phi_yy(y,k);
            K_n(i,k) = EJ*integral(fun,0,b) + k_eq*Phi(b/2,i)*Phi(b/2,k);
        end
    end
    % K_wt & K_tw : bending-torsion coupling through strut eccentricity xB
    for i = 1:Nw
        for j = 1:Nt
            K_n(i,Nw+j) = -xB_init*k_eq*Phi(b/2,i)*Psi(b/2,j);
            K_n(Nw+j,i) = K_n(i,Nw+j);
        end
    end
    % K_tt : torsional stiffness (GJ) + strut moment arm squared
    for j = 1:Nt
        for k = 1:Nt
            fun = @(y) Psi_y(y,j).*Psi_y(y,k);
            K_n(Nw+j,Nw+k) = GJ*integral(fun,0,b) + xB_init^2*k_eq*Psi(b/2,j)*Psi(b/2,k);
        end
    end
    
    % Solve generalised eigenvalue problem K q = omega^2 M q
    [~, D] = eig(K_n, M_n);
    f_Hz   = sort(sqrt(abs(diag(D)))) / (2*pi);
    
    freq_bend(n) = f_Hz(1);   % 1st bending frequency [Hz]
    freq_tors(n) = f_Hz(2);   % 1st torsion frequency [Hz]
end

% Relative convergence error between successive N values
err_bend = zeros(N_max,1);
err_tors = zeros(N_max,1);
for n = 2:N_max
    err_bend(n) = abs((freq_bend(n)-freq_bend(n-1))/freq_bend(n));
    err_tors(n) = abs((freq_tors(n)-freq_tors(n-1))/freq_tors(n));
end

fprintf('  1st bending frequency (N=%d): %.4f Hz\n', N_max, freq_bend(N_max));
fprintf('  1st torsion frequency (N=%d): %.4f Hz\n', N_max, freq_tors(N_max));

% Figure 1 : Convergence
fig1 = figure('Name','1a_Convergence','Position',[80 80 860 520]);
ax1 = subplot(2,1,1);
semilogy(2:N_max, err_bend(2:end), '-o', ...
    'Color',C1,'MarkerFaceColor',C1,'MarkerSize',6); hold on;
yline(1e-3,'--','Color',C2,'LineWidth',1.3,'Label','  \epsilon = 10^{-3}', ...
    'LabelHorizontalAlignment','left');
ylabel('\epsilon_w  [–]','Interpreter','tex');
title('Convergence of the First Natural Frequencies — Bending & Torsion', ...
      'FontWeight','bold','FontSize',12);
legend('Relative error — bending','Tolerance 10^{-3}','Location','northeast');
ylim([1e-4 1e-1]); grid on; set(ax1,'XTickLabel','');

ax2 = subplot(2,1,2);
semilogy(2:N_max, err_tors(2:end), '-s', ...
    'Color',C2,'MarkerFaceColor',C2,'MarkerSize',6); hold on;
yline(1e-3,'--','Color',C2,'LineWidth',1.3,'Label','  \epsilon = 10^{-3}', ...
    'LabelHorizontalAlignment','left');
xlabel('Number of shape functions  N  [–]','Interpreter','tex');
ylabel('\epsilon_\theta  [–]','Interpreter','tex');
legend('Relative error — torsion','Tolerance 10^{-3}','Location','northeast');
ylim([1e-4 1e0]); grid on;

linkaxes([ax1 ax2],'x'); xlim([2 N_max]);
set([ax1 ax2],'XTick',2:N_max);


%  POINT 1.b — AEROELASTIC FRF   Beta -> w_P   (Bode plot)

fprintf('\n POINT 1.b — Aeroelastic FRF \n');

Nw = 4;              % selected after convergence study
Nt = 8;
N_tot = Nw + Nt;

% Reassemble M and K with final Nw, Nt
M = zeros(N_tot);
K = zeros(N_tot);

for i = 1:Nw
    for k = 1:Nw
        fun = @(y) Phi(y,i).*Phi(y,k);
        M(i,k) = m*integral(fun,0,b) + MT*Phi(yT,i)*Phi(yT,k);
    end
end
for i = 1:Nw
    for j = 1:Nt
        fun = @(y) Phi(y,i).*Psi(y,j);
        M(i,Nw+j) = -m*d*integral(fun,0,b);
        M(Nw+j,i) = M(i,Nw+j);
    end
end
for j = 1:Nt
    for k = 1:Nt
        fun = @(y) Psi(y,j).*Psi(y,k);
        M(Nw+j,Nw+k) = J*integral(fun,0,b) + MT*zT^2*Psi(yT,j)*Psi(yT,k);
    end
end

for i = 1:Nw
    for k = 1:Nw
        fun = @(y) Phi_yy(y,i).*Phi_yy(y,k);
        K(i,k) = EJ*integral(fun,0,b) + k_eq*Phi(b/2,i)*Phi(b/2,k);
    end
end
for i = 1:Nw
    for j = 1:Nt
        K(i,Nw+j) = -xB_init*k_eq*Phi(b/2,i)*Psi(b/2,j);
        K(Nw+j,i) = K(i,Nw+j);
    end
end
for j = 1:Nt
    for k = 1:Nt
        fun = @(y) Psi_y(y,j).*Psi_y(y,k);
        K(Nw+j,Nw+k) = GJ*integral(fun,0,b) + xB_init^2*k_eq*Psi(b/2,j)*Psi(b/2,k);
    end
end

% Structural natural frequencies (no aerodynamics)
[~, D_struct] = eig(K, M);
f_struct = sort(sqrt(abs(diag(D_struct))));
omega_w1     = f_struct(1);   
omega_theta1 = f_struct(2);   

fprintf('  Structural 1st bending   ω₁ = %.4f rad/s  (%.4f Hz)\n', ...
    omega_w1, omega_w1/(2*pi));
fprintf('  Structural 1st torsion   ω₂ = %.4f rad/s  (%.4f Hz)\n', ...
    omega_theta1, omega_theta1/(2*pi));

% Aerodynamic matrices at U_inf = 25 m/s
U_bode = 25;                        
q_bode = 0.5*rho*U_bode^2;         
K_aero = zeros(N_tot);
C_aero = zeros(N_tot);

for i = 1:Nw
    for j = 1:Nt
        fun = @(y) Phi(y,i).*Psi(y,j);
        K_aero(i,Nw+j) = q_bode * c * CL_alpha * integral(fun,0,b);
    end
end
for j = 1:Nt
    for k = 1:Nt
        fun = @(y) Psi(y,j).*Psi(y,k);
        K_aero(Nw+j,Nw+k) = q_bode * c * CL_alpha * e * integral(fun,0,b);
    end
end

for i = 1:Nw
    for k = 1:Nw
        fun = @(y) Phi(y,i).*Phi(y,k);
        C_aero(i,k) = 0.5 * rho * U_bode * c * CL_alpha * integral(fun,0,b);
    end
    for j = 1:Nt
        fun = @(y) Phi(y,i).*Psi(y,j);
        C_aero(i,Nw+j) = 0.5*rho*U_bode*c*CL_alpha*e*integral(fun,0,b);
        C_aero(Nw+j,i) = C_aero(i,Nw+j);
    end
end
for j = 1:Nt
    for k = 1:Nt
        fun = @(y) Psi(y,j).*Psi(y,k);
        C_aero(Nw+j,Nw+k) = 0.5 * rho * U_bode * c * CL_alpha * e^2 * integral(fun,0,b);
    end
end

% Aileron forcing vector F_beta
F_beta = zeros(N_tot,1);
for i = 1:Nw
    fun = @(y) Phi(y,i);
    F_beta(i) = q_bode*c*CL_beta*integral(fun, b-ba, b);
end
for j = 1:Nt
    fun = @(y) Psi(y,j);
    F_beta(Nw+j) = (q_bode*c*CL_beta*e + q_bode*c^2*CMAC_beta) * integral(fun, b-ba, b);
end

% Output selector: w_P = w(b) + (c/2) theta(b) 
H_P = zeros(1,N_tot);
for i = 1:Nw 
    H_P(i)    = Phi(b,i);         
end
for j = 1:Nt
    H_P(Nw+j) = (c/2)*Psi(b,j); 
end

% FRF sweep
freq_Hz = logspace(log10(1/(2*pi)), log10(100/(2*pi)), 300);
Omega   = 2*pi*freq_Hz;
H_Wp    = zeros(length(Omega),1);
K_eff   = K - K_aero;

for k = 1:length(Omega)
    Z = -Omega(k)^2 * M + 1j*Omega(k)*C_aero + K_eff;
    H_Wp(k) = H_P * (Z \ F_beta);
end

% Figure 2 : Bode plot
fig2 = figure('Name','1b_Bode','Position',[100 80 900 580]);
ax3 = subplot(2,1,1);
semilogx(Omega, 20*log10(abs(H_Wp)), 'Color',C1); hold on;
xline(omega_w1,     '--','Color',C2, 'LineWidth',1.2, ...
    'Label',sprintf('  \\omega_{w1} = %.2f rad/s', omega_w1), ...
    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
xline(omega_theta1, '--','Color',C3, 'LineWidth',1.2, ...
    'Label',sprintf('  \\omega_{\\theta1} = %.2f rad/s', omega_theta1), ...
    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
ylabel('|H(\Omega)|  [dB]','Interpreter','tex');
title('Bode Plot — Aeroelastic FRF:  \beta  \rightarrow  w_P  (U_\infty = 25 m/s)', ...
      'Interpreter','tex','FontWeight','bold','FontSize',12);
legend('Aeroelastic FRF','1st Bending \omega_{w1}','1st Torsion \omega_{\theta1}', ...
       'Location','southwest');
grid on; xlim([Omega(1) Omega(end)]);

ax4 = subplot(2,1,2);
phase_deg = unwrap(angle(H_Wp))*180/pi;
semilogx(Omega, phase_deg, 'Color',C1); hold on;
xline(omega_w1,     '--','Color',C2,'LineWidth',1.2);
xline(omega_theta1, '--','Color',C3,'LineWidth',1.2);
xlabel('\Omega  [rad/s]','Interpreter','tex');
ylabel('\angle H(\Omega)  [deg]','Interpreter','tex');
legend('\angle H(\Omega)','Location','southwest');
grid on; xlim([Omega(1) Omega(end)]);
linkaxes([ax3 ax4],'x');


%  POINT 1.c — FLUTTER SPEED vs STRUT CHORDWISE POSITION


fprintf('\n POINT 1.c — Flutter Analysis \n');

xb_vec     = linspace(-0.5*c, 0.5*c, 50);
Uf_vec     = NaN(length(xb_vec),1);
U_scan     = 10:0.5:1000;         
n_bisect   = 12;                 
InvM       = M \ eye(N_tot);

for idx_xb = 1:length(xb_vec)
    xB = xb_vec(idx_xb);
    
    K_curr = K;
    for i = 1:Nw
        for j = 1:Nt
            K_curr(i,Nw+j) = -xB*k_eq*Phi(b/2,i)*Psi(b/2,j);
            K_curr(Nw+j,i) = K_curr(i,Nw+j);
        end
    end
    for j = 1:Nt
        for k = 1:Nt
            fun = @(y) Psi_y(y,j).*Psi_y(y,k);
            K_curr(Nw+j,Nw+k) = GJ*integral(fun,0,b) + xB^2*k_eq*Psi(b/2,j)*Psi(b/2,k);
        end
    end
    
    for idx_u = 1:length(U_scan)
        u     = U_scan(idx_u);
        q_dyn = 0.5*rho*u^2;
        KA = K_aero * (q_dyn / q_bode);
        CA = C_aero * (u    / U_bode );
        A_ss = [zeros(N_tot),           eye(N_tot);  -InvM*(K_curr - KA),    -InvM*CA  ];
        
        if any(real(eig(A_ss)) > 0)
            u_lo = U_scan(max(1,idx_u-1));
            u_hi = u;
            for iter = 1:n_bisect
                u_mid   = (u_lo + u_hi)/2;
                q_mid   = 0.5*rho*u_mid^2;
                KA_mid  = K_aero * (q_mid / q_bode);
                CA_mid  = C_aero * (u_mid / U_bode);
                A_mid   = [zeros(N_tot),             eye(N_tot);  -InvM*(K_curr - KA_mid),  -InvM*CA_mid];
                
                if any(real(eig(A_mid)) > 0)
                    u_hi = u_mid;
                else
                    u_lo = u_mid;
                end
            end
            Uf_vec(idx_xb) = (u_lo + u_hi)/2;
            break;
        end
    end
end

[Uf_min, idx_min] = min(Uf_vec);
xb_worst = xb_vec(idx_min);

fprintf('  Minimum flutter speed   V_F      = %.3f m/s\n',   Uf_min);
fprintf('  Most critical position  x_B/c    = %.4f\n\n',     xb_worst/c);

% Figure 3 : Flutter speed curve
fig3 = figure('Name','1c_Flutter','Position',[120 80 860 500]);
plot(xb_vec/c, Uf_vec, '-', 'Color',C1); hold on;
plot(xb_worst/c, Uf_min, 'o', 'MarkerEdgeColor',C2, ...
    'MarkerFaceColor',C2,'MarkerSize',9,'DisplayName','Critical point');
xline(xb_worst/c,'--','Color',C2,'LineWidth',1.1, ...
    'Label',sprintf('  x_B/c = %.3f',xb_worst/c), ...
    'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom');
yline(Uf_min,    '--','Color',C2,'LineWidth',1.1, ...
    'Label',sprintf('  V_F = %.1f m/s',Uf_min), ...
    'LabelHorizontalAlignment','right');

% Modificata la posizione di questo riquadro (in alto a sinistra)
annotation('textbox',[0.15 0.70 0.28 0.15], ...
    'String',{sprintf('Minimum flutter speed'), ...
              sprintf('$V_F = %.2f$ m/s', Uf_min), ...
              sprintf('$x_B/c = %.3f$',   xb_worst/c)}, ...
    'Interpreter','latex','FontSize',11, ...
    'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.8);

xlabel('Strut chordwise position  x_B/c  [–]','Interpreter','latex','FontSize',13);
ylabel('Flutter speed  $U_F$  [m/s]','Interpreter','latex','FontSize',13);
title('\textbf{Flutter Speed vs.\ Strut Chordwise Position}  ($N_w=4,\; N_t=8$)', ...
    'Interpreter','latex','FontSize',13);
legend({'Flutter speed $U_F$','Critical point (min $V_F$)'}, ...
    'Interpreter','latex','Location','best');
xlim([-0.5 0.5]);
y_pad = 5;
ylim([floor(min(Uf_vec)/10)*10 - y_pad, ceil(max(Uf_vec)/10)*10 + y_pad]);
grid on;


%  POINT 1.d — STEP GUST RESPONSE  (direct recovery of M_x at root)

fprintf('\n POINT 1.d — Step Gust Response \n');

U_gust   = 0.75 * Uf_min;          
q_gust   = 0.5*rho*U_gust^2;       
fprintf('  Flight speed  U_∞ = %.3f m/s  (75%% of V_F)\n', U_gust);

K_gust = K;
for i = 1:Nw
    for j = 1:Nt
        K_gust(i,Nw+j) = -xb_worst*k_eq*Phi(b/2,i)*Psi(b/2,j);
        K_gust(Nw+j,i) = K_gust(i,Nw+j);
    end
end
for j = 1:Nt
    for k = 1:Nt
        fun = @(y) Psi_y(y,j).*Psi_y(y,k);
        K_gust(Nw+j,Nw+k) = GJ*integral(fun,0,b) + xb_worst^2*k_eq*Psi(b/2,j)*Psi(b/2,k);
    end
end

KA_gust = K_aero * (q_gust / q_bode);
CA_gust = C_aero * (U_gust / U_bode);

A_state = [zeros(N_tot),  eye(N_tot); -InvM*(K_gust - KA_gust),   -InvM*CA_gust];

F_gust = zeros(N_tot,1);
for i = 1:Nw
    fun = @(y) Phi(y,i);
    F_gust(i)    = 0.5*rho*U_gust*c*CL_alpha*integral(fun,0,b);
end
for j = 1:Nt
    fun = @(y) Psi(y,j);
    F_gust(Nw+j) = 0.5*rho*U_gust*c*CL_alpha*e*integral(fun,0,b);
end
B_state = [zeros(N_tot,1); InvM*F_gust];

C_state = zeros(1, 2*N_tot);
for i = 1:Nw
    C_state(i) = EJ * Phi_yy(0, i);
end
D_state = 0;

sys_gust = ss(A_state, B_state, C_state, D_state);

% Ridotto il tempo di simulazione per dilatare la risposta in ascissa
t_sim    = linspace(0, 3, 1000);
[Mx_root, t_out] = step(sys_gust, t_sim);

fprintf('  Peak |M_x(0,t)| = %.3e N·m\n\n', max(abs(Mx_root)));

% Figure 4 : Step gust response
fig4 = figure('Name','1d_GustResponse','Position',[140 80 860 440]);
plot(t_out, Mx_root, 'Color',C3); hold on;
yline(0,'k-','LineWidth',0.8);
xlabel('Time  [s]');
ylabel('M_x(0,t)  [N·m]');
title(sprintf(['Root Bending Moment — Step Gust Response\\newline', ...
    'U_{\\infty} = %.1f m/s (75%% V_F),  v_G = step(t),  1 m/s'], U_gust), ...
    'FontWeight','bold','FontSize',12,'Interpreter','tex');
% Aggiunto xlim per vedere più chiaramente la risposta transitoria
xlim([0 3]);


%  POINT 1.e — STOCHASTIC GUST ANALYSIS  (Lyapunov equation)

fprintf(' POINT 1.e — Stochastic (White-Noise) Gust Analysis \n');
P_state = lyap(A_state, B_state * B_state');
sigma2_Mx = real( C_state * P_state * C_state' );

fprintf('  Variance of M_x(0,t)   σ²  = %.4e  (N·m)²\n',  sigma2_Mx);
fprintf('  Standard deviation      σ   = %.4e  N·m\n\n',   sqrt(sigma2_Mx));


%  FINAL OUTPUT

fprintf('  Shape functions    Nw = %d,  Nt = %d \n', Nw, Nt);
fprintf(' ω₁ bending  = %7.4f rad/s  (%6.4f Hz) \n', omega_w1,     omega_w1/(2*pi));
fprintf(' ω₁ torsion  = %7.4f rad/s  (%6.4f Hz) \n', omega_theta1, omega_theta1/(2*pi));
fprintf(' Min flutter speed   V_F   = %6.2f m/s \n', Uf_min);
fprintf(' Critical strut pos  xB/c  = %7.4f    \n', xb_worst/c);
fprintf('  Gust speed (75%% VF)  U∞   = %6.2f m/s \n', U_gust);
fprintf('   σ² (Mx at root)     = %.4e (N·m)² \n', sigma2_Mx);
