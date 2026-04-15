%% ====
%  J_numerical_method_DFB.m
%
%  Nonlinear ODE simulation (ode45), steady-state finding, Jacobian
%  construction, eigenvalue stability analysis, and numerical influence
%  matrix (sign pattern of adj(-J)) for the dichotomous feedback system.
%
%  State vector: x = [H, U, U_P, M, S, R, V, V_P, M_V]'  (n = 9)
%
%  New species vs. original:
%    V    = unphosphorylated secondary response regulator
%    V_P  = phosphorylated secondary response regulator
%    M_V  = mRNA for V (subject to same sRNA S)
%  ====

clear; clc; close all;

%% ---- 1. Parameters ----
p = struct();

% Growth / dilution
p.mu       = 0.0234;

% Histidine kinase
p.HK_ss    = 50;
p.beta_hk  = p.mu * p.HK_ss;

% Phosphotransfer (U)
p.k_f      = 6.12e3;
p.k_d      = 4.0e-3;
p.k_ap_max = 0.02;
p.I        = 1000;
p.K_da     = 2000;
p.k_ap     = p.k_ap_max * (p.I / (p.I + p.K_da));

% Phosphotransfer (V) — same as U by default
p.k_f_V    = p.k_f;
p.k_d_V    = p.k_d;

% Decoy (QSS)
p.k_Dp     = 40.0;
p.k_Dm     = 0.0126 * p.k_Dp;
p.D_sum    = 20.0;

% Translation
p.k_tl     = 1.0;
p.k_tl_V   = 1.0;

% Transcription (Hill)
p.k_tx     = 0.4;
p.G_0      = 20.0;
p.K_tx     = 0.0126;
p.n_hill   = 2.0794;

% mRNA-sRNA interaction
p.k_ms     = 0.01344;

% mRNA degradation
p.delta_M  = 0.246;

% sRNA
p.beta_S   = 100;
p.delta_S  = 0.048;

% RNA sponge
p.beta_r   = 10.0;
p.delta_R  = 0.048;
p.k_rs     = 951.36576;

% Dichotomous feedback strength
p.P_fb     = 0.4;

%% ---- 2. Simulate with ode45 ----
x0 = [1; 1; 0; 0; 0; 0; 0; 0; 0];   % [H, U, U_P, M, S, R, V, V_P, M_V]
tspan = [0 200];
opts  = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[t, x] = ode45(@(t,x) ode_system(t, x, p), tspan, x0, opts);

x_ss = x(end, :)';
var_names = {'H','U','U_P','M','S','R','V','V_P','M_V'};
n_vars = 9;

fprintf('=== Steady State from ode45 ===\n');
for k = 1:n_vars
    fprintf('  %s* = %.6f\n', var_names{k}, x_ss(k));
end

% Residual check
res = ode_system(0, x_ss, p);
fprintf('\n=== Residual |dx/dt| at steady state ===\n');
for k = 1:n_vars
    fprintf('  |d%s/dt| = %.2e\n', var_names{k}, abs(res(k)));
end

%% ---- 3. Jacobian at steady state ----
J = compute_jacobian(x_ss, p);

fprintf('\n=== Jacobian Matrix J ===\n');
disp(J);

%% ---- 4. Eigenvalue stability analysis ----
eigenvalues = eig(J);

fprintf('=== Eigenvalues of J ===\n');
for k = 1:n_vars
    fprintf('  lambda_%d = %+.6f %+.6fi\n', k, ...
        real(eigenvalues(k)), imag(eigenvalues(k)));
end

max_re = max(real(eigenvalues));
fprintf('\nMax real part: %.6f\n', max_re);
if max_re < -1e-10
    fprintf('RESULT: Steady state is ASYMPTOTICALLY STABLE.\n\n');
elseif abs(max_re) < 1e-10
    fprintf('RESULT: Steady state is MARGINALLY STABLE.\n\n');
else
    fprintf('RESULT: Steady state is UNSTABLE.\n\n');
end

%% ---- 5. Numerical influence matrix: sign(adj(-J)) ----
negJ      = -J;
det_negJ  = det(negJ);

fprintf('=== Numerical Influence Matrix ===\n');
fprintf('det(-J) = %.6e\n\n', det_negJ);

if abs(det_negJ) < 1e-15
    warning('det(-J) near zero — J may be singular.');
end

adj_negJ = det_negJ * inv(negJ); %#ok<MINV>
tol = 0;

M_sign = char(zeros(n_vars));
for i = 1:n_vars
    for j = 1:n_vars
        val = adj_negJ(i, j);
        if abs(val) == tol
            M_sign(i, j) = '0';
        elseif val > tol
            M_sign(i, j) = '+';
        else
            M_sign(i, j) = '-';
        end
    end
end

fprintf('  M_ij: influence on variable i due to perturbation of variable j\n');
fprintf('  +: positive,  -: negative,  0: no influence\n\n');
fprintf('%10s', '');
for j = 1:n_vars; fprintf('%7s', var_names{j}); end
fprintf('\n');
for i = 1:n_vars
    fprintf('%8s |', var_names{i});
    for j = 1:n_vars; fprintf('%6s ', M_sign(i,j)); end
    fprintf('\n');
end

%% ---- 6. Visualisation ----
figure('Position', [100 100 1400 900], 'Color', 'w');

% Time traces — split into two panels for clarity
subplot(2, 3, 1);
plot(t, x(:,1:3), 'LineWidth', 1.5);
legend({'H','U','U_P'}, 'Location','best');
xlabel('Time (min)'); ylabel('Conc. (\muM)');
title('Kinase & Response Regulator'); grid on;

subplot(2, 3, 2);
plot(t, x(:,4), 'm-', 'LineWidth',1.5); hold on;
plot(t, x(:,9), 'm--','LineWidth',1.5);
plot(t, x(:,5), 'g-', 'LineWidth',1.5);
legend({'M','M_V','S'}, 'Location','best');
xlabel('Time (min)'); ylabel('Conc. (\muM)');
title('RNA Species'); grid on;

subplot(2, 3, 3);
plot(t, x(:,6), 'c-', 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('Conc. (\muM)');
title('R (RNA Sponge)'); grid on;

subplot(2, 3, 4);
plot(t, x(:,7), 'Color',[0.6 0.2 0.8], 'LineWidth',1.5); hold on;
plot(t, x(:,8), 'Color',[0.9 0.5 0.9], 'LineWidth',1.5, 'LineStyle','--');
legend({'V','V_P'}, 'Location','best');
xlabel('Time (min)'); ylabel('Conc. (\muM)');
title('Secondary Regulator (Dichotomous FB)'); grid on;

% Eigenvalue plot
subplot(2, 3, 5);
plot(real(eigenvalues), imag(eigenvalues), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold on; xline(0,'k--'); yline(0,'k--');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Eigenvalues of J'); grid on;

% Jacobian heatmap
subplot(2, 3, 6);
imagesc(J); colorbar;
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names, ...
         'YTick', 1:n_vars, 'YTickLabel', var_names);
title('Jacobian J'); xlabel('Variable j'); ylabel('Equation i');

% Influence matrix heatmap
figure;
M_numeric = zeros(n_vars);
for i = 1:n_vars
    for j = 1:n_vars
        if strcmp(M_sign(i,j), '+');  M_numeric(i,j) =  1; end
        if strcmp(M_sign(i,j), '-');  M_numeric(i,j) = -1; end
    end
end
imagesc(M_numeric); colorbar; clim([-1 1]);
%colormap(gca, default(256));
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names, ...
         'YTick', 1:n_vars, 'YTickLabel', var_names);
title('Numerical Influence Matrix  sign(adj(-J))');
xlabel('Input perturbation j'); ylabel('Effect on i');
for i = 1:n_vars
    for j = 1:n_vars
        text(j, i, M_sign(i,j), 'HorizontalAlignment','center', ...
            'FontSize', 11, 'FontWeight','bold');
    end
end


%% ====
%  LOCAL FUNCTIONS
%  ====

function dxdt = ode_system(~, x, p)
    H   = x(1);  U   = x(2);  U_P = x(3);
    M   = x(4);  S   = x(5);  R   = x(6);
    V   = x(7);  V_P = x(8);  M_V = x(9);

    % Shared denominator (both U and V compete for HK active site)
    Delta = p.k_ap + p.k_f * U + p.k_f_V * V + p.mu;

    % Phosphorylation / dephosphorylation fluxes
    phospho_U   = p.k_f   * H * (p.k_ap / Delta) * U;
    dephospho_U = p.k_d   * H * ((p.k_f * U + p.k_f_V * V + p.mu) / Delta) * U_P;
    phospho_V   = p.k_f_V * H * (p.k_ap / Delta) * V;
    dephospho_V = p.k_d_V * H * ((p.k_f * U + p.k_f_V * V + p.mu) / Delta) * V_P;

    % Decoy (QSS)
    D_B = (p.k_Dp * U_P * p.D_sum) / (p.k_Dm + p.k_Dp * U_P);

    % Hill function
    Hill = p.k_tx * p.G_0 * (U_P^p.n_hill) / (p.K_tx^p.n_hill + U_P^p.n_hill);

    % Sponge absorption
    sponge_eff = p.k_rs * p.delta_R / (p.delta_R + p.k_rs * S);

    % ODEs
    dH   = p.beta_hk - p.mu * H;
    dU   = -phospho_U + dephospho_U + p.k_tl * M - p.mu * U;
    dU_P = phospho_U - dephospho_U ...
           - p.k_Dp * U_P * (p.D_sum - D_B) + p.k_Dm * D_B ...
           - p.mu * U_P;
    dM   = Hill - p.k_ms * M * S - p.delta_M * M;
    dS   = p.beta_S - p.k_ms * M * S - p.k_ms * M_V * S ...
           - sponge_eff * R * S - p.delta_S * S;
    dR   = p.beta_r - p.delta_R * R;
    dV   = p.k_tl_V * M_V - phospho_V + dephospho_V - p.mu * V;
    dV_P = phospho_V - dephospho_V - p.mu * V_P;
    dM_V = p.P_fb * Hill - p.k_ms * M_V * S - p.delta_M * M_V;

    dxdt = [dH; dU; dU_P; dM; dS; dR; dV; dV_P; dM_V];
end

function J = compute_jacobian(x_ss, p)
    H   = x_ss(1);  U   = x_ss(2);  U_P = x_ss(3);
    M   = x_ss(4);  S   = x_ss(5);  R   = x_ss(6);
    V   = x_ss(7);  V_P = x_ss(8);  %M_V = x_ss(9); % not needed explicitly

    % Shared denominator and its partial derivatives
    Delta   = p.k_ap + p.k_f * U + p.k_f_V * V + p.mu;
    dDelta_dU = p.k_f;
    dDelta_dV = p.k_f_V;

    % Shorthand fluxes at steady state
    A_U  = p.k_ap / Delta;
    B_U  = (p.k_f * U + p.k_f_V * V + p.mu) / Delta;

    % Hill function and its derivative w.r.t. U_P
    Hill_val   = p.k_tx * p.G_0 * (U_P^p.n_hill) / (p.K_tx^p.n_hill + U_P^p.n_hill);
    H_frac     = (U_P^p.n_hill) / (p.K_tx^p.n_hill + U_P^p.n_hill);
    dHill_dUP  = (Hill_val / U_P) * p.n_hill * (1 - H_frac);

    % Sponge terms
    sponge_denom = p.delta_R + p.k_rs * S;
    sponge_eff   = p.k_rs * p.delta_R / sponge_denom;
    dsponge_dS   = -p.k_rs^2 * p.delta_R / sponge_denom^2;

    % ---- Row 1: dH/dt ----
    J11 = -p.mu;

    % ---- Row 2: dU/dt ----
    % d(phospho_U)/dH = k_f * A_U * U
    % d(dephospho_U)/dH = k_d * B_U * U_P
    J21 = -p.k_f * A_U * U + p.k_d * B_U * U_P;

    % d(phospho_U)/dU: A_U changes via Delta
    dA_dU = -p.k_ap * dDelta_dU / Delta^2;
    dB_dU = (p.k_f * Delta - (p.k_f * U + p.k_f_V * V + p.mu) * dDelta_dU) / Delta^2;
    J22 = -(p.k_f * H * (A_U + U * dA_dU)) + p.k_d * H * U_P * dB_dU - p.mu;

    J23 = p.k_d * H * B_U;   % d(dephospho_U)/dU_P
    J24 = p.k_tl;             % d(k_tl*M)/dM

    % d(phospho_U)/dV and d(dephospho_U)/dV via Delta
    dA_dV = -p.k_ap * dDelta_dV / Delta^2;
    dB_dV = (p.k_f_V * Delta - (p.k_f * U + p.k_f_V * V + p.mu) * dDelta_dV) / Delta^2;
    % Note: dB_dV simplifies since d(k_f*U + k_f_V*V + mu)/dV = k_f_V
    J27 = -(p.k_f * H * U * dA_dU * (dDelta_dV/dDelta_dU)) ...  % chain rule via Delta
          + p.k_d * H * U_P * dB_dV;
    % Cleaner direct form:
    J27 = -p.k_f * H * U * dA_dV + p.k_d * H * U_P * dB_dV;

    % ---- Row 3: dU_P/dt ----
    J31 = p.k_f * A_U * U - p.k_d * B_U * U_P;
    J32 = p.k_f * H * (A_U + U * dA_dU) - p.k_d * H * U_P * dB_dU;
    J33 = -p.k_d * H * B_U - p.mu;   % decoy terms cancel under QSS
    J37 = p.k_f * H * U * dA_dV - p.k_d * H * U_P * dB_dV;

    % ---- Row 4: dM/dt ----
    J43 = dHill_dUP;
    J44 = -p.k_ms * S - p.delta_M;
    J45 = -p.k_ms * M;

    % ---- Row 5: dS/dt ----
    J54 = -p.k_ms * S;
    J55 = -p.k_ms * M - p.k_ms * x_ss(9) ...   % M_V contribution
          - (sponge_eff + dsponge_dS * S) * R - p.delta_S;
    J56 = -sponge_eff * S;
    J59 = -p.k_ms * S;   % d(dS)/dM_V

    % ---- Row 6: dR/dt ----
    J66 = -p.delta_R;

    % ---- Row 7: dV/dt ----
    % d(phospho_V)/dH = k_f_V * A_U * V  (same A_U since same Delta)
    % d(dephospho_V)/dH = k_d_V * B_U * V_P
    J71 = -p.k_f_V * A_U * V + p.k_d_V * B_U * V_P;

    % d/dU: via Delta
    dA_dU_forV = dA_dU;   % same Delta
    dB_dU_forV = dB_dU;
    J72 = -p.k_f_V * H * V * dA_dU_forV + p.k_d_V * H * V_P * dB_dU_forV;

    % d/dV
    J77 = p.k_tl_V * 0 ...   % M_V not in this row
          - p.k_f_V * H * (A_U + V * dA_dV) ...
          + p.k_d_V * H * V_P * dB_dV - p.mu;
    J77 = -p.k_f_V * H * (A_U + V * dA_dV) + p.k_d_V * H * V_P * dB_dV - p.mu;

    J78 = p.k_d_V * H * B_U;   % d(dephospho_V)/dV_P
    J79 = p.k_tl_V;             % d(k_tl_V * M_V)/dM_V

    % ---- Row 8: dV_P/dt ----
    J81 = p.k_f_V * A_U * V - p.k_d_V * B_U * V_P;
    J82 = p.k_f_V * H * V * dA_dU_forV - p.k_d_V * H * V_P * dB_dU_forV;
    J87 = p.k_f_V * H * (A_U + V * dA_dV) - p.k_d_V * H * V_P * dB_dV;
    J88 = -p.k_d_V * H * B_U - p.mu;

    % ---- Row 9: dM_V/dt ----
    J93 = p.P_fb * dHill_dUP;
    J95 = -p.k_ms * x_ss(9);   % d(-k_ms * M_V * S)/dS
    J99 = -p.k_ms * S - p.delta_M;

    % ---- Assemble 9x9 Jacobian ----
    %      H    U    U_P  M    S    R    V    V_P  M_V
    J = zeros(9);
    J(1,1) = J11;
    J(2,1) = J21;  J(2,2) = J22;  J(2,3) = J23;  J(2,4) = J24;  J(2,7) = J27;
    J(3,1) = J31;  J(3,2) = J32;  J(3,3) = J33;  J(3,7) = J37;
    J(4,3) = J43;  J(4,4) = J44;  J(4,5) = J45;
    J(5,4) = J54;  J(5,5) = J55;  J(5,6) = J56;  J(5,9) = J59;
    J(6,6) = J66;
    J(7,1) = J71;  J(7,2) = J72;  J(7,7) = J77;  J(7,8) = J78;  J(7,9) = J79;
    J(8,1) = J81;  J(8,2) = J82;  J(8,7) = J87;  J(8,8) = J88;
    J(9,3) = J93;  J(9,5) = J95;  J(9,9) = J99;
end