%% ========================================================================
%  
%  Nonlinear ODE simulation (ode45), steady-state finding, Jacobian
%  construction, eigenvalue stability analysis, and numerical influence
%  matrix (sign pattern of adj(-J)).
%  
%  State vector: x = [h, u, p, m, s, r]'
%  ========================================================================

clear; clc; close all;

%% ---- 1. Define model parameters ----------------------------------------

params = readstruct("parameters.json");

%% ---- 2. Simulate the full nonlinear system with ode45 ------------------

% Initial conditions: [h0, u0, p0, m0, s0, r0]
x0 = [1; 1; 0.0; 0.0; 0.0; 0];

% Time span
tspan = [0 200];  % long enough to reach steady state

% Solve ODE
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, x] = ode45(@(t, x) ode_system(t, x, params), tspan, x0, opts);

% Extract steady state from final time point
x_ss = x(end, :)';
fprintf('=== Steady State from ode45 ===\n');
var_names = {'h', 'u', 'p', 'm', 's', 'r'};
for k = 1:6
    fprintf('  %s* = %.6f\n', var_names{k}, x_ss(k));
end
fprintf('\n');

% Verify: check that derivatives are near zero at steady state
dxdt_check = ode_system(0, x_ss, params);
fprintf('=== Residual |dx/dt| at steady state ===\n');
for k = 1:6
    fprintf('  |d%s/dt| = %.2e\n', var_names{k}, abs(dxdt_check(k)));
end
fprintf('\n');

%% ---- 3. Evaluate the Jacobian at steady state --------------------------

J = compute_jacobian(x_ss, params);

fprintf('=== Jacobian Matrix J ===\n');
disp(J);

%% ---- 4. Eigenvalue stability analysis ----------------------------------

eigenvalues = eig(J);

fprintf('=== Eigenvalues of J ===\n');
for k = 1:length(eigenvalues)
    fprintf('  lambda_%d = %.6f + %.6fi\n', k, ...
        real(eigenvalues(k)), imag(eigenvalues(k)));
end

max_real_part = max(real(eigenvalues));
fprintf('\nMax real part: %.6f\n', max_real_part);

if max_real_part < -1e-10
    fprintf('RESULT: Steady state is ASYMPTOTICALLY STABLE.\n\n');
elseif abs(max_real_part) < 1e-10
    fprintf('RESULT: Steady state is MARGINALLY STABLE (bifurcation point).\n\n');
else
    fprintf('RESULT: Steady state is UNSTABLE.\n\n');
end

%% ---- 5. Numerical influence matrix: sign(adj(-J)) ----------------------

negJ = -J;
det_negJ = det(negJ);
fprintf('=== Numerical Influence Matrix ===\n');
fprintf('det(-J) = %.6e\n', det_negJ);

if abs(det_negJ) < 1e-15
    warning('det(-J) near zero — J may be singular.');
end

adj_negJ = det_negJ * inv(negJ);  %#ok<MINV>

% Sign pattern
n_vars = 6;
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

fprintf('\n  M_ij: influence on variable i due to input on variable j\n');
fprintf('  +: positive, -: negative, 0: no influence\n\n');
fprintf('%8s', '');
for j = 1:n_vars; fprintf('%6s', var_names{j}); end
fprintf('\n');
for i = 1:n_vars
    fprintf('%6s |', var_names{i});
    for j = 1:n_vars; fprintf('%5c ', M_sign(i, j)); end
    fprintf('\n');
end
fprintf('\n');

%% ---- 6. Visualisation --------------------------------------------------

figure('Position', [100 100 1400 800]);

% Time traces
subplot(2, 3, 1);
plot(t, x, 'LineWidth', 1.5);
legend(var_names, 'Location', 'best');
xlabel('Time'); ylabel('Concentration');
title('ode45 Time Traces');
grid on;

% Eigenvalue plot
subplot(2, 3, 2);
plot(real(eigenvalues), imag(eigenvalues), 'rx', 'MarkerSize', 12, ...
    'LineWidth', 2);
hold on;
xline(0, 'k--'); yline(0, 'k--');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Eigenvalues of J');
grid on; axis equal;

% Jacobian heatmap
subplot(2, 3, 3);
imagesc(J); colorbar;
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names);
set(gca, 'YTick', 1:n_vars, 'YTickLabel', var_names);
title('Jacobian J'); xlabel('Variable j'); ylabel('Equation i');

% Influence matrix heatmap
subplot(2, 3, 4);
M_numeric = zeros(n_vars);
for i = 1:n_vars
    for j = 1:n_vars
        if M_sign(i,j) == '+'; M_numeric(i,j) = 1;
        elseif M_sign(i,j) == '-'; M_numeric(i,j) = -1;
        end
    end
end
imagesc(M_numeric); colorbar; caxis([-1 1]);
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names);
set(gca, 'YTick', 1:n_vars, 'YTickLabel', var_names);
title('Numerical Influence Matrix');
for i = 1:n_vars
    for j = 1:n_vars
        text(j, i, M_sign(i,j), 'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Phase portrait: u vs p
subplot(2, 3, 5);
plot(x(:,2), x(:,3), 'b-', 'LineWidth', 1.5);
hold on;
plot(x_ss(2), x_ss(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('u'); ylabel('p');
title('Phase Portrait: u vs p');
grid on;

% Phase portrait: m vs s
subplot(2, 3, 6);
plot(x(:,4), x(:,5), 'b-', 'LineWidth', 1.5);
hold on;
plot(x_ss(4), x_ss(5), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('m'); ylabel('s');
title('Phase Portrait: m vs s');
grid on;

sgtitle('ODE45 Simulation & Stability Analysis', 'FontSize', 14);

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function dxdt = ode_system(~, x, p)
    % Unpack state
    h = x(1); u = x(2); pp = x(3); m = x(4); s = x(5); r = x(6);
    
    % Auxiliary
    Phi = p.k_ap + p.k_f*u + p.mu;
    A   = p.k_ap / Phi;
    B   = (p.k_f*u + p.mu) / Phi;
    Hill = pp^p.n / (p.K_tx^p.n + pp^p.n);
    
    % Effective sponge absorption rate
    sponge_eff = p.k_rs * p.delta_R / (p.delta_R + p.k_rs * s);
    
    % ODEs (decoy terms cancel under QSS)
    dh = p.beta_hk - p.mu * h;
    du = -p.k_f * h * A * u + p.k_d * h * B * pp + p.k_tl * m - p.mu * u;
    dp = p.k_f * h * A * u - p.k_d * h * B * pp - p.mu * pp;
    dm = p.k_tx * p.G_0 * Hill - p.k_ms * m * s - p.delta_M * m;
    ds = p.beta_S - p.k_ms * m * s - sponge_eff * r * s - p.delta_S * s;
    dr = p.beta_r - p.delta_R * r;
    
    dxdt = [dh; du; dp; dm; ds; dr];
end

function J = compute_jacobian(x_ss, p)
    h = x_ss(1); u = x_ss(2); pp = x_ss(3);
    m = x_ss(4); s = x_ss(5); r = x_ss(6);
    
    Phi = p.k_ap + p.k_f*u + p.mu;
    A   = p.k_ap / Phi;
    B   = (p.k_f*u + p.mu) / Phi;
    H_val = pp^p.n / (p.K_tx^p.n + pp^p.n);
    
    % Row 1 (h): decoupled
    J11 = -p.mu;
    
    % Row 2 (u):
    J21 = -p.k_f * A * u + p.k_d * B * pp;
    J22 = -p.k_f * h * A * (1 - p.k_f*u/Phi) ...
          + (p.k_f * A / Phi) * p.k_d * h * pp - p.mu;
    J23 = p.k_d * h * B;
    J24 = p.k_tl;
    
    % Row 3 (p):
    J31 = p.k_f * A * u - p.k_d * B * pp;
    J32 = p.k_f * h * A * (1 - p.k_f*u/Phi) ...
          - (p.k_f * A / Phi) * p.k_d * h * pp;
    J33 = -p.k_d * h * B - p.mu;
    
    % Row 4 (m):
    dHill_dp = p.n * H_val * (1 - H_val) / pp;
    J43 = p.k_tx * p.G_0 * dHill_dp;
    J44 = -p.k_ms * s - p.delta_M;
    J45 = -p.k_ms * m;
    
    % Row 5 (s):
    J54 = -p.k_ms * s;
    J55 = -p.k_ms * m ...
          - p.k_rs * r * p.delta_R^2 / (p.delta_R + p.k_rs*s)^2 ...
          - p.delta_S;
    J56 = -p.k_rs * p.delta_R * s / (p.delta_R + p.k_rs * s);
    
    % Row 6 (r): decoupled
    J66 = -p.delta_R;
    
    % Assemble
    J = [J11,  0,    0,    0,    0,    0;
         J21,  J22,  J23,  J24,  0,    0;
         J31,  J32,  J33,  0,    0,    0;
         0,    0,    J43,  J44,  J45,  0;
         0,    0,    0,    J54,  J55,  J56;
         0,    0,    0,    0,    0,    J66];
end