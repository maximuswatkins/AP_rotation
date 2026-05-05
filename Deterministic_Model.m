%% Deterministic_Model_Dichotomous.m
% Extended model with dichotomous feedback via phosphorylation sequestration
%
% STATE VARIABLES (9):
%   y(1) = HK_sum   - Total histidine kinase
%   y(2) = U        - Unphosphorylated response regulator
%   y(3) = U_P      - Phosphorylated response regulator
%   y(4) = M        - mRNA for U
%   y(5) = S        - sRNA
%   y(6) = R_sum    - RNA sponge
%   y(7) = V       - Unphosphorylated second response regulator (NEW)
%   y(8) = V_P     - Phosphorylated second response regulator (NEW)
%   y(9) = M_V     - mRNA for SR
%
% KEY MODELLING CHOICES:
%   - M_V is transcribed from the same promoter as M, scaled by P_fb
%   - M_V is degraded by the SAME sRNA S (same k_ms rate)
%   - This ensures sRNA regulates both U and V production symmetrically
%   - V is translated from M_V, then phosphorylated by HK_P
%   - V_P has no downstream transcriptional activity
%   - Production enters V (unphosphorylated), matching U/U_P convention

%% ====
%  1. PARAMETERS
%% ====
p = readstruct("parameters.json");

% *** Dichotomous feedback strength ***
% P = 0 gives the non-dichotomous model
p.P_fb = 1.0;

%% ====
%  2. INITIAL CONDITIONS
%% ====
y0 = [1.0;   % HK_sum
      1.0;   % U
      0.0;   % U_P
      0.0;   % M
      0.0;   % S
      0.0;   % R_sum
      0.0;   % V
      0.0;   % V_P
      0.0];  % M_V

%% ====
%  3. TIME SPAN
%% ====
tspan = [0, 300];

%% ====
%  4. SOLVE
%% ====
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, y] = ode45(@(t,y) odes(t, y, p), tspan, y0, opts);

% Extract state variables
HK_sum = y(:,1);
U      = y(:,2);
U_P    = y(:,3);
M      = y(:,4);
S      = y(:,5);
R_sum  = y(:,6);
V     = y(:,7);
V_P   = y(:,8);
M_V   = y(:,9);

% Bound DNA decoy
D_B = (p.k_Dplus .* U_P .* p.D_tot) ./ (p.k_Dminus + p.k_Dplus .* U_P);

% Derived: free decoy and occupancy fractions (QSS)
D_free = p.D_tot - D_B;

D_bound_frac = D_B    ./ max(p.D_tot, 1e-12);
D_free_frac  = D_free ./ max(p.D_tot, 1e-12);

% Derived sponge occupancy (QSS interpretation)
R_free  = (p.delta_R ./ (p.delta_R + p.k_rs .* S)) .* R_sum;
R_bound = R_sum - R_free;

%% ====
%  5. PLOT
%% ====
figure('Name','Signalling System with Dichotomous Feedback','NumberTitle','off','Color','w');

% --- Main panel: U ---
subplot(3,2,[1 2]);
plot(t, U_P, '-', 'Color', [216 91 0]./255, 'LineWidth', 4);
xlabel('Time (minutes)');
ylabel('Concentration (\muM)');
ylim([0 1.5])
title('U_P – Phosphorylated Response Regulator (Output)', 'FontSize', 12);
grid on;

% --- All protein species together ---
subplot(3,2,3);
plot(t, HK_sum, '-', 'Color', [255 200 135]./255,  'LineWidth',4); hold on;
plot(t, U,      '-', 'Color', [255 145 53]./255, 'LineWidth',4);
plot(t, U_P,    '-', 'Color', [230 91 0]./255, 'LineWidth',4);
plot(t, V,     '--', 'Color', [160 58 0]./255, 'LineWidth',4);
plot(t, V_P,   '--', 'Color', [0 0 0]./255, 'LineWidth',4);
legend('HK_{sum}','U','U_P','V','V_P','Location','northeast');
xlabel('Time (min.)'); ylabel('Conc. (\muM)'); title('Kinase & Response Regulators', 'FontSize', 12); grid on;
ylim([0 500])

% --- Bound decoy ---
subplot(3,2,4);
plot(t, D_B,    'Color',[255 145 53]./255, 'LineWidth',4); hold on;
plot(t, D_free, 'Color',[0 0 0], 'LineWidth',4);
xlabel('Time (min.)');
ylabel('Conc. (\muM)');
title('DNA Buffer Occupancy (QSS)', 'FontSize', 12);
legend('D_{bound}','D_{free}','Location','northeast');
grid on;

% --- RNA species ---
subplot(3,2,[5 6]);
plot(t, M,    '-', 'Color', [255 200 135]./255,'LineWidth',4); hold on;
plot(t, M_V, '--','Color', [255 145 53]./255,'LineWidth',4);
plot(t, S,    '-', 'Color', [230 91 0]./255,'LineWidth',4);
plot(t, R_bound, '-', 'Color', [160 58 0]./255,'LineWidth',4);
plot(t, R_free,  '--', 'Color', [0 0 0],'LineWidth',4);
legend('M','M_{V}','S','R_{bound}','R_{free}','Location','northeast');
xlabel('Time (min.)'); ylabel('Conc. (\muM)'); title('RNA Species', 'FontSize', 12); grid on;
ylim([0 175])


%% ====
%  5b. FIGURE 2 — Feedback strength comparison (cf. Fig. 4c)
%% ====
figure('Name','Feedback Strength Comparison','NumberTitle','off','Color','w');

P_fb_values = [0, 0.1, 0.4, 0.7, 1.0];
colors = [255 200 135;  
          255 145 53;   
          230 91 0;   
          160 58 0;   
          0.0 0.0 0.0]./255;  
leg_entries = cell(1, numel(P_fb_values));

for i = 1:numel(P_fb_values)
    p_i = p;
    p_i.P_fb = P_fb_values(i);
    y0_i = y0;
    if P_fb_values(i) == 0
        y0_i(7:9) = 0;
    end
    [t_i, y_i] = ode45(@(t,y) odes(t, y, p_i), tspan, y0_i, opts);
    plot(t_i, y_i(:,3), '-', 'Color', colors(i,:), 'LineWidth', 4); hold on;
    if P_fb_values(i) == 0
        leg_entries{i} = 'P_{fb} = 0 (no FB)';
    else
        leg_entries{i} = sprintf('P_{fb} = %.1f', P_fb_values(i));
    end
end

legend(leg_entries, 'Location', 'best');
xlabel('Time (minutes)');
ylabel('[U_P] (\muM)');
title('Effect of Dichotomous Feedback Strength on U_P', 'FontSize',12);
grid on;
box on;
set(gca, 'FontSize', 12);

%% ====
%  6. STEADY-STATE SUMMARY
%% ====
fprintf('\n=== Approximate Steady-State Values (last time point) ===\n');
fprintf('  HK_sum : %.4f\n', HK_sum(end));
fprintf('  U      : %.4f\n', U(end));
fprintf('  U_P    : %.4f\n', U_P(end));
fprintf('  M      : %.4f\n', M(end));
fprintf('  M_V   : %.4f\n', M_V(end));
fprintf('  S      : %.4f\n', S(end));
fprintf('  R_sum  : %.4f\n', R_sum(end));
fprintf('  D_B    : %.4f  (QSS)\n', D_B(end));
fprintf('  V     : %.4f\n', V(end));
fprintf('  V_P   : %.4f\n', V_P(end));
fprintf('  P_fb   : %.4f\n', p.P_fb);
fprintf('====\n\n');

%% ====
%  LOCAL FUNCTION: ODE RIGHT-HAND SIDE
%% ====
function dydt = odes(~, y, p)
    HK_sum = y(1);
    U      = y(2);
    U_P    = y(3);
    M      = y(4);
    S      = y(5);
    R_sum  = y(6);
    V     = y(7);
    V_P   = y(8);
    M_V   = y(9);

    HK_tot = HK_sum;

    % Bound decoy (QSS)
    D_B = (p.k_Dplus * U_P * p.D_tot) / (p.k_Dminus + p.k_Dplus * U_P);

    % *** MODIFIED denominator: V competes for HK_P ***
    denom = p.k_ap + p.k_f * U + p.k_f_V * V + p.mu;

    % Hill function
    Hill = p.k_tx * p.G_0 * (U_P^p.n) / (p.K_tx^p.n + U_P^p.n);

    % --- (1) HK_sum ---
    dHK_sum = p.beta_hk - p.mu * HK_sum;

    % --- (2) U ---
    phospho_U   = p.k_f * HK_tot * (p.k_ap / denom) * U;
    dephospho_U = p.k_d * HK_tot * ((p.k_f * U + p.k_f_V * V + p.mu) / denom) * U_P;
    dU = -phospho_U + dephospho_U + p.k_tl * M - p.mu * U;

    % --- (3) U_P ---
    dU_P = phospho_U - dephospho_U ...
           - p.k_Dplus * U_P * (p.D_tot - D_B) ...
           + p.k_Dminus * D_B ...
           - p.mu * U_P;

    % --- (4) M (mRNA for U) — unchanged ---
    dM = Hill - p.k_ms * M * S - p.delta_M * M;

    % --- (5) S — NOW TITRATES AGAINST BOTH M AND M_V ---
    sponge_S = p.k_rs * (p.delta_R / (p.delta_R + p.k_rs * S)) * R_sum * S;
    dS = p.beta_S - p.k_ms * M * S - p.k_ms * M_V * S - sponge_S - p.delta_S * S;

    % --- (6) R_sum ---
    dR_sum = p.beta_r - p.delta_R * R_sum;

    % --- (7) V ---
    phospho_V   = p.k_f_V * HK_tot * (p.k_ap / denom) * V;
    dephospho_V = p.k_d_V * HK_tot * ((p.k_f * U + p.k_f_V * V + p.mu) / denom) * V_P;
    dV = p.k_tl_V * M_V - phospho_V + dephospho_V - p.mu * V;

    % --- (8) V_P ---
    dV_P = phospho_V - dephospho_V - p.mu * V_P;

    % --- (9) M_V (mRNA for V) — NEW ---
    %   Same promoter as M, scaled by P_fb
    %   Degraded by same sRNA S at same rate k_ms
    dM_V = p.P_fb * Hill - p.k_ms * M_V * S - p.delta_M * M_V;

    dydt = [dHK_sum; dU; dU_P; dM; dS; dR_sum; dV; dV_P; dM_V];
end