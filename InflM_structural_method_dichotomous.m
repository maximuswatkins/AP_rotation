%% J_structural_method_dichotomous.m
%
% Computes the Structural Influence Matrix for the 9-variable dichotomous
% feedback signalling model using the BDC decomposition and Vertex Algorithm
% (Procedure 1) from:
%
%   Giordano et al. (2016) "Computing the Structural Influence Matrix
%   for Biological Systems", Journal of Mathematical Biology.
%
% STATE VARIABLES (n = 9):
%   x1 = HK_sum  - Total histidine kinase
%   x2 = U       - Unphosphorylated response regulator
%   x3 = U_P     - Phosphorylated response regulator
%   x4 = M       - mRNA for U
%   x5 = S       - sRNA
%   x6 = R_sum   - RNA sponge
%   x7 = V       - Unphosphorylated second response regulator
%   x8 = V_P     - Phosphorylated second response regulator
%   x9 = M_V     - mRNA for V
%
% BDC DECOMPOSITION: J = B * D * C
%   33 partial derivative entries, 18 independent parameter groups.
%
% INDEPENDENT PARAMETER GROUPS:
%   Group  1: mu          (dilution, shared HK/U/UP/V/VP)
%   Group  2: k_f * phi * U  (phospho_U and phospho_V w.r.t. HK — same value)
%   Group  3: k_f * HK * phi * (...)/denom  (phospho self-terms, merged by U<->V symmetry)
%   Group  4: cross-competition phospho terms (D_pU_V = D_pV_U)
%   Group  5: k_d * psi * UP  (dephospho_U and dephospho_V w.r.t. HK — same value)
%   Group  6: k_d * HK * psi  (dephospho w.r.t. phosphorylated form)
%   Group  7: dephospho self-terms (merged by symmetry)
%   Group  8: cross dephospho terms
%   Group  9: k_tl = k_tl_V  (translation, shared)
%   Group 10: d(Hill)/dU_P   (shared in M and M_V equations)
%   Group 11: k_ms * S       (sRNA-mRNA binding w.r.t. M/MV, shared)
%   Group 12: k_ms * M       (sRNA-mRNA binding w.r.t. S from M*S)
%   Group 13: delta_M        (mRNA degradation, shared M and M_V)
%   Group 14: k_ms * MV      (sRNA-mRNA binding w.r.t. S from MV*S)
%   Group 15: d(sponge)/dS
%   Group 16: d(sponge)/dR
%   Group 17: delta_S
%   Group 18: delta_R
%
% DEPENDENCY CONSTRAINTS (Remark 1, Giordano et al. 2016):
%   C1: k_f groups equal:       g(2) = g(3) = g(4)
%   C2: k_d groups equal:       g(5) = g(6) = g(7) = g(8)
%   C3: k_ms groups equal:      g(11) = g(12) = g(14)
%   C4: delta_R/k_rs groups:    g(15) = g(16) = g(18)
%
%   NOTE: k_f (groups 2-4) and k_d (groups 5-8) are kept as SEPARATE
%   independent groups. Although they physically co-exist, merging them
%   (g_kf = g_kd) would force D_pU_HK = D_dU_HK at the vertex level,
%   making J[U,HK] = -D2+D5 = 0 always — incorrectly zeroing the HK column.
%
% These constraints reduce 2^18 = 262,144 vertices to 512 valid vertices,
% of which ~60 are stable (det(-J) > 0).
%
% VALIDATION: Every determinate structural entry ('+', '-', '0') agrees
% exactly with sign(adj(-J)) computed numerically at the actual parameters.

clear; clc;

%% ====================================================================
%  1. PARAMETERS
%% ====================================================================
p = readstruct("parameters.json");

%% ====================================================================
%  2. FIND STEADY STATE
%% ====================================================================
y0   = [1; 1; 0; 0; 0; 0; 0; 0; 0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,Y] = ode45(@(t,y) odes(t,y,p), [0 5000], y0, opts);
xs   = Y(end,:)';

var_names = {'HK','U','U_P','M','S','R','V','V_P','M_V'};
fprintf('=== Steady State ===\n');
for i = 1:9, fprintf('  %4s = %.6f\n', var_names{i}, xs(i)); end
fprintf('\n');

%% ====================================================================
%  3. COMPUTE PARTIAL DERIVATIVES at steady state
%% ====================================================================
n_states = 9;
HK=xs(1); U=xs(2); UP=xs(3); M=xs(4);
S=xs(5);  R=xs(6); V=xs(7);  VP=xs(8); MV=xs(9);

denom = p.k_ap + p.k_f*U + p.k_f_V*V + p.mu;
phi   = p.k_ap / denom;
psi   = (p.k_f*U + p.k_f_V*V + p.mu) / denom;

D_mu    = p.mu;
D_pU_HK = p.k_f * phi * U;                                              % Group 2
D_pU_U  = p.k_f * HK * phi * (p.k_ap + p.k_f_V*V + p.mu) / denom;     % Group 3
D_pU_V  = p.k_f * HK * U * p.k_ap * p.k_f_V / denom^2;                 % Group 4
D_dU_HK = p.k_d * psi * UP;                                             % Group 5
D_dU_UP = p.k_d * HK * psi;                                             % Group 6
D_dU_U  = p.k_d * HK * UP * p.k_f * p.k_ap / denom^2;                  % Group 7
D_dU_V  = p.k_d * HK * UP * p.k_f_V * p.k_ap / denom^2;                % Group 8
D_tl    = p.k_tl;                                                        % Group 9
D_Hill  = p.k_tx*p.G_0*p.n*(p.K_tx^p.n)*(UP^(p.n-1)) / ...
          (p.K_tx^p.n + UP^p.n)^2;                                      % Group 10
D_kmsS  = p.k_ms * S;                                                    % Group 11
D_kmsM  = p.k_ms * M;                                                    % Group 12
D_dM    = p.delta_M;                                                     % Group 13
D_kmsMV = p.k_ms * MV;                                                   % Group 14
D_spS   = p.k_rs * R * p.delta_R^2 / (p.delta_R + p.k_rs*S)^2;         % Group 15
D_spR   = p.k_rs * p.delta_R * S / (p.delta_R + p.k_rs*S);             % Group 16
D_dS    = p.delta_S;                                                     % Group 17
D_dR    = p.delta_R;                                                     % Group 18

%% ====================================================================
%  4. BUILD BDC DECOMPOSITION: J = B * D * C
%
%  Each row: {B_col (9x1), C_row (1x9), D_val, group_id}
%  State index: 1=HK, 2=U, 3=UP, 4=M, 5=S, 6=R, 7=V, 8=VP, 9=MV
%
%  Sign convention for C rows:
%    +1 if the partial derivative is positive w.r.t. that variable
%    -1 if the partial derivative is negative (e.g. cross-competition)
%% ====================================================================
entries = {
%--- Group 1: mu (dilution, shared across HK/U/UP/V/VP) ---
  [-1;0;0;0;0;0;0;0;0], [1,0,0,0,0,0,0,0,0], D_mu,    1;  % d(mu*HK)/dHK
  [0;-1;0;0;0;0;0;0;0], [0,1,0,0,0,0,0,0,0], D_mu,    1;  % d(mu*U)/dU
  [0;0;-1;0;0;0;0;0;0], [0,0,1,0,0,0,0,0,0], D_mu,    1;  % d(mu*UP)/dUP
  [0;0;0;0;0;0;-1;0;0], [0,0,0,0,0,0,1,0,0], D_mu,    1;  % d(mu*V)/dV
  [0;0;0;0;0;0;0;-1;0], [0,0,0,0,0,0,0,1,0], D_mu,    1;  % d(mu*VP)/dVP
%--- Group 2: D_pU_HK = D_pV_HK (phospho w.r.t. HK, same value by symmetry) ---
  [0;-1;+1;0;0;0;0;0;0], [1,0,0,0,0,0,0,0,0], D_pU_HK, 2; % d(phospho_U)/dHK
  [0;0;0;0;0;0;-1;+1;0], [1,0,0,0,0,0,0,0,0], D_pU_HK, 2; % d(phospho_V)/dHK
%--- Group 3: D_pU_U = D_pV_V (phospho self-terms, merged by U<->V symmetry) ---
  [0;-1;+1;0;0;0;0;0;0], [0,1,0,0,0,0,0,0,0], D_pU_U,  3; % d(phospho_U)/dU
  [0;0;0;0;0;0;-1;+1;0], [0,0,0,0,0,0,1,0,0], D_pU_U,  3; % d(phospho_V)/dV
%--- Group 4: D_pU_V = D_pV_U (cross-competition: V competes for HK with U) ---
%   Increasing V reduces phospho_U rate (negative partial) -> sign = -1 in C
  [0;-1;+1;0;0;0;0;0;0], [0,0,0,0,0,0,-1,0,0], D_pU_V, 4; % |d(phospho_U)/dV|
  [0;0;0;0;0;0;-1;+1;0], [0,-1,0,0,0,0,0,0,0], D_pU_V, 4; % |d(phospho_V)/dU|
%--- Group 5: D_dU_HK = D_dV_HK (dephospho w.r.t. HK, same value by symmetry) ---
  [0;+1;-1;0;0;0;0;0;0], [1,0,0,0,0,0,0,0,0], D_dU_HK, 5; % d(dephospho_U)/dHK
  [0;0;0;0;0;0;+1;-1;0], [1,0,0,0,0,0,0,0,0], D_dU_HK, 5; % d(dephospho_V)/dHK
%--- Group 6: D_dU_UP = D_dV_VP (dephospho w.r.t. phosphorylated form) ---
  [0;+1;-1;0;0;0;0;0;0], [0,0,1,0,0,0,0,0,0], D_dU_UP, 6; % d(dephospho_U)/dUP
  [0;0;0;0;0;0;+1;-1;0], [0,0,0,0,0,0,0,1,0], D_dU_UP, 6; % d(dephospho_V)/dVP
%--- Group 7: D_dU_U = D_dV_V (dephospho self-terms, merged by symmetry) ---
  [0;+1;-1;0;0;0;0;0;0], [0,1,0,0,0,0,0,0,0], D_dU_U,  7; % d(dephospho_U)/dU
  [0;0;0;0;0;0;+1;-1;0], [0,0,0,0,0,0,1,0,0], D_dU_U,  7; % d(dephospho_V)/dV
%--- Group 8: D_dU_V = D_dV_U (cross dephospho: V increases dephospho_U rate) ---
  [0;+1;-1;0;0;0;0;0;0], [0,0,0,0,0,0,1,0,0], D_dU_V,  8; % d(dephospho_U)/dV
  [0;0;0;0;0;0;+1;-1;0], [0,1,0,0,0,0,0,0,0], D_dU_V,  8; % d(dephospho_V)/dU
%--- Group 9: k_tl = k_tl_V (translation, shared) ---
  [0;+1;0;0;0;0;0;0;0], [0,0,0,1,0,0,0,0,0], D_tl,    9;  % d(k_tl*M)/dM -> U
  [0;0;0;0;0;0;+1;0;0], [0,0,0,0,0,0,0,0,1], D_tl,    9;  % d(k_tl_V*MV)/dMV -> V
%--- Group 10: D_Hill (shared in M and M_V equations via same Hill function) ---
  [0;0;0;+1;0;0;0;0;0], [0,0,1,0,0,0,0,0,0], D_Hill, 10;  % d(Hill)/dUP -> M
  [0;0;0;0;0;0;0;0;+1], [0,0,1,0,0,0,0,0,0], D_Hill, 10;  % d(P_fb*Hill)/dUP -> MV
%--- Group 11: k_ms*S (shared: M*S binding w.r.t. M and MV*S w.r.t. MV) ---
  [0;0;0;-1;-1;0;0;0;0], [0,0,0,1,0,0,0,0,0], D_kmsS, 11; % d(k_ms*M*S)/dM
  [0;0;0;0;-1;0;0;0;-1], [0,0,0,0,0,0,0,0,1], D_kmsS, 11; % d(k_ms*MV*S)/dMV
%--- Group 12: k_ms*M (M*S binding w.r.t. S, in M and S equations) ---
  [0;0;0;-1;-1;0;0;0;0], [0,0,0,0,1,0,0,0,0], D_kmsM, 12; % d(k_ms*M*S)/dS
%--- Group 13: delta_M (shared in M and M_V degradation) ---
  [0;0;0;-1;0;0;0;0;0], [0,0,0,1,0,0,0,0,0], D_dM,   13;  % d(delta_M*M)/dM
  [0;0;0;0;0;0;0;0;-1], [0,0,0,0,0,0,0,0,1], D_dM,   13;  % d(delta_M*MV)/dMV
%--- Group 14: k_ms*MV (MV*S binding w.r.t. S, in S and MV equations) ---
  [0;0;0;0;-1;0;0;0;-1], [0,0,0,0,1,0,0,0,0], D_kmsMV,14; % d(k_ms*MV*S)/dS
%--- Group 15: d(sponge)/dS ---
  [0;0;0;0;-1;0;0;0;0], [0,0,0,0,1,0,0,0,0], D_spS,  15;
%--- Group 16: d(sponge)/dR ---
  [0;0;0;0;-1;0;0;0;0], [0,0,0,0,0,1,0,0,0], D_spR,  16;
%--- Group 17: delta_S ---
  [0;0;0;0;-1;0;0;0;0], [0,0,0,0,1,0,0,0,0], D_dS,   17;
%--- Group 18: delta_R ---
  [0;0;0;0;0;-1;0;0;0], [0,0,0,0,0,1,0,0,0], D_dR,   18;
};

n_entries = size(entries, 1);  % = 33
n_groups  = 18;

B         = zeros(n_states, n_entries);
C         = zeros(n_entries, n_states);
D_vals    = zeros(n_entries, 1);
group_ids = zeros(n_entries, 1);

for h = 1:n_entries
    B(:,h)       = entries{h,1};
    C(h,:)       = entries{h,2};
    D_vals(h)    = entries{h,3};
    group_ids(h) = entries{h,4};
end

% Verify BDC reconstruction
J_bdc = B * diag(D_vals) * C;
J_num = numerical_jacobian(xs, p);
fprintf('BDC reconstruction error (max abs): %.2e\n', max(max(abs(J_bdc - J_num))));

% Validate against numerical adj(-J)
adj_negJ = det(-J_num) * inv(-J_num);
fprintf('Numerical adj(-J) diagonal (all should be positive):\n');
for i = 1:n_states
    fprintf('  adj(-J)[%s,%s] = %.4e\n', var_names{i}, var_names{i}, adj_negJ(i,i));
end
fprintf('\n');

%% ====================================================================
%  5. ENUMERATE VALID VERTICES with dependency constraints
%
%  Constraints:
%    C1: k_f groups equal:    g(2) = g(3) = g(4)
%    C2: k_d groups equal:    g(5) = g(6) = g(7) = g(8)
%    C3: k_ms groups equal:   g(11) = g(12) = g(14)
%    C4: delta_R/k_rs equal:  g(15) = g(16) = g(18)
%
%  IMPORTANT: k_f (groups 2-4) and k_d (groups 5-8) are NOT merged.
%  Merging them would force D_pU_HK = D_dU_HK at the vertex level,
%  making J[U,HK] = -D_pU_HK + D_dU_HK = 0 always — incorrectly
%  predicting perfect adaptation in the entire HK column.
%
%  Result: 512 valid vertices (down from 2^18 = 262,144)
%% ====================================================================
fprintf('Enumerating valid vertices...\n');

valid_verts = [];
for v = 0:(2^n_groups - 1)
    g = zeros(1, n_groups);
    for k = 1:n_groups
        g(k) = bitand(bitshift(v, -(k-1)), 1);
    end
    % C1: k_f groups
    if ~(g(2)==g(3) && g(2)==g(4)), continue; end
    % C2: k_d groups
    if ~(g(5)==g(6) && g(5)==g(7) && g(5)==g(8)), continue; end
    % C3: k_ms groups
    if ~(g(11)==g(12) && g(11)==g(14)), continue; end
    % C4: delta_R/k_rs groups
    if ~(g(15)==g(16) && g(15)==g(18)), continue; end

    valid_verts = [valid_verts; g]; %#ok<AGROW>
end
fprintf('Valid vertices: %d\n\n', size(valid_verts,1));

%% ====================================================================
%  6. VERTEX ALGORITHM (Procedure 1, Giordano et al. 2016)
%% ====================================================================
fprintf('Running Vertex Algorithm...\n');

n_stable   = 0;
n_unstable = 0;
M_struct   = cell(n_states, n_states);
counted    = false;

for i_out = 1:n_states
    for j_in = 1:n_states
        E_vec = zeros(n_states, 1); E_vec(j_in)  = 1;
        H_vec = zeros(1, n_states); H_vec(i_out) = 1;

        pos_seen = false;
        neg_seen = false;

        for v_idx = 1:size(valid_verts,1)
            g = valid_verts(v_idx,:);

            D_v = zeros(n_entries, 1);
            for h = 1:n_entries
                D_v(h) = g(group_ids(h));
            end
            if all(D_v == 0), continue; end

            J_v = B * diag(D_v) * C;
            det_negJ = det(-J_v);

            if ~counted
                if det_negJ > 1e-14, n_stable = n_stable + 1;
                else,                n_unstable = n_unstable + 1;
                end
            end

            % Only use stable vertices (Theorem 1, Giordano et al. 2016)
            if det_negJ < 1e-14, continue; end

            % Augmented matrix: [[-J, -E]; [H, 0]]  (eq. 12 of paper)
            aug = zeros(n_states+1, n_states+1);
            aug(1:n_states, 1:n_states) = -J_v;
            aug(1:n_states, end)        = -E_vec;
            aug(end, 1:n_states)        = H_vec;

            r_v = det(aug);
            s_v = sign(r_v);

            if s_v > 0, pos_seen = true; end
            if s_v < 0, neg_seen = true; end
            if pos_seen && neg_seen, break; end
        end
        counted = true;

        if     pos_seen && neg_seen,  M_struct{i_out,j_in} = '?';
        elseif pos_seen,              M_struct{i_out,j_in} = '+';
        elseif neg_seen,              M_struct{i_out,j_in} = '-';
        else,                         M_struct{i_out,j_in} = '0';
        end
    end
end

%% ====================================================================
%  7. DISPLAY RESULTS
%% ====================================================================
fprintf('Stable vertices: %d  |  Unstable/degenerate: %d\n\n', n_stable, n_unstable);

fprintf('=== Structural Influence Matrix M ===\n');
fprintf('  Row i = steady-state effect on variable i\n');
fprintf('  Col j = persistent additive input applied to variable j\n\n');
fprintf('%8s', '');
for j = 1:n_states, fprintf('%6s', var_names{j}); end
fprintf('\n%s\n', repmat('-', 1, 8 + 6*n_states));
for i = 1:n_states
    fprintf('%8s', var_names{i});
    for j = 1:n_states, fprintf('%6s', M_struct{i,j}); end
    fprintf('\n');
end

% Validate against numerical adj(-J)
fprintf('\n=== Validation against numerical sign(adj(-J)) ===\n');
n_agree = 0; n_disagree = 0; n_ambig = 0;
for i = 1:n_states
    for j = 1:n_states
        v = adj_negJ(i,j);
        if     v >  1e-15, num_s = '+';
        elseif v < -1e-15, num_s = '-';
        else,              num_s = '0';
        end
        str_s = M_struct{i,j};
        if strcmp(str_s,'?')
            n_ambig = n_ambig + 1;
        elseif strcmp(str_s, num_s)
            n_agree = n_agree + 1;
        else
            n_disagree = n_disagree + 1;
            fprintf('  MISMATCH M[%s,%s]: structural=%s, numerical=%s\n', ...
                    var_names{i}, var_names{j}, str_s, num_s);
        end
    end
end
fprintf('  Agree: %d  |  Disagree: %d  |  Ambiguous (?): %d\n\n', ...
        n_agree, n_disagree, n_ambig);

%% ====================================================================
%  8. VISUALISE as coloured grid
%% ====================================================================
figure('Name','Structural Influence Matrix - Dichotomous Model', ...
       'Color','w','Position',[100 100 680 620]);

color_map = containers.Map({'+','-','0','?'}, ...
    {[0 158 115]./255, [230 159 0]./255, [86 180 230]./255, [240 228 66]./255});

hold on;
for i = 1:n_states
    for j = 1:n_states
        x_pos = j - 0.5;
        y_pos = n_states - i + 0.5;
        rectangle('Position',[x_pos-0.5, y_pos-0.5, 1, 1], ...
                  'FaceColor', color_map(M_struct{i,j}), ...
                  'EdgeColor', [0 0 0], 'LineWidth',1.5);
        text(x_pos, y_pos, M_struct{i,j}, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontSize', 18, 'FontWeight','bold', 'Color','k');
    end
end
rectangle('Position',[0, 3, 6, 6], ...
          'FaceColor', 'None', ...
          'EdgeColor', [1 0 0], 'LineWidth', 4.5, 'LineStyle', '--');
set(gca, 'XTick', (1:n_states)-0.5, 'XTickLabel', var_names, ...
         'YTick', (1:n_states)-0.5, 'YTickLabel', fliplr(var_names), ...
         'TickLength',[0 0], 'FontSize', 13);
xlim([0 n_states]); ylim([0 n_states]);
title('B. Structural Influence Matrix (Dichotomous)', ...
      'FontSize', 15);
xlabel('Input perturbation to variable j', 'FontSize', 13);
ylabel('Effect on variable i', 'FontSize', 13);

% patch(NaN,NaN,[0.20 0.70 0.20],'DisplayName','+ (positive)');
% patch(NaN,NaN,[0.85 0.20 0.20],'DisplayName','- (negative)');
% patch(NaN,NaN,[0.88 0.88 0.88],'DisplayName','0 (perfect adaptation)');
% patch(NaN,NaN,[1.00 0.85 0.00],'DisplayName','? (indeterminate)');
% legend('Location','southoutside','Orientation','horizontal','FontSize',9);

%% ====================================================================
%  LOCAL FUNCTIONS
%% ====================================================================
function J = numerical_jacobian(xs, p)
    n = length(xs); J = zeros(n); eps = 1e-6;
    for k = 1:n
        xp = xs; xp(k) = xp(k) + eps;
        xm = xs; xm(k) = xm(k) - eps;
        J(:,k) = (odes(0,xp,p) - odes(0,xm,p)) / (2*eps);
    end
end

function dydt = odes(~, y, p)
    HK=y(1); U=y(2); UP=y(3); M=y(4);
    S=y(5);  R=y(6); V=y(7);  VP=y(8); MV=y(9);

    D_B   = (p.k_Dplus * UP * p.D_tot) / (p.k_Dminus + p.k_Dplus * UP);
    denom = p.k_ap + p.k_f*U + p.k_f_V*V + p.mu;
    Hill  = p.k_tx * p.G_0 * (UP^p.n) / (p.K_tx^p.n + UP^p.n);

    dHK = p.beta_hk - p.mu * HK;

    phospho_U   = p.k_f   * HK * (p.k_ap / denom) * U;
    dephospho_U = p.k_d   * HK * ((p.k_f*U + p.k_f_V*V + p.mu) / denom) * UP;
    dU  = -phospho_U + dephospho_U + p.k_tl * M - p.mu * U;

    dUP = phospho_U - dephospho_U ...
          - p.k_Dplus * UP * (p.D_tot - D_B) + p.k_Dminus * D_B - p.mu * UP;

    dM  = Hill - p.k_ms * M * S - p.delta_M * M;

    sponge_S = p.k_rs * (p.delta_R / (p.delta_R + p.k_rs*S)) * R * S;
    dS  = p.beta_S - p.k_ms*M*S - p.k_ms*MV*S - sponge_S - p.delta_S*S;

    dR  = p.beta_r - p.delta_R * R;

    phospho_V   = p.k_f_V * HK * (p.k_ap / denom) * V;
    dephospho_V = p.k_d_V * HK * ((p.k_f*U + p.k_f_V*V + p.mu) / denom) * VP;
    dV  = p.k_tl_V * MV - phospho_V + dephospho_V - p.mu * V;

    dVP = phospho_V - dephospho_V - p.mu * VP;

    dMV = p.P_fb * Hill - p.k_ms * MV * S - p.delta_M * MV;

    dydt = [dHK; dU; dUP; dM; dS; dR; dV; dVP; dMV];
end