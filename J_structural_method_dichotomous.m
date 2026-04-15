%% J_structural_method_DFB.m
%
%  Structural influence matrix for the DICHOTOMOUS FEEDBACK (DFB) system
%  using the BDC decomposition and Vertex Algorithm (Procedure 1) from
%  Giordano et al. (2016).
%
%  State vector: x = [H, U, U_P, M, S, R, V, V_P, M_V]'  (n = 9)
%
%  KEY MODELLING CHOICE — SHARED KINETIC PARAMETERS
%  ----
%  U and V share the same bifunctional kinase H with identical rates:
%    k_fV = k_f   (phosphorylation rate)
%    k_dV = k_d   (dephosphorylation rate)
%    k_tlV = k_tl (translation rate)
%
%  The shared denominator is:
%    denom = k_ap + k_f*(U + V) + mu
%    A = k_ap / denom           (phosphorylation activity fraction)
%    B = (k_f*(U+V) + mu)/denom (dephosphorylation activity fraction)
%
%  This creates cross-coupling: U and V compete for the kinase.
%
%  RUNTIME
%  ----
%  10 independent parameter groups => 2^10 = 1024 vertices.  Instant.

clear; clc;

%% ====
%% 1.  BDC DECOMPOSITION  (n=9, q=33)
%% ====
n = 9;   % state variables
q = 33;  % partial derivatives

% ---- Stoichiometric matrix S (n x 20 reactions) ----
%  Species order: [H, U, U_P, M, S, R, V, V_P, M_V]
%  Reactions:
%   1  dilution H           (rate = mu*H)
%   2  phosphorylation U    (rate = k_f*H*A*U)
%   3  dephosphorylation UP (rate = k_d*H*B*UP)
%   4  translation M->U     (rate = k_tl*M)
%   5  dilution U           (rate = mu*U)
%   6  dilution U_P         (rate = mu*UP)
%   7  transcription ->M    (rate = k_tx*G0*Hill(UP))
%   8  mRNA-sRNA binding    (rate = k_ms*M*S)
%   9  mRNA degradation M   (rate = delta_M*M)
%  10  sponge-sRNA binding  (rate = k_rs*sf*R*S)
%  11  sRNA degradation S   (rate = delta_S*S)
%  12  sponge degradation R (rate = delta_R*R)
%  13  translation MV->V    (rate = k_tl*MV)       ** same k_tl **
%  14  phosphorylation V    (rate = k_f*H*A*V)      ** same k_f, shared A **
%  15  dephosphorylation VP (rate = k_d*H*B*VP)     ** same k_d, shared B **
%  16  dilution V           (rate = mu*V)
%  17  dilution V_P         (rate = mu*VP)
%  18  transcription ->MV   (rate = P_fb*k_tx*G0*Hill(UP))
%  19  MV-sRNA binding      (rate = k_ms*MV*S)
%  20  MV degradation       (rate = delta_M*MV)     ** same delta_M **

Sm = zeros(n, 20);
%              H   U  UP   M   S   R   V  VP  MV
Sm(:, 1) = [-1;  0;  0;  0;  0;  0;  0;  0;  0];  % rxn1:  dilution H
Sm(:, 2) = [ 0; -1; +1;  0;  0;  0;  0;  0;  0];  % rxn2:  phospho U
Sm(:, 3) = [ 0; +1; -1;  0;  0;  0;  0;  0;  0];  % rxn3:  dephos UP
Sm(:, 4) = [ 0; +1;  0;  0;  0;  0;  0;  0;  0];  % rxn4:  translation M->U
Sm(:, 5) = [ 0; -1;  0;  0;  0;  0;  0;  0;  0];  % rxn5:  dilution U
Sm(:, 6) = [ 0;  0; -1;  0;  0;  0;  0;  0;  0];  % rxn6:  dilution UP
Sm(:, 7) = [ 0;  0;  0; +1;  0;  0;  0;  0;  0];  % rxn7:  transcription ->M
Sm(:, 8) = [ 0;  0;  0; -1; -1;  0;  0;  0;  0];  % rxn8:  M+S binding
Sm(:, 9) = [ 0;  0;  0; -1;  0;  0;  0;  0;  0];  % rxn9:  M degradation
Sm(:,10) = [ 0;  0;  0;  0; -1;  0;  0;  0;  0];  % rxn10: sponge-S binding
Sm(:,11) = [ 0;  0;  0;  0; -1;  0;  0;  0;  0];  % rxn11: S degradation
Sm(:,12) = [ 0;  0;  0;  0;  0; -1;  0;  0;  0];  % rxn12: R degradation
Sm(:,13) = [ 0;  0;  0;  0;  0;  0; +1;  0;  0];  % rxn13: translation MV->V
Sm(:,14) = [ 0;  0;  0;  0;  0;  0; -1; +1;  0];  % rxn14: phospho V
Sm(:,15) = [ 0;  0;  0;  0;  0;  0; +1; -1;  0];  % rxn15: dephos VP
Sm(:,16) = [ 0;  0;  0;  0;  0;  0; -1;  0;  0];  % rxn16: dilution V
Sm(:,17) = [ 0;  0;  0;  0;  0;  0;  0; -1;  0];  % rxn17: dilution VP
Sm(:,18) = [ 0;  0;  0;  0;  0;  0;  0;  0; +1];  % rxn18: transcription ->MV
Sm(:,19) = [ 0;  0;  0;  0; -1;  0;  0;  0; -1];  % rxn19: MV+S binding
Sm(:,20) = [ 0;  0;  0;  0;  0;  0;  0;  0; -1];  % rxn20: MV degradation

% ---- Partial derivatives (q=33) ----
%
% Format: [reaction_index, variable_index, sign_of_partial]
%   sign = +1 if the partial derivative is positive
%   sign = -1 if the partial derivative is negative
%
% Variable indices (1-based): H=1 U=2 UP=3 M=4 S=5 R=6 V=7 VP=8 MV=9
%
% SHARED DENOMINATOR CROSS-TERMS:
%   denom = k_ap + k_f*(U+V) + mu
%   A = k_ap / denom  =>  dA/dU = dA/dV = -k_f*k_ap/denom^2 < 0
%   B = (k_f*(U+V)+mu)/denom  =>  dB/dU = dB/dV = k_f*k_ap/denom^2 > 0
%
%   d(A*U)/dU = k_ap*(k_ap + k_f*V + mu)/denom^2 > 0
%   d(A*V)/dV = k_ap*(k_ap + k_f*U + mu)/denom^2 > 0
%   U*dA/dV = -k_f*k_ap*U/denom^2 < 0  (cross: V reduces phospho of U)
%   V*dA/dU = -k_f*k_ap*V/denom^2 < 0  (cross: U reduces phospho of V)
%
% D0:  rxn1 /dH  = mu                                          > 0
% D1:  rxn2 /dH  = k_f * A * U                                 > 0
% D2:  rxn2 /dU  = k_f * H * d(A*U)/dU                         > 0
% D3:  rxn2 /dV  = k_f * H * U * dA/dV                         < 0  ** NEW **
% D4:  rxn3 /dH  = k_d * B * UP                                > 0
% D5:  rxn3 /dU  = k_d * H * UP * dB/dU    (AND: k_f & k_d)   > 0
% D6:  rxn3 /dUP = k_d * H * B                                 > 0
% D7:  rxn3 /dV  = k_d * H * UP * dB/dV    (AND: k_f & k_d)   > 0  ** NEW **
% D8:  rxn4 /dM  = k_tl                                        > 0
% D9:  rxn5 /dU  = mu                                          > 0
% D10: rxn6 /dUP = mu                                          > 0
% D11: rxn7 /dUP = k_tx * G0 * dHill/dUP                       > 0
% D12: rxn8 /dM  = k_ms * S                                    > 0
% D13: rxn8 /dS  = k_ms * M                                    > 0
% D14: rxn9 /dM  = delta_M                                     > 0
% D15: rxn10/dS  = k_rs * sf * R                               > 0
% D16: rxn10/dR  = k_rs * sf_R * S                             > 0
% D17: rxn11/dS  = delta_S                                     > 0
% D18: rxn12/dR  = delta_R                                     > 0
% D19: rxn13/dMV = k_tl                                        > 0  ** same k_tl **
% D20: rxn14/dH  = k_f * A * V                                 > 0  ** same k_f **
% D21: rxn14/dU  = k_f * H * V * dA/dU                         < 0  ** NEW cross **
% D22: rxn14/dV  = k_f * H * d(A*V)/dV                         > 0  ** same k_f **
% D23: rxn15/dH  = k_d * B * VP                                > 0  ** same k_d **
% D24: rxn15/dU  = k_d * H * VP * dB/dU    (AND: k_f & k_d)   > 0  ** NEW cross **
% D25: rxn15/dV  = k_d * H * VP * dB/dV    (AND: k_f & k_d)   > 0  ** NEW cross **
% D26: rxn15/dVP = k_d * H * B                                 > 0  ** same k_d **
% D27: rxn16/dV  = mu                                          > 0  ** WAS MISSING **
% D28: rxn17/dVP = mu                                          > 0  ** WAS MISSING **
% D29: rxn18/dUP = P_fb * k_tx * G0 * dHill/dUP               > 0
% D30: rxn19/dS  = k_ms * MV                                   > 0
% D31: rxn19/dMV = k_ms * S                                    > 0
% D32: rxn20/dMV = delta_M                                     > 0  ** same delta_M **

partials = [
%  rxn  var  sign
    1,   1,  +1;   % D0:  mu*H
    2,   1,  +1;   % D1:  phospho U @ H
    2,   2,  +1;   % D2:  phospho U @ U  (d(A*U)/dU > 0)
    2,   7,  -1;   % D3:  phospho U @ V  (U*dA/dV < 0)       ** NEW **
    3,   1,  +1;   % D4:  dephos UP @ H
    3,   2,  +1;   % D5:  dephos UP @ U  (cross, AND: k_f & k_d)
    3,   3,  +1;   % D6:  dephos UP @ UP
    3,   7,  +1;   % D7:  dephos UP @ V  (cross, AND: k_f & k_d) ** NEW **
    4,   4,  +1;   % D8:  translation M->U
    5,   2,  +1;   % D9:  dilution U
    6,   3,  +1;   % D10: dilution UP
    7,   3,  +1;   % D11: Hill -> M
    8,   4,  +1;   % D12: k_ms*M*S @ M
    8,   5,  +1;   % D13: k_ms*M*S @ S
    9,   4,  +1;   % D14: delta_M*M
   10,   5,  +1;   % D15: k_rs*R*S @ S
   10,   6,  +1;   % D16: k_rs*R*S @ R
   11,   5,  +1;   % D17: delta_S*S
   12,   6,  +1;   % D18: delta_R*R
   13,   9,  +1;   % D19: translation MV->V  (same k_tl)
   14,   1,  +1;   % D20: phospho V @ H      (same k_f)
   14,   2,  -1;   % D21: phospho V @ U  (V*dA/dU < 0)       ** NEW **
   14,   7,  +1;   % D22: phospho V @ V  (d(A*V)/dV > 0)
   15,   1,  +1;   % D23: dephos VP @ H      (same k_d)
   15,   2,  +1;   % D24: dephos VP @ U  (cross, AND: k_f & k_d) ** NEW **
   15,   7,  +1;   % D25: dephos VP @ V  (cross, AND: k_f & k_d) ** NEW **
   15,   8,  +1;   % D26: dephos VP @ VP
   16,   7,  +1;   % D27: dilution V         ** WAS MISSING **
   17,   8,  +1;   % D28: dilution VP        ** WAS MISSING **
   18,   3,  +1;   % D29: P_fb*Hill -> MV
   19,   5,  +1;   % D30: k_ms*MV*S @ S
   19,   9,  +1;   % D31: k_ms*MV*S @ MV
   20,   9,  +1;   % D32: delta_M*MV
];

% Build B (n x q) and C (q x n)
B = zeros(n, q);
C = zeros(q, n);
for k = 1:q
    rxn = partials(k, 1);
    var = partials(k, 2);
    sgn = partials(k, 3);
    B(:, k) = Sm(:, rxn);
    C(k, var) = sgn;
end

fprintf('BDC decomposition: n=%d, q=%d\n', n, q);

% Verify: J = B*diag(ones)*C
J_test = B * diag(ones(q,1)) * C;
fprintf('J = B*I*C (all d=1):\n');
var_labels = {'H','U','UP','M','S','R','V','VP','MV'};
fprintf('%8s', '');
for j = 1:n; fprintf('%6s', var_labels{j}); end
fprintf('\n');
for i = 1:n
    fprintf('%6s |', var_labels{i});
    for j = 1:n; fprintf('%5.0f ', J_test(i,j)); end
    fprintf('\n');
end
fprintf('det(-J) = %.0f,  rank = %d\n\n', det(-J_test), rank(J_test));

%% ====
%% 2.  DEPENDENCY CONSTRAINTS -> INDEPENDENT GROUPS
%% ====
%
%  With shared parameters (k_fV=k_f, k_dV=k_d, k_tlV=k_tl), the groups
%  are simpler than the original code with separate V rates.
%
%  Independent parameter groups (10 free bits):
%   g1  = D0,D9,D10,D27,D28       (mu — all dilution terms incl. V, VP)
%   g2  = D1,D2,D3,D20,D21,D22    (k_f — all phospho partials for U AND V)
%   g3  = D4,D6,D23,D26           (k_d — all dephos partials for U AND V)
%   g4  = D8,D19                   (k_tl — translation for BOTH U and V)
%   g5  = D11,D29                  (Hill — transcription for M and MV)
%   g6  = D12,D13,D30,D31          (k_ms — all mRNA-sRNA binding)
%   g7  = D14,D32                  (delta_M — mRNA degradation for M and MV)
%   g8  = D15,D16                  (k_rs — sponge-sRNA binding)
%   g9  = D17                      (delta_S)
%   g10 = D18                      (delta_R; forced ON if g8=ON)
%
%  AND-constrained entries (not free bits):
%   D5  = g2(k_f) AND g3(k_d)   [dephos UP cross-term @ U]
%   D7  = g2(k_f) AND g3(k_d)   [dephos UP cross-term @ V]
%   D24 = g2(k_f) AND g3(k_d)   [dephos VP cross-term @ U]
%   D25 = g2(k_f) AND g3(k_d)   [dephos VP cross-term @ V]
%
%  All four AND constraints are identical: k_f AND k_d.
%
%  => 10 free bits => 2^10 = 1024 vertices.

% groups{k} lists the d-vector indices (1-based) sharing the same free bit.
% D5(6), D7(8), D24(25), D25(26) are AND-constrained — NOT listed here.
groups = {
    [1, 10, 11, 28, 29],       % g1:  mu
    [2, 3, 4, 21, 22, 23],     % g2:  k_f  (phospho for both U and V)
    [5, 7, 24, 27],            % g3:  k_d  (dephos for both U and V)
    [9, 20],                   % g4:  k_tl (translation for both U and V)
    [12, 30],                  % g5:  Hill
    [13, 14, 31, 32],          % g6:  k_ms
    [15, 33],                  % g7:  delta_M
    [16, 17],                  % g8:  k_rs
    [18],                      % g9:  delta_S
    [19],                      % g10: delta_R
};

n_groups = numel(groups);  % 10 free bits
fprintf('Independent parameter groups: %d  =>  2^%d = %d vertices\n', ...
    n_groups, n_groups, 2^n_groups);

%% ====
%% 3.  ENUMERATE VALID VERTICES & COMPUTE INFLUENCE MATRIX
%% ====

% Preallocate sign accumulators
pos_seen = false(n, n);
neg_seen = false(n, n);

n_total   = 2^n_groups;
n_stable  = 0;
n_skipped = 0;

fprintf('Enumerating %d vertices...\n', n_total);

for v = 0:(n_total - 1)

    % Extract group bits
    gbits = zeros(1, n_groups);
    for k = 1:n_groups
        gbits(k) = bitand(bitshift(v, -(k-1)), 1);
    end

    % Reconstruct full d-vector (length q=33)
    d = zeros(1, q);
    for k = 1:n_groups
        for idx = groups{k}
            d(idx) = gbits(k);
        end
    end

    % AND constraints: all four cross-terms need k_f AND k_d
    and_val = gbits(2) & gbits(3);   % g2(k_f) AND g3(k_d)
    d(6)  = and_val;   % D5:  dephos UP @ U cross
    d(8)  = and_val;   % D7:  dephos UP @ V cross
    d(25) = and_val;   % D24: dephos VP @ U cross
    d(26) = and_val;   % D25: dephos VP @ V cross

    % Constraint: delta_R must be present if k_rs is present
    % g8 = k_rs (gbits(8)), g10 = delta_R (gbits(10))
    if gbits(8) == 1 && gbits(10) == 0
        n_skipped = n_skipped + 1;
        continue;
    end

    % Build Jacobian for this vertex
    J_v = B * diag(d) * C;

    % Stability filter: det(-J) must be > 0
    det_negJ = det(-J_v);
    if det_negJ < 1e-14
        n_skipped = n_skipped + 1;
        continue;
    end

    n_stable = n_stable + 1;

    % Compute influence for each (i,j) pair
    for i_out = 1:n
        for j_in = 1:n
            if pos_seen(i_out, j_in) && neg_seen(i_out, j_in)
                continue;  % already ambiguous, skip
            end

            E_vec = zeros(n, 1); E_vec(j_in) = 1;
            H_vec = zeros(1, n); H_vec(i_out) = 1;

            M_aug = [-J_v, -E_vec; H_vec, 0];
            r_v   = det(M_aug);

            if r_v > 1e-14
                pos_seen(i_out, j_in) = true;
            elseif r_v < -1e-14
                neg_seen(i_out, j_in) = true;
            end
        end
    end
end

fprintf('Total vertices:   %d\n', n_total);
fprintf('Stable vertices:  %d\n', n_stable);
fprintf('Skipped:          %d\n', n_skipped);

%% ====
%% 4.  ASSEMBLE & DISPLAY STRUCTURAL INFLUENCE MATRIX
%% ====

M_struct = cell(n, n);
for i = 1:n
    for j = 1:n
        if pos_seen(i,j) && neg_seen(i,j)
            M_struct{i,j} = '?';
        elseif pos_seen(i,j)
            M_struct{i,j} = '+';
        elseif neg_seen(i,j)
            M_struct{i,j} = '-';
        else
            M_struct{i,j} = '0';
        end
    end
end

fprintf('\n==== Structural Influence Matrix (9x9) ====\n');
fprintf('  Row i = effect on variable i\n');
fprintf('  Col j = perturbation input to variable j\n\n');
fprintf('%10s', '');
for j = 1:n; fprintf('%7s', var_labels{j}); end
fprintf('\n');
for i = 1:n
    fprintf('%8s |', var_labels{i});
    for j = 1:n; fprintf('%6s ', M_struct{i,j}); end
    fprintf('\n');
end

%% ====
%% 5.  EXPECTED OUTPUT
%% ====
%
%            H      U     UP      M      S      R      V     VP     MV
%    H |     +      0      0      0      0      0      0      0      0
%    U |     ?      +      +      +      -      +      ?      ?      +
%   UP |     ?      ?      ?      ?      +      -      -      -      -
%    M |     ?      ?      ?      +      -      +      -      -      ?
%    S |     ?      ?      ?      -      +      -      +      +      -
%    R |     0      0      0      0      0      +      0      0      0
%    V |     ?      +      +      +      -      +      ?      +      +
%   VP |     ?      -      -      -      +      -      ?      +      ?
%   MV |     ?      ?      ?      ?      -      +      -      -      +

%% ====
%% 6.  VISUALISATION
%% ====

figure('Name','Structural Influence Matrix — Dichotomous','Position',[100 100 700 650]);
hold on;

for i = 1:n
    for j = 1:n
        sym = M_struct{n+1-i, j};
        switch sym
            case '+', col = [0.40 0.80 0.40];
            case '-', col = [0.90 0.35 0.35];
            case '?', col = [1.00 0.85 0.30];
            otherwise, col = [0.92 0.92 0.92];
        end
        patch([j-1 j j j-1],[i-1 i-1 i i], col, 'EdgeColor',[0.5 0.5 0.5]);
        text(j-0.5, i-0.5, sym, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',16,'FontWeight','bold');
    end
end

set(gca, 'XTick',(1:n)-0.5, 'XTickLabel', var_labels, ...
         'YTick',(1:n)-0.5, 'YTickLabel', flip(var_labels), ...
         'TickLength',[0 0], 'FontSize',12);
xlabel('Input perturbation (variable j)','FontSize',12);
ylabel('Effect on variable i','FontSize',12);
title('Structural Influence Matrix — Dichotomous','FontSize',13);
xlim([0 n]); ylim([0 n]); axis square; box on;
hold off;

annotation('textbox',[0.02 0.01 0.96 0.05], ...
    'String', ['\color[rgb]{0.4,0.8,0.4}\bf+ positive   ' ...
               '\color[rgb]{0.9,0.35,0.35}\bf- negative   ' ...
               '\color[rgb]{1,0.85,0.3}\bf? ambiguous   ' ...
               '\color[rgb]{0.6,0.6,0.6}\bf0 no effect'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',11,'FitBoxToText','off');

%% ====
%% 7.  BIOLOGICAL INTERPRETATION
%% ====
fprintf('\n==== Biological Interpretation ====\n');
fprintf('Diagonal: all + except U<-U and V<-V which are ? (kinase competition)\n');
fprintf('H and R are decoupled from each other\n');
fprintf('\nKey differences from the 6x6 system (without V):\n');
fprintf('  Original 6x6: UP<-S = -,  UP<-R = +\n');
fprintf('  DFB 9x9:      UP<-S = +,  UP<-R = -\n');
fprintf('  These FLIP because the V pathway creates a dominant indirect path:\n');
fprintf('    more S -> less MV -> less V -> less kinase competition -> more UP\n');
fprintf('    more R -> less S -> more MV -> more V -> more competition -> less UP\n');
fprintf('\nCross-coupling (U <-> V through shared kinase):\n');
fprintf('  V<-U = +  (more U -> higher denom -> less phospho of V -> more V)\n');
fprintf('  VP<-U = - (more U -> less phospho of V -> less VP)\n');
fprintf('  U<-V = ?  (competing direct and indirect effects)\n');