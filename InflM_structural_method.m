%% structural_influence.m
% Computes the structural influence matrix using the BDC decomposition
% and vertex algorithm (Procedure 1) from Giordano et al. (2016).
%
% Accounts for parameter dependencies in the phosphorylation/
% dephosphorylation cycle and restricts to stable vertices.
%
% State vector: x = [h; u; p; m; s; r]  (n = 6)

clear; clc;

%% ===== Build BDC Decomposition =====
% 12 reactions, 17 distinct partial derivatives (q = 17)

n = 6;   % number of state variables
q = 17;  % number of free parameters D_0 ... D_16

% Stoichiometric matrix S (n x 12)
%   Columns: reactions
%   Rows: species [h, u, p, m, s, r]
S = zeros(n, 12);
S(:, 1)  = [-1;  0;  0;  0;  0;  0];  % Rxn 1:  dilution of h
S(:, 2)  = [ 0; -1; +1;  0;  0;  0];  % Rxn 2:  phosphorylation
S(:, 3)  = [ 0; +1; -1;  0;  0;  0];  % Rxn 3:  dephosphorylation
S(:, 4)  = [ 0; +1;  0;  0;  0;  0];  % Rxn 4:  translation
S(:, 5)  = [ 0; -1;  0;  0;  0;  0];  % Rxn 5:  dilution of u
S(:, 6)  = [ 0;  0; -1;  0;  0;  0];  % Rxn 6:  dilution of p
S(:, 7)  = [ 0;  0;  0; +1;  0;  0];  % Rxn 7:  transcription (Hill)
S(:, 8)  = [ 0;  0;  0; -1; -1;  0];  % Rxn 8:  mRNA-sRNA binding
S(:, 9)  = [ 0;  0;  0; -1;  0;  0];  % Rxn 9:  mRNA degradation
S(:, 10) = [ 0;  0;  0;  0; -1;  0];  % Rxn 10: sponge-sRNA binding
S(:, 11) = [ 0;  0;  0;  0; -1;  0];  % Rxn 11: sRNA degradation
S(:, 12) = [ 0;  0;  0;  0;  0; -1];  % Rxn 12: sponge degradation

% Partial derivatives: (reaction_index, variable_index, sign)
%   D0:  d(rxn1)/dh   = mu                          > 0
%   D1:  d(rxn2)/dh   = k_f * A(u) * u              > 0
%   D2:  d(rxn2)/du   = k_f * h * ...               > 0
%   D3:  d(rxn3)/dh   = k_d * B(u) * p              > 0
%   D4:  d(rxn3)/du   = k_d * h * p * ...           > 0
%   D5:  d(rxn3)/dp   = k_d * h * B                 > 0
%   D6:  d(rxn4)/dm   = k_tl                        > 0
%   D7:  d(rxn5)/du   = mu                          > 0
%   D8:  d(rxn6)/dp   = mu                          > 0
%   D9:  d(rxn7)/dp   = k_tx * G0 * dHill/dp        > 0
%   D10: d(rxn8)/dm   = k_ms * s                    > 0
%   D11: d(rxn8)/ds   = k_ms * m                    > 0
%   D12: d(rxn9)/dm   = delta_M                     > 0
%   D13: d(rxn10)/ds  = k_rs * sponge_frac * r      > 0
%   D14: d(rxn10)/dr  = k_rs * sponge_frac_r * s    > 0
%   D15: d(rxn11)/ds  = delta_S                     > 0
%   D16: d(rxn12)/dr  = delta_R                     > 0

partials = [
%   rxn  var  sign
     1,   1,  +1;   % D0
     2,   1,  +1;   % D1
     2,   2,  +1;   % D2
     3,   1,  +1;   % D3
     3,   2,  +1;   % D4
     3,   3,  +1;   % D5
     4,   4,  +1;   % D6
     5,   2,  +1;   % D7
     6,   3,  +1;   % D8
     7,   3,  +1;   % D9
     8,   4,  +1;   % D10
     8,   5,  +1;   % D11
     9,   4,  +1;   % D12
    10,   5,  +1;   % D13
    10,   6,  +1;   % D14
    11,   5,  +1;   % D15
    12,   6,  +1;   % D16
];

% Build B (n x q) and C (q x n)
B = zeros(n, q);
C = zeros(q, n);
for k = 1:q
    rxn = partials(k, 1);
    var = partials(k, 2);
    sgn = partials(k, 3);
    B(:, k) = S(:, rxn);
    C(k, var) = sgn;
end

fprintf('B matrix (n x q = %d x %d):\n', n, q);
disp(B);
fprintf('C matrix (q x n = %d x %d):\n', q, n);
disp(C);

% Verify: J = B * diag(ones) * C should give correct sign pattern
J_test = B * diag(ones(q,1)) * C;
fprintf('J = B*I*C (all params = 1):\n');
disp(J_test);

%% ===== Vertex Algorithm with Dependency Constraints =====
% Enumerate all 2^17 = 131072 vertices of {0,1}^17
% Apply dependency constraints to filter valid vertices
% Restrict to stable vertices (det(-J) > 0)

fprintf('Enumerating vertices...\n');

valid_vertices = [];
n_total = 2^q;

for v = 0:(n_total - 1)
    % Extract bits
    d = zeros(1, q);
    for k = 1:q
        d(k) = bitand(bitshift(v, -(k-1)), 1);
    end

    % ---- Dependency constraints ----
    % (indices are 1-based: D1=d(1), D2=d(2), ..., D17=d(17))

    % 1) mu constraint: D0 = D7 = D8  =>  d(1) = d(8) = d(9)
    if ~(d(1) == d(8) && d(1) == d(9))
        continue;
    end

    % 2) k_f constraint: D1 <-> D2  =>  d(2) = d(3)
    if d(2) ~= d(3)
        continue;
    end

    % 3) k_d constraint: D3 <-> D5  =>  d(4) = d(6)
    if d(4) ~= d(6)
        continue;
    end

    % 4) D4 needs both k_f > 0 AND k_d > 0:
    %    d(5) = 1 => d(2) = 1 AND d(4) = 1
    %    d(2) = 1 AND d(4) = 1 => d(5) = 1
    if d(5) == 1 && (d(2) == 0 || d(4) == 0)
        continue;
    end
    if d(2) == 1 && d(4) == 1 && d(5) == 0
        continue;
    end

    % 5) k_ms constraint: D10 <-> D11  =>  d(11) = d(12)
    if d(11) ~= d(12)
        continue;
    end

    % 6) k_rs constraint: D13 <-> D14  =>  d(14) = d(15)
    if d(14) ~= d(15)
        continue;
    end

    % 7) delta_R constraint: D13=1 or D14=1 => D16=1
    %    d(14)=1 or d(15)=1 => d(17)=1
    if (d(14) == 1 || d(15) == 1) && d(17) == 0
        continue;
    end

    valid_vertices = [valid_vertices; d]; %#ok<AGROW>
end

n_valid = size(valid_vertices, 1);
fprintf('Total vertices: %d\n', n_total);
fprintf('Valid vertices (satisfying constraints): %d\n', n_valid);

%% ===== Compute Structural Influence Matrix =====
var_labels = {'h', 'u', 'p', 'm', 's', 'r'};
M_struct = cell(n, n);

n_stable = 0;
n_unstable = 0;
stability_counted = false;

for i_out = 1:n
    for j_in = 1:n
        % E = e_{j_in},  H = e_{i_out}^T
        E_vec = zeros(n, 1);
        E_vec(j_in) = 1;
        H_vec = zeros(1, n);
        H_vec(i_out) = 1;

        pos_seen = false;
        neg_seen = false;

        for v_idx = 1:n_valid
            d = valid_vertices(v_idx, :);
            D_mat = diag(d);
            J_v = B * D_mat * C;

            det_negJ = det(-J_v);

            % Skip unstable / degenerate vertices
            if det_negJ < 1e-14
                if ~stability_counted
                    n_unstable = n_unstable + 1;
                end
                continue;
            end

            if ~stability_counted
                n_stable = n_stable + 1;
            end

            % Augmented matrix: [[-J, -E]; [H, 0]]
            M_aug = [-J_v, -E_vec; H_vec, 0];
            r_v = det(M_aug);

            if abs(r_v) < 1e-14
                s_v = 0;
            else
                s_v = sign(r_v);
            end

            if s_v > 0
                pos_seen = true;
            elseif s_v < 0
                neg_seen = true;
            end

            % Early termination
            if pos_seen && neg_seen
                break;
            end
        end

        stability_counted = true;

        if pos_seen && neg_seen
            M_struct{i_out, j_in} = '?';
        elseif pos_seen
            M_struct{i_out, j_in} = '+';
        elseif neg_seen
            M_struct{i_out, j_in} = '-';
        else
            M_struct{i_out, j_in} = '0';
        end
    end
end

%% ===== Display Results =====
fprintf('\nStable vertices: %d\n', n_stable);
fprintf('Unstable/degenerate vertices: %d\n', n_unstable);

fprintf('\n===== Structural Influence Matrix M =====\n');
fprintf('  (Row i = effect on variable i, Column j = input to variable j)\n\n');
fprintf('%8s', '');
for j = 1:n
    fprintf('%8s', var_labels{j});
end
fprintf('\n');
for i = 1:n
    fprintf('%8s', var_labels{i});
    for j = 1:n
        fprintf('%8s', M_struct{i, j});
    end
    fprintf('\n');
end

%% ===== Biological Interpretation =====
fprintf('\n===== Biological Interpretation =====\n');
fprintf('Diagonal: all + (sanity check: input to x_i increases x_i)\n');
fprintf('h and r are decoupled (independent ODEs)\n');
fprintf('u, p, m form a mutually positive influence group\n');
fprintf('s (sRNA) negatively influences u, p, m and vice versa\n');
fprintf('r (sponge) positively influences u, p, m; negatively influences s\n');
fprintf('Column h has ? entries: h catalyzes both phosphorylation\n');
fprintf('  and dephosphorylation, so net effect is parameter-dependent\n');

% %% ===== Expected Output =====
% fprintf('\n===== Expected Structural Influence Matrix =====\n');
% fprintf('       h     u     p     m     s     r\n');
% fprintf('  h    +     0     0     0     0     0\n');
% fprintf('  u    ?     +     +     +     -     +\n');
% fprintf('  p    ?     +     +     +     -     +\n');
% fprintf('  m    ?     +     +     +     -     +\n');
% fprintf('  s    ?     -     -     -     +     -\n');
% fprintf('  r    0     0     0     0     0     +\n');

%% visualize_influence.m
% Visualizes the structural influence matrix as a colored grid.
% Run structural_influence.m first to generate M_struct.

figure('Name', 'Structural Influence Matrix', 'Position', [100 100 500 450]);

n = 6;
var_labels = {'HK','U','U_P','M','S','R'};

% Map symbols to numeric values for coloring
%   +1 = positive (green), -1 = negative (red), 0 = zero (white), NaN = ambiguous (yellow)
num_map = zeros(n, n);
for i = 1:n
    for j = 1:n
        switch M_struct{i, j}
            case '+', num_map(i, j) = 1;
            case '-', num_map(i, j) = -1;
            case '0', num_map(i, j) = 0;
            case '?', num_map(i, j) = 0.5;  % placeholder for ambiguous
        end
    end
end

% Flip vertically so row 1 is at the top
num_plot = flipud(num_map);

% Draw colored patches
hold on;
for i = 1:n
    for j = 1:n
        val = num_plot(i, j);
        if val == 1
            col = [0 158 115]./255;      % green
        elseif val == -1
            col = [230 159 0]./255;    % orange
        elseif val == 0.5
            col = [240 228 66]./255;     % yellow (ambiguous)
        else
            col = [86 180 230]./255;   % blue (zero)
        end
        patch([j-1 j j j-1], [i-1 i-1 i i], col, 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
    end
end

% Overlay sign symbols
for i = 1:n
    for j = 1:n
        row_flipped = n + 1 - i;  % map back to original row
        txt = M_struct{row_flipped, j};
        text(j - 0.5, i - 0.5, txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 18, 'FontWeight', 'bold');
    end
end

% Axis labels
set(gca, 'XTick', (1:n) - 0.5, 'XTickLabel', var_labels, ...
         'YTick', (1:n) - 0.5, 'YTickLabel', flip(var_labels), ...
         'TickLength', [0 0], 'FontSize', 13);
xlabel('Input perturbation to variable j', 'FontSize', 13);
ylabel('Effect on variable i', 'FontSize', 13);
title('A. Structural Influence Matrix (Non-Dichotomous)', 'FontSize', 15);

xlim([0 n]); ylim([0 n]);
axis square;
box on;
hold off;

% Add legend
% annotation('textbox', [0.02 0.01 0.96 0.05], ...
%     'String', '\color[rgb]{0.4,0.8,0.4}\bf+ positive   \color[rgb]{0.9,0.35,0.35}\bf- negative   \color[rgb]{1,0.85,0.3}\bf? ambiguous   \color[rgb]{0.6,0.6,0.6}\bf0 no effect (P.A.)', ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
%     'FontSize', 11, 'FitBoxToText', 'off');

