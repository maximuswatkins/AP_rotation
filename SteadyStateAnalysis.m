%% =========================================================
%  LINEARISATION — compute Jacobian entries at steady state
%% =========================================================

% RUN THIS ONLY AFTER RUNNING THE ORIGINAL MODEL!

% Extract steady-state values (last time point)
U_ss    = U(end);
U_P_ss  = U_P(end);
M_ss    = M(end);
S_ss    = S(end);
R_ss    = R_sum(end);
p.HK_tot  = HK_sum(end);

% --- Shared quantities ---
Delta_ss  = p.k_ap + p.k_f * U_ss + p.mu;
D_B_ss    = (p.k_Dplus * U_P_ss * p.D_sum) / (p.k_Dminus + p.k_Dplus * U_P_ss);
phi_ss    = (p.k_rs * p.delta_R * R_ss) / (p.delta_R + p.k_rs * S_ss)^2;

% --- d[U]/dt entries ---
% Partial w.r.t. U
a11 = - p.k_f * p.HK_tot * (p.k_ap / Delta_ss) ...
      + p.k_d * p.HK_tot * (p.k_f * U_P_ss / Delta_ss) ...
      - p.k_d * p.HK_tot * ((p.k_f * U_ss + p.mu) * p.k_f * U_P_ss) / Delta_ss^2 ...
      - p.mu;

% Partial w.r.t. U_P
a12 = p.k_d * p.HK_tot * (p.k_f * U_ss + p.mu) / Delta_ss;

% Partial w.r.t. M
a13 = p.k_tl;

% --- d[U_P]/dt entries ---
% Partial w.r.t. U  (note: a21 = -(a11 + mu), included explicitly for clarity)
a21 = p.k_f * p.HK_tot * (p.k_ap / Delta_ss) ...
      - p.k_d * p.HK_tot * (p.k_f * U_P_ss / Delta_ss) ...
      + p.k_d * p.HK_tot * ((p.k_f * U_ss + p.mu) * p.k_f * U_P_ss) / Delta_ss^2;

% Partial w.r.t. U_P (includes linearised decoy QSS term)
a22 = - p.k_d * p.HK_tot * (p.k_f * U_ss + p.mu) / Delta_ss ...
      - p.mu ...
      - p.k_Dplus * (p.D_sum - D_B_ss) ...
      + (p.k_Dplus * p.k_Dminus * p.D_sum) / (p.k_Dminus + p.k_Dplus * U_P_ss)^2;

% --- d[M]/dt entries ---
% Partial w.r.t. U_P (linearised Hill function)
H_prime = p.k_tx * p.G_0 * p.n * p.K_tx^p.n * U_P_ss^(p.n-1) ...
          / (p.K_tx^p.n + U_P_ss^p.n)^2;
a32 = H_prime;

% Partial w.r.t. M
a33 = -(p.k_ms * S_ss + p.k_tl + p.delta_M);

% Partial w.r.t. S
a34 = -p.k_ms * M_ss;

% --- d[S]/dt entries ---
% Partial w.r.t. M
a43 = -p.k_ms * S_ss;

% Partial w.r.t. S (includes linearised sponge term)
a44 = -(p.k_ms * M_ss + p.delta_S + phi_ss);

% --- d[R_sum]/dt ---
% Only entry: partial w.r.t. R_sum
a55 = -p.delta_R;

%% =========================================================
%  ASSEMBLE STATE SPACE AND COMPUTE TRANSFER FUNCTIONS
%% =========================================================
A_mat = [a11, a12, a13,   0,    0;
         a21, a22,   0,   0,    0;
           0, a32, a33, a34,    0;
           0,   0, a43, a44,    0;
           0,   0,   0,   0,  a55];

B_mat = [0, 0;
         0, 0;
         0, 0;
         1, 0;   % beta_S input
         0, 1];  % beta_r input

C_mat = [1, 0, 0, 0, 0];  % output: U
D_mat = [0, 0];

sys       = ss(A_mat, B_mat, C_mat, D_mat);
tf_betaS  = tf(sys(1,1));   % G(s): beta_S -> U
tf_betar  = tf(sys(1,2));   % G(s): beta_r -> U (expect ~0)

fprintf('\n=== Transfer function: beta_S -> U ===\n'); tf_betaS
fprintf('\n=== Transfer function: beta_r -> U ===\n'); tf_betar

% Bode plot
figure;
%bode(tf_betaS);
margin(tf_betaS);
title('Transfer function: \beta_S \rightarrow U');
grid on;

% Sanity check: real parts of eigenvalues should all be negative (stable SS)
fprintf('\n=== Eigenvalues of A (stability check) ===\n');
disp(eig(A_mat));