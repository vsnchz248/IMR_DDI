function DDI_example()
%DDI_KINEMATICS  Data-driven identification with noisy kinematic measurements
%
% Implements the staggered minimisation of Ulloa's working notes.
%
% Objective function:
%   L = sum_alpha sum_n w_n [ f^2(z_n) + d^2(z_n, z_bar_k) ]
%     + beta_F * sum_alpha sum_{n=1}^{N-1} w_n * dt * (F_{n+1} - F_n)^2
%
% Force law (state-dependent, shared across all experiments):
%   F(x,v) = -k3*x^3 - cd*v*|v|
%   Both terms are strictly restoring/dissipative -- ODE is unconditionally
%   stable for any initial condition.
%
% Observations:
%   - x is measured directly with noise
%   - v is NOT observed; estimated by central differences on Gaussian-smoothed x_obs
%
% Force panel shows two estimates:
%   - Trivial (algebraic): F = m*a + c*v + k*x  (requires known m,c,k)
%   - DDI: recovered without assuming a force model

clear; close all; clc;
rng(42);   % controls measurement noise only

%% =========================================================================
%  USER PARAMETERS
%% =========================================================================

% -- Structural parameters (used in ODE and KKT -- DDI requires these) --
m = 1.0;
c = 0.1;
k = 4.0;

% -- Nonlinear force law (used only for data generation) --
% -k3*x^3 : hardening cubic spring  (always restoring)
% -cd*v|v| : quadratic drag          (always dissipative)
k3 = 2.0;
cd = 0.05;
F_law = @(x,v) -k3*x.^3 - cd*v.*abs(v);

% -- Time grid --
T     = 8.0;
T_buf = 2.0;
T_ext = T + T_buf;
N_ext = 126;
t_ext = linspace(0, T_ext, N_ext);
dt    = t_ext(2) - t_ext(1);
N_out = round(T/T_ext * (N_ext-1)) + 1;
t_out = t_ext(1:N_out);

% -- Experiments -- change only Nexp; ICs are drawn automatically --
Nexp = 3;
rng(13);                         % IC seed (independent of noise seed above)
x_ini = 0.5*randn(Nexp, 1);
v_ini = 0.3*randn(Nexp, 1);

% -- Noise (applied to x only) --
noise_frac = 0.1;

% -- DDI weights --
w_obs_x = 50.0;   % trust x measurements strongly
w_obs_v =  0.5;   % v estimated numerically -- trust weakly
wx = 1.0;
wv = 1.0;
wF = 3.0;
wn = ones(N_ext, 1);

% -- beta_F annealing --
beta_F0      = 80.0;
beta_F1      = 15.0;
anneal_iters = 50;

% -- Centroids and solver --
Nbar     = 25 * Nexp;
max_iter = 200;
tol      = 1e-9;

%% =========================================================================
%  GENERATE SYNTHETIC DATA
%% =========================================================================
fprintf('=== Generating synthetic data ===\n');

x_true_ext = zeros(Nexp, N_ext);
v_true_ext = zeros(Nexp, N_ext);

ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
for a = 1:Nexp
    rhs_f = @(~,y) [y(2); (F_law(y(1),y(2)) - c*y(2) - k*y(1))/m];
    [t_sol, Y] = ode45(rhs_f, t_ext, [x_ini(a); v_ini(a)], ode_opts);
    if length(t_sol) ~= N_ext
        error(['ODE failed for experiment %d at t=%.4f. ' ...
               'Reduce IC amplitudes or check F_law stability.'], a, t_sol(end));
    end
    x_true_ext(a,:) = Y(:,1)';
    v_true_ext(a,:) = Y(:,2)';
end

sig   = std(x_true_ext(:,1:N_out), [], 'all');
x_obs = x_true_ext + noise_frac*sig*randn(Nexp, N_ext);
% v is NOT observed -- constructed below from x_obs

%% =========================================================================
%  GAUSSIAN SMOOTHING
%  Used for: (1) estimating v_obs, (2) initialising trajectories
%% =========================================================================
smooth_sigma = 3;
hw = ceil(3*smooth_sigma);
gw = exp(-(-hw:hw).^2 / (2*smooth_sigma^2));
gw = gw / sum(gw);

x_smooth = zeros(Nexp, N_ext);
for a = 1:Nexp
    x_smooth(a,:) = conv(x_obs(a,:), gw, 'same');
    x_smooth(a,1) = x_ini(a);   % enforce known IC exactly
end

%% =========================================================================
%  VELOCITY ESTIMATE FROM SMOOTHED DISPLACEMENT
%  Central differences on x_smooth -- only velocity information available.
%% =========================================================================
v_obs = zeros(Nexp, N_ext);
for a = 1:Nexp
    v_obs(a,2:N_ext-1) = (x_smooth(a,3:N_ext) - x_smooth(a,1:N_ext-2)) / (2*dt);
    v_obs(a,1)         = v_ini(a);
    v_obs(a,N_ext)     = v_obs(a,N_ext-1);
end

%% =========================================================================
%  INITIALISE TRAJECTORIES
%% =========================================================================
fprintf('=== Initialising ===\n');

x_traj = x_obs;
v_traj = v_obs;
F_traj = zeros(Nexp, N_ext);
for a = 1:Nexp
    for n = 2:N_ext-1
        F_traj(a,n) = m*(v_traj(a,n+1)-v_traj(a,n-1))/(2*dt) ...
                    + c*v_traj(a,n) + k*x_smooth(a,n);
    end
    F_traj(a,1)     = F_traj(a,2);
    F_traj(a,N_ext) = F_traj(a,N_ext-1);
end

%% =========================================================================
%  K-MEANS++ CENTROID INITIALISATION
%% =========================================================================
fprintf('=== K-means++ centroid seeding ===\n');

all_pts   = [x_traj(:), v_traj(:), F_traj(:)];
pts_scale = std(all_pts) + eps;
pts_n     = all_pts ./ pts_scale;

Npts  = size(all_pts, 1);
ci    = zeros(Nbar, 1);
ci(1) = randi(Npts);
min_d2 = inf(Npts, 1);
for kk = 2:Nbar
    d2c    = sum((pts_n - pts_n(ci(kk-1),:)).^2, 2);
    min_d2 = min(min_d2, d2c);
    prob   = min_d2 / sum(min_d2);
    ci(kk) = find(rand < cumsum(prob), 1);
end
xbar = all_pts(ci, 1)';
vbar = all_pts(ci, 2)';
Fbar = all_pts(ci, 3)';

assign = ones(Nexp, N_ext, 'int32');

%% =========================================================================
%  DDI MAIN LOOP
%% =========================================================================
fprintf('=== DDI alternating minimisation ===\n');
fprintf('  %5s  %14s  %14s  %10s\n','Iter','Objective L','Rel. dL','beta_F');

L_hist = nan(max_iter, 1);

for iter = 1:max_iter

    % beta_F geometric annealing
    if iter <= anneal_iters
        frac   = (iter-1) / (anneal_iters-1);
        beta_F = beta_F0 * (beta_F1/beta_F0)^frac;
    else
        beta_F = beta_F1;
    end

    % Assignment: nearest centroid in weighted (x,v,F) space
    pts_all = [x_traj(:), v_traj(:), F_traj(:)];
    cen_all = [xbar(:),   vbar(:),   Fbar(:)  ];
    D2 = wx*(pts_all(:,1) - cen_all(:,1)').^2 + ...
         wv*(pts_all(:,2) - cen_all(:,2)').^2 + ...
         wF*(pts_all(:,3) - cen_all(:,3)').^2;
    [~, assign_vec] = min(D2, [], 2);
    assign = reshape(int32(assign_vec), Nexp, N_ext);

    % Step 1: KKT solve per experiment (Eqs. 12-17)
    for a = 1:Nexp
        [x_traj(a,:), v_traj(a,:), F_traj(a,:)] = ...
            solve_trajectory(m, c, k, dt, N_ext, wn, ...
                             w_obs_x, w_obs_v, wx, wv, wF, beta_F, ...
                             x_obs(a,:), v_obs(a,:), ...
                             x_ini(a), v_ini(a), ...
                             xbar, vbar, Fbar, assign(a,:));
    end

    % Step 2: Centroid update (Eqs. 19-21)
    num_x = zeros(1, Nbar);
    num_v = zeros(1, Nbar);
    num_F = zeros(1, Nbar);
    den   = zeros(1, Nbar);
    for a = 1:Nexp
        for n = 1:N_ext
            kk        = assign(a,n);
            w         = wn(n);
            num_x(kk) = num_x(kk) + w*x_traj(a,n);
            num_v(kk) = num_v(kk) + w*v_traj(a,n);
            num_F(kk) = num_F(kk) + w*F_traj(a,n);
            den(kk)   = den(kk)   + w;
        end
    end
    for kk = 1:Nbar
        if den(kk) > 0
            xbar(kk) = num_x(kk)/den(kk);
            vbar(kk) = num_v(kk)/den(kk);
            Fbar(kk) = num_F(kk)/den(kk);
        end
    end

    % Convergence check
    L = eval_objective(Nexp, N_ext, wn, dt, ...
                       w_obs_x, w_obs_v, wx, wv, wF, beta_F, ...
                       x_traj, v_traj, F_traj, x_obs, v_obs, ...
                       xbar, vbar, Fbar, assign);
    L_hist(iter) = L;

    if iter > 1
        dL = abs(L - L_hist(iter-1)) / (abs(L_hist(iter-1)) + eps);
        fprintf('  %5d  %14.6e  %14.6e  %10.3f\n', iter, L, dL, beta_F);
        if dL < tol && iter > anneal_iters
            fprintf('  --> Converged at iteration %d.\n', iter);
            break;
        end
    else
        fprintf('  %5d  %14.6e  %14s  %10.3f\n', iter, L, '---', beta_F);
    end
end

%% =========================================================================
%  TRIM TO [0,T]
%% =========================================================================
fprintf('=== Plotting ===\n');

x_traj_out = x_traj(:,    1:N_out);
v_traj_out = v_traj(:,    1:N_out);
F_traj_out = F_traj(:,    1:N_out);
x_obs_out  = x_obs(:,     1:N_out);
x_true_out = x_true_ext(:,1:N_out);
v_true_out = v_true_ext(:,1:N_out);

% Trivial (algebraic) force: F = m*a + c*v + k*x using DDI trajectories.
% Requires known m, c, k. Serves as a baseline to compare DDI against.
F_trivial = zeros(Nexp, N_out);
for a = 1:Nexp
    for n = 2:N_out-1
        F_trivial(a,n) = m*(v_traj_out(a,n+1)-v_traj_out(a,n-1))/(2*dt) ...
                       + c*v_traj_out(a,n) + k*x_traj_out(a,n);
    end
    F_trivial(a,1)     = F_trivial(a,2);
    F_trivial(a,N_out) = F_trivial(a,N_out-1);
end

%% =========================================================================
%  FIGURE 1: Trajectories
%% =========================================================================
cBlue = [0.20 0.45 0.75];
cRed  = [0.80 0.15 0.15];

figure('Name','DDI - Trajectories','NumberTitle','off', ...
       'Position',[50 50 1300 260*Nexp]);
for a = 1:Nexp

    % -- Displacement --
    subplot(Nexp, 3, (a-1)*3 + 1);
    plot(t_out, x_obs_out(a,:),  'ko', 'MarkerSize',4); hold on;
    plot(t_out, x_true_out(a,:), 'k-', 'LineWidth',1.5);
    plot(t_out, x_traj_out(a,:), '--', 'Color',cRed, 'LineWidth',1.5);
    xlabel('t'); ylabel('x'); grid on; box on;
    title(sprintf('Exp %d - Displacement',a));
    if a==1, legend('Observed','True','DDI','Location','best'); end

    % -- Velocity --
    subplot(Nexp, 3, (a-1)*3 + 2);
    plot(t_out, v_true_out(a,:), 'k-', 'LineWidth',1.5); hold on;
    plot(t_out, v_traj_out(a,:), '--', 'Color',cRed, 'LineWidth',1.5);
    xlabel('t'); ylabel('v'); grid on; box on;
    title(sprintf('Exp %d - Velocity',a));
    if a==1, legend('True','DDI','Location','best'); end

    % -- Force: trivial algebraic baseline vs DDI --
    subplot(Nexp, 3, (a-1)*3 + 3);
    plot(t_out, F_trivial(a,:),  'k-', 'LineWidth',1.5); hold on;
    plot(t_out, F_traj_out(a,:), '--', 'Color',cRed, 'LineWidth',1.5);
    xlabel('t'); ylabel('F'); grid on; box on;
    title(sprintf('Exp %d - Force',a));
    if a==1, legend('Trivial case','DDI','Location','best'); end
end

%% =========================================================================
%  FIGURE 2: Phase space, centroids, convergence
%% =========================================================================
figure('Name','DDI - Phase Space and Convergence','NumberTitle','off', ...
       'Position',[100 100 950 720]);

subplot(2,1,1);
cols = lines(Nexp);
for a = 1:Nexp
    plot3(x_traj_out(a,:), v_traj_out(a,:), F_traj_out(a,:), '-', ...
          'Color',cols(a,:), 'LineWidth',1, 'DisplayName',sprintf('Traj %d',a));
    hold on;
end
scatter3(xbar, vbar, Fbar, 60, 'r', 'filled', 'DisplayName','Centroids');
xlabel('x'); ylabel('v'); zlabel('F');
title('Phase-space trajectories and identified centroids');
legend('Location','best'); grid on; view(30,20);

subplot(2,1,2);
valid = ~isnan(L_hist);
semilogy(find(valid), L_hist(valid), '-o', 'Color',cBlue, ...
         'LineWidth',1.5, 'MarkerSize',4, 'MarkerFaceColor',cBlue);
xlabel('Iteration'); ylabel('Objective  L');
title('DDI convergence'); grid on; box on;

%% =========================================================================
%  ERROR SUMMARY (x and v only -- F has no unambiguous ground truth)
%% =========================================================================
fprintf('\n=== Recovery error (window [0,T]) ===\n');
fprintf('  %6s  %10s  %10s\n','Exp','RMSE x','RMSE v');
for a = 1:Nexp
    fprintf('  %6d  %10.4f  %10.4f\n', a, ...
        rms(x_traj_out(a,:) - x_true_out(a,:)), ...
        rms(v_traj_out(a,:) - v_true_out(a,:)));
end
fprintf('=== Done ===\n');
end


%% =========================================================================
%  LOCAL FUNCTION: solve_trajectory
%% =========================================================================
function [x, v, F] = solve_trajectory(m, c, k, dt, N, wn, ...
    w_obs_x, w_obs_v, wx, wv, wF, beta_F, ...
    x_obs, v_obs, x_ini, v_ini, ...
    xbar, vbar, Fbar, assign)

col_x   = @(n) n;
col_v   = @(n) N   + n;
col_F   = @(n) 2*N + n;
col_lam = @(n) 3*N + (n-1);
col_mu  = @(n) 4*N - 2 + (n-1);

Ndof   = 5*N - 4;
nz_max = 30*N;
TI  = zeros(nz_max,1);
TJ  = zeros(nz_max,1);
TS  = zeros(nz_max,1);
rhs = zeros(Ndof,1);
nz  = 0;
row = 0;

% Eq. (17): initial conditions
row=row+1; nz=nz+1; TI(nz)=row; TJ(nz)=col_x(1); TS(nz)=1; rhs(row)=x_ini;
row=row+1; nz=nz+1; TI(nz)=row; TJ(nz)=col_v(1); TS(nz)=1; rhs(row)=v_ini;

% Interior nodes n = 2,...,N-1
for n = 2:N-1
    kk = assign(n);

    % Eq. (12): discrete dynamics
    row=row+1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_v(n+1); TS(nz)= m/(2*dt);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_v(n-1); TS(nz)=-m/(2*dt);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_v(n);   TS(nz)= c;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_x(n);   TS(nz)= k;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_F(n);   TS(nz)=-1;

    % Eq. (13): kinematics
    row=row+1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_v(n);   TS(nz)= 1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_x(n+1); TS(nz)=-1/(2*dt);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_x(n-1); TS(nz)= 1/(2*dt);

    % Eq. (14): stationarity w.r.t. x_n
    row=row+1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_x(n);   TS(nz)= 2*wn(n)*(w_obs_x+wx);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(n); TS(nz)= k;
    if n<=N-2, nz=nz+1; TI(nz)=row; TJ(nz)=col_mu(n+1); TS(nz)= 1/(2*dt); end
    if n>=3,   nz=nz+1; TI(nz)=row; TJ(nz)=col_mu(n-1); TS(nz)=-1/(2*dt); end
    rhs(row) = 2*wn(n)*(w_obs_x*x_obs(n) + wx*xbar(kk));

    % Eq. (15): stationarity w.r.t. v_n
    row=row+1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_v(n);   TS(nz)= 2*wn(n)*(w_obs_v+wv);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(n); TS(nz)= c;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_mu(n);  TS(nz)= 1;
    if n>=3,   nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(n-1); TS(nz)= m/(2*dt); end
    if n<=N-2, nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(n+1); TS(nz)=-m/(2*dt); end
    rhs(row) = 2*wn(n)*(w_obs_v*v_obs(n) + wv*vbar(kk));

    % Eq. (16): stationarity w.r.t. F_n (extended with beta_F smoothing)
    diag_F = 2*wn(n)*wF + 2*beta_F*dt*(wn(n-1)+wn(n));
    row=row+1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_F(n);   TS(nz)= diag_F;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(n); TS(nz)=-1;
    nz=nz+1; TI(nz)=row; TJ(nz)=col_F(n-1); TS(nz)=-2*beta_F*dt*wn(n-1);
    nz=nz+1; TI(nz)=row; TJ(nz)=col_F(n+1); TS(nz)=-2*beta_F*dt*wn(n);
    rhs(row) = 2*wn(n)*wF*Fbar(kk);
end

% Boundary stationarity at n=1 and n=N
% dL/dF(1)
row=row+1;
nz=nz+1; TI(nz)=row; TJ(nz)=col_F(1); TS(nz)= 2*wn(1)*wF + 2*beta_F*dt*wn(1);
nz=nz+1; TI(nz)=row; TJ(nz)=col_F(2); TS(nz)=-2*beta_F*dt*wn(1);
rhs(row) = 2*wn(1)*wF*Fbar(assign(1));

% dL/dx(N)
row=row+1;
nz=nz+1; TI(nz)=row; TJ(nz)=col_x(N);    TS(nz)= 2*wn(N)*(w_obs_x+wx);
nz=nz+1; TI(nz)=row; TJ(nz)=col_mu(N-1); TS(nz)=-1/(2*dt);
rhs(row) = 2*wn(N)*(w_obs_x*x_obs(N) + wx*xbar(assign(N)));

% dL/dv(N)
row=row+1;
nz=nz+1; TI(nz)=row; TJ(nz)=col_v(N);     TS(nz)= 2*wn(N)*(w_obs_v+wv);
nz=nz+1; TI(nz)=row; TJ(nz)=col_lam(N-1); TS(nz)= m/(2*dt);
rhs(row) = 2*wn(N)*(w_obs_v*v_obs(N) + wv*vbar(assign(N)));

% dL/dF(N)
row=row+1;
nz=nz+1; TI(nz)=row; TJ(nz)=col_F(N);   TS(nz)= 2*wn(N)*wF + 2*beta_F*dt*wn(N-1);
nz=nz+1; TI(nz)=row; TJ(nz)=col_F(N-1); TS(nz)=-2*beta_F*dt*wn(N-1);
rhs(row) = 2*wn(N)*wF*Fbar(assign(N));

if row ~= Ndof
    error('KKT assembly error: %d rows assembled, %d expected.', row, Ndof);
end

sol = sparse(TI(1:nz), TJ(1:nz), TS(1:nz), Ndof, Ndof) \ rhs;
x = sol(1    :  N)';
v = sol(N+1  :2*N)';
F = sol(2*N+1:3*N)';
end


%% =========================================================================
%  LOCAL FUNCTION: eval_objective
%% =========================================================================
function L = eval_objective(Nexp, N, wn, dt, ...
    w_obs_x, w_obs_v, wx, wv, wF, beta_F, ...
    x_traj, v_traj, F_traj, x_obs, v_obs, ...
    xbar, vbar, Fbar, assign)
L = 0;
for a = 1:Nexp
    for n = 1:N
        kk = assign(a,n);
        f2 = w_obs_x*(x_traj(a,n)-x_obs(a,n))^2 ...
           + w_obs_v*(v_traj(a,n)-v_obs(a,n))^2;
        d2 = wx*(x_traj(a,n)-xbar(kk))^2 ...
           + wv*(v_traj(a,n)-vbar(kk))^2 ...
           + wF*(F_traj(a,n)-Fbar(kk))^2;
        L = L + wn(n)*(f2+d2);
    end
    for n = 1:N-1
        L = L + beta_F*wn(n)*dt*(F_traj(a,n+1)-F_traj(a,n))^2;
    end
end
end