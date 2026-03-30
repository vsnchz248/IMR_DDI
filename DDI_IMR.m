function DDI_IMR(use_gt_init, validation_mode)
%DDI_IMR_FIXED  Full joint (R*, V*, S*) DDI via 5N KKT system.
%
% CALLS
%   DDI_IMR_fixed()               production-style R-only inverse
%   DDI_IMR_fixed(false,false)    production-style R-only inverse
%   DDI_IMR_fixed(false,true)     synthetic validation of the same inverse
%
% IMPORTANT
%   This file is a synthetic harness: it loads synthetic_data.mat from
%   IMRv2 and assumes those columns are already nondimensional outputs
%   (dimout = 0). The inverse solve itself only treats R*(t*) as observed.
%   GT curves (R_gt, V_gt, S_gt) are used only for comparison/diagnostics.
%
% THEORY/SETTINGS
%   - fixed gas pressure state pb from observed R only (polytropic law)
%   - wV = 0, beta_R = 0 (campaign settings)
%   - no constitutive GT stress in the solve
%   - native per-experiment grids; no common-grid resampling
%   - S smoothness scaled exactly like beta_S*dt*sum((dS)^2)

if nargin < 1, use_gt_init = false; end
if nargin < 2, validation_mode = false; end
close all; clc;

%% =========================================================================
%  PHYSICAL CONSTANTS
%% =========================================================================
rho_val   = 1064.0;
P_inf_val = 101325.0;
gamma_s   = 0.07;
c_inf     = 1484.0;
c_sp      = c_inf / sqrt(P_inf_val / rho_val);
kappa     = 1.47;
pv_sat    = 0.0;
mu_gt     = 0.1;
G_gt      = 1e4;
Ca_gt     = P_inf_val / G_gt;

%% =========================================================================
%  DATA LOADING
%% =========================================================================
load('synthetic_data.mat');
lambda_vs_Rmax_fig; close all;
fprintf('  Rmax=%.1f um, lambda=%.2f\n', Rmax*1e6, lambda);

Nexp_gen       = numel(synthetic_data);
Rmax_range_gen = linspace(Rmax-0.5*Rmax_std, Rmax+0.5*Rmax_std, Nexp_gen);

Nraw  = size(synthetic_data{1}, 1);
t_raw = zeros(Nraw, Nexp_gen);
R_raw = zeros(Nraw, Nexp_gen);
V_raw = zeros(Nraw, Nexp_gen);
for a = 1:Nexp_gen
    t_raw(:,a) = synthetic_data{a}(:,1);
    R_raw(:,a) = synthetic_data{a}(:,2);
    V_raw(:,a) = synthetic_data{a}(:,3);
end

fprintf('  Assuming synthetic_data columns are already nondimensional IMR outputs (dimout=0).\n');
fprintf('  GT curves from synthetic_data are used only for comparison/diagnostics.\n');

Nexp        = Nexp_gen;
N           = Nraw;
Rmax_actual = Rmax_range_gen;
We_vec      = P_inf_val .* Rmax_actual ./ (2*gamma_s);
Re_gt_vec   = Rmax_actual .* sqrt(rho_val*P_inf_val) ./ mu_gt;

t_grid = zeros(Nexp, N);
for a = 1:Nexp
    t_grid(a,:) = t_raw(:,a)' - t_raw(1,a);
end
dt_vec = t_grid(:,2) - t_grid(:,1);

R_gt  = max(R_raw', 0.02);
V_gt  = V_raw';
R_obs = R_gt;

fprintf('  We: ');  fprintf('%.1f  ', We_vec);   fprintf('\n');
fprintf('  GT-only Re: ');  fprintf('%.1f  ', Re_gt_vec); fprintf('\n');
fprintf('  GT-only Ca: %.4f\n\n', Ca_gt);

%% =========================================================================
%  GT diagnostics only
%% =========================================================================
Req_range_gen = Rmax_range_gen ./ lambda;
Req_star_gt   = Req_range_gen ./ Rmax_actual;
fprintf('  R*_eq used in GT comparison: '); fprintf('%.4f  ', Req_star_gt); fprintf('\n');

S_gt = zeros(Nexp, N);
for a = 1:Nexp
    Rstar     = max(R_gt(a,:), 1e-3);
    ratio     = Req_star_gt(a) ./ Rstar;
    S_gt(a,:) = -1/(2*Ca_gt) .* (5 - ratio.^4 - 4.*ratio) ...
               -4/Re_gt_vec(a) .* V_gt(a,:) ./ Rstar;
end
fprintf('  GT S* range: [%.2f, %.2f]\n\n', min(S_gt(:)), max(S_gt(:)));

%% =========================================================================
%  Native observed grids
%% =========================================================================
fprintf('  Solving on native per-experiment normalized grids (no common-grid resampling).\n\n');

%% =========================================================================
%  R*_eq estimate used in the solve
%% =========================================================================
last10        = max(1, floor(0.1*N));
Req_star_solve = mean(R_obs(:, end-last10+1:end), 2)';
Req_star_solve = max(Req_star_solve, 0.02);
fprintf('  R*_eq used in solve: '); fprintf('%.4f  ', Req_star_solve); fprintf('\n\n');

%% =========================================================================
%  Estimate V* from observed R* only
%% =========================================================================
V_est = zeros(Nexp, N);
for a = 1:Nexp
    R_seed    = prefilter_trace(R_obs(a,:));
    V_est(a,:) = differentiate_fd(R_seed, dt_vec(a));
end
fprintf('  V_est warm start: fixed 5-point binomial prefilter + finite differences.\n\n');

%% =========================================================================
%  pb fixed from observed R only via the polytropic gas law
%% =========================================================================
pb0_vec = zeros(1, Nexp);
pb      = zeros(Nexp, N);
dpb     = zeros(Nexp, N);
for a = 1:Nexp
    Req_a      = Req_star_solve(a);
    gas_coeff  = 1 + 1/(We_vec(a)*Req_a) - pv_sat;
    Rstar       = max(R_obs(a,:), 0.02);
    pb0_vec(a)  = pv_sat + gas_coeff * Req_a^(3*kappa);
    pb(a,:)     = pv_sat + gas_coeff * (Req_a ./ Rstar).^(3*kappa);
    dpb(a,:)    = differentiate_fd(pb(a,:), dt_vec(a));
end
fprintf('  pb0: '); fprintf('%.4f  ', pb0_vec); fprintf('\n\n');

%% =========================================================================
%  Constitutive-free warm start from the native-grid KM closure
%% =========================================================================
S_init = zeros(Nexp, N);
for a = 1:Nexp
    S_init(a,:) = infer_S_warm_start( ...
        N, dt_vec(a), c_sp, We_vec(a), R_obs(a,:), V_est(a,:), pb(a,:), dpb(a,:));
end
fprintf('  S_init (KM-consistent constitutive-free warm start) range: [%.2f, %.2f]\n', ...
    min(S_init(:)), max(S_init(:)));
fprintf('\n\n');

%% =========================================================================
%  Weights and DDI settings
%% =========================================================================
wn_all = ones(Nexp, N);

w_obs_R      = 1.0;
wR           = 1.0;
wV           = 1.0;
wS           = 1.0;
beta_S0      = 1.0;
beta_S1      = 1.0;
beta_R       = 1.0;
anneal_iters = 50;
Npts_total   = Nexp*N;
if validation_mode
    Nbar = Npts_total;
else
    Nbar = min(50*Nexp, Npts_total);
end
max_iter = 200;
tol      = 1e-8;
max_nr   = 50;

fprintf('  Using Nbar = %d centroids (of %d total samples).\n\n', Nbar, Npts_total);
fprintf('  Notes/campaign settings: fixed pb, wV=0, beta_R=0, no extra priors.\n\n');
if validation_mode
    wR_man = 0.0;
    wV_man = 0.0;
    wS_man = 0.0;
    fprintf('  Synthetic validation mode: same R-only inverse, centroid compression disabled.\n');
    fprintf('  GT curves are shown only for comparison.\n\n');
else
    wR_man = wR;
    wV_man = wV;
    wS_man = wS;
    fprintf('  Production-style synthetic harness: same R-only inverse, centroid compression enabled.\n');
    fprintf('  GT curves are shown only for comparison.\n\n');
end

%% =========================================================================
%  Initialization
%% =========================================================================
if use_gt_init
    fprintf('  [DEBUG] Initializing from GT (R_gt, V_gt, S_gt).\n\n');
    R_traj = R_gt;
    V_traj = V_gt;
    S_traj = S_gt;
else
    if validation_mode
        fprintf('  [VALIDATION] Initializing from (R_obs, V_est, S_init).\n\n');
    else
        fprintf('  [PRODUCTION] Initializing from (R_obs, V_est, S_init).\n\n');
    end
    R_traj = R_obs;
    V_traj = V_est;
    S_traj = S_init;
end

fprintf('  Initial RMSE R*: ');
for a = 1:Nexp, fprintf('Exp%d=%.5f  ', a, rms(R_traj(a,:)-R_gt(a,:))); end
fprintf('\n  Initial RMSE V*: ');
for a = 1:Nexp, fprintf('Exp%d=%.5f  ', a, rms(V_traj(a,:)-V_gt(a,:))); end
fprintf('\n  Initial RMSE S*: ');
for a = 1:Nexp, fprintf('Exp%d=%.5f  ', a, rms(S_traj(a,:)-S_gt(a,:))); end
fprintf('\n\n');

%% =========================================================================
%  K-means++ with the same metric as the objective
%% =========================================================================
rng(42);
pts_all = [R_traj(:), V_traj(:), S_traj(:)];
if validation_mode
    Rbar   = pts_all(:,1)';
    Vbar   = pts_all(:,2)';
    Sbar   = pts_all(:,3)';
    assign = reshape(int32(1:Npts_total), Nexp, N);
else
    pts_seed = [sqrt(wR)*R_traj(:), sqrt(wV)*V_traj(:), sqrt(wS)*S_traj(:)];
    ptsn     = pts_seed ./ (std(pts_seed)+eps);
    ci       = zeros(Nbar,1);
    ci(1)    = randi(size(pts_all,1));
    d2       = inf(size(pts_all,1),1);
    for k = 2:Nbar
        d2      = min(d2, sum((ptsn-ptsn(ci(k-1),:)).^2,2));
        mass    = sum(d2);
        if mass <= eps || ~isfinite(mass)
            ci(k) = ci(k-1);
            continue;
        end
        cdf     = cumsum(d2 / mass);
        next_ci = find(rand <= cdf, 1, 'first');
        if isempty(next_ci), next_ci = ci(k-1); end
        ci(k)   = next_ci;
    end
    Rbar   = pts_all(ci,1)';
    Vbar   = pts_all(ci,2)';
    Sbar   = pts_all(ci,3)';
    assign = ones(Nexp, N, 'int32');
end

%% =========================================================================
%  Main loop
%% =========================================================================
fprintf('  %5s  %12s  %10s  %8s  |  %-9s  %-9s  %-9s\n', ...
    'Iter','Objective','Rel.dL','beta_S','RMSE R1','RMSE R2','RMSE R3');
L_hist = nan(max_iter,1);

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

for iter = 1:max_iter
    if iter <= anneal_iters
        beta_S = beta_S0*(beta_S1/beta_S0)^((iter-1)/(anneal_iters-1));
    else
        beta_S = beta_S1;
    end

    if ~validation_mode
        for a = 1:Nexp
            for n = 1:N
                d_assign = wR*(R_traj(a,n)-Rbar).^2 + ...
                           wV*(V_traj(a,n)-Vbar).^2 + ...
                           wS*(S_traj(a,n)-Sbar).^2;
                [~, assign(a,n)] = min(d_assign);
            end
        end
    end

    for a = 1:Nexp
        [R_traj(a,:), V_traj(a,:), S_traj(a,:)] = solve_KKT_fixed( ...
            N, dt_vec(a), c_sp, We_vec(a), R_obs(a,:), pb(a,:), dpb(a,:), ...
            wn_all(a,:)', ...
            w_obs_R, wR_man, wV_man, wS_man, beta_S, Rbar, Vbar, Sbar, ...
            assign(a,:), max_nr, R_traj(a,:)', V_traj(a,:)', S_traj(a,:)');
    end

    if validation_mode
        Rbar = R_traj(:)';
        Vbar = V_traj(:)';
        Sbar = S_traj(:)';
    else
        nR = zeros(1,Nbar);
        nV = zeros(1,Nbar);
        nS = zeros(1,Nbar);
        dn = zeros(1,Nbar);
        for a = 1:Nexp
            for n = 1:N
                k = assign(a,n);
                nR(k) = nR(k) + R_traj(a,n);
                nV(k) = nV(k) + V_traj(a,n);
                nS(k) = nS(k) + S_traj(a,n);
                dn(k) = dn(k) + 1;
            end
        end
        mk = dn > 0;
        Rbar(mk) = nR(mk) ./ dn(mk);
        Vbar(mk) = nV(mk) ./ dn(mk);
        Sbar(mk) = nS(mk) ./ dn(mk);
    end

    L = 0;
    for a = 1:Nexp
        for n = 1:N
            k = assign(a,n);
            L = L + (w_obs_R*(R_traj(a,n)-R_obs(a,n))^2 ...
                   +  wR_man*(R_traj(a,n)-Rbar(k))^2 ...
                   +  wV_man*(V_traj(a,n)-Vbar(k))^2 ...
                   +  wS_man*(S_traj(a,n)-Sbar(k))^2);
        end
        for n = 1:N-1
            L = L + beta_S*dt_vec(a)*(S_traj(a,n+1)-S_traj(a,n))^2;
        end
    end
    L_hist(iter) = L;

    rmse_R = arrayfun(@(a) rms(R_traj(a,:)-R_gt(a,:)), 1:Nexp);
    if iter > 1
        dL = abs(L-L_hist(iter-1)) / (abs(L_hist(iter-1))+eps);
        fprintf('  %5d  %12.4e  %10.2e  %8.3f  |  %-9.5f  %-9.5f  %-9.5f\n', ...
            iter, L, dL, beta_S, rmse_R(1), rmse_R(2), rmse_R(3));
        if dL < tol && iter > anneal_iters
            fprintf('  Converged.\n');
            break;
        end
        if iter >= anneal_iters + 6
            rc = L_hist(iter-4:iter);
            if std(diff(rc)) / (abs(mean(rc))+eps) < 1e-3
                fprintf('  Stagnated.\n');
                break;
            end
        end
    else
        fprintf('  %5d  %12.4e  %10s  %8.3f  |  %-9.5f  %-9.5f  %-9.5f\n', ...
            iter, L, '---', beta_S, rmse_R(1), rmse_R(2), rmse_R(3));
    end
end

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');

%% =========================================================================
%  Plots
%% =========================================================================
cDDI = [0.8 0.15 0.15];
cGT  = [0.05 0.55 0.15];
cols = lines(Nexp);

figure('Name','DDI: R*, V*, S* trajectories','Position',[30 30 1600 320*Nexp]);
for a = 1:Nexp
    subplot(Nexp,3,(a-1)*3+1);
    plot(t_grid(a,:), R_gt(a,:), '-', 'Color', [0 0 0], 'LineWidth', 1.2, 'DisplayName', 'R* GT (synthetic)'); hold on;
    plot(t_grid(a,:), R_traj(a,:), '-', 'Color', cDDI, 'LineWidth', 1.5, 'DisplayName', 'R* DDI');
    yline(Req_star_solve(a), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'R*_{eq} used');
    xlabel('t*'); ylabel('R*'); grid on; box on;
    title(sprintf('Exp %d  R*  We=%.0f  RMSE=%.4f', a, We_vec(a), rms(R_traj(a,:)-R_gt(a,:))));
    if a == 1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+2);
    plot(t_grid(a,:), V_gt(a,:), 'k-', 'LineWidth', 1, 'DisplayName', 'V* GT (synthetic)'); hold on;
    plot(t_grid(a,:), V_traj(a,:), '-', 'Color', cDDI, 'LineWidth', 1.5, 'DisplayName', 'V* DDI');
    xlabel('t*'); ylabel('V*'); grid on; box on;
    title(sprintf('Exp %d  V*  RMSE=%.4f', a, rms(V_traj(a,:)-V_gt(a,:))));
    if a == 1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+3);
    plot(t_grid(a,:), S_gt(a,:), '-', 'Color', cGT, 'LineWidth', 2.5, 'DisplayName', 'S* GT (stress integral)'); hold on;
    plot(t_grid(a,:), S_traj(a,:), '-', 'Color', cDDI, 'LineWidth', 1.5, 'DisplayName', 'S* DDI');
    xlabel('t*'); ylabel('S*'); grid on; box on;
    ylim([min(S_gt(a,:))*1.5-1, max(S_gt(a,:))*1.5+1]);
    title(sprintf('Exp %d  S*  RMSE=%.3f', a, rms(S_traj(a,:)-S_gt(a,:))));
    if a == 1, legend('Location','best','FontSize',7); end
end

figure('Name','DDI: Phase space, manifold, convergence','Position',[100 50 1300 800]);

subplot(2,2,1);
for a = 1:Nexp
    lc = cols(a,:);
    plot(R_gt(a,:), S_gt(a,:), '-', 'Color', 0.5*lc+0.4, 'LineWidth', 2, ...
        'DisplayName', sprintf('GT exp%d', a)); hold on;
    plot(R_traj(a,:), S_traj(a,:), '--', 'Color', lc, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('DDI exp%d', a));
end
xlabel('R*'); ylabel('S*'); title('S*(R*): GT (solid) vs DDI (dashed)');
ylim([min(S_gt(:))*1.5-1, max(S_gt(:))*1.5+1]);
grid on; legend('Location','best','FontSize',7); box on;

subplot(2,2,2);
for a = 1:Nexp
    lc = cols(a,:);
    plot3(R_traj(a,:), V_traj(a,:), S_traj(a,:), '-', 'Color', lc, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('DDI exp%d', a)); hold on;
end
occupied = find(arrayfun(@(k) any(assign(:)==k), 1:Nbar));
scatter3(Rbar(occupied), Vbar(occupied), Sbar(occupied), 20, 'r', 'filled', ...
    'DisplayName', 'Centroids');
xlabel('R*'); ylabel('V*'); zlabel('S*');
title('3D Phase space (R*, V*, S*) with centroids');
legend('Location','best','FontSize',7); grid on; view(35,20); box on;

subplot(2,2,3);
v = ~isnan(L_hist);
semilogy(find(v), L_hist(v), '-o', 'Color', [0.2 0.45 0.75], 'LineWidth', 1.5, ...
    'MarkerSize', 3, 'MarkerFaceColor', [0.2 0.45 0.75]);
xlabel('Iteration'); ylabel('L'); title('Convergence'); grid on;

subplot(2,2,4);
for a = 1:Nexp
    lc = cols(a,:);
    plot(t_grid(a,:), S_gt(a,:), '-', 'Color', 0.5*lc+0.4, 'LineWidth', 2, ...
        'DisplayName', sprintf('GT exp%d', a)); hold on;
    plot(t_grid(a,:), S_traj(a,:), '--', 'Color', lc, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('DDI exp%d', a));
end
xlabel('t*'); ylabel('S*'); title('S*(t*): GT vs DDI');
ylim([min(S_gt(:))*1.5-1, max(S_gt(:))*1.5+1]);
grid on; legend('Location','best','FontSize',7); box on;

fprintf('\n  %-5s  %-8s  %-12s  %-12s  %-12s\n', 'Exp', 'We', 'RMSE R*', 'RMSE V*', 'RMSE S*');
for a = 1:Nexp
    fprintf('  %-5d  %-8.1f  %-12.5f  %-12.5f  %-12.4f\n', a, ...
        We_vec(a), ...
        rms(R_traj(a,:)-R_gt(a,:)), ...
        rms(V_traj(a,:)-V_gt(a,:)), ...
        rms(S_traj(a,:)-S_gt(a,:)));
end

fprintf('\n  DDI discrete constraint residuals on each native experiment grid\n');
fprintf('  %-5s  %-12s  %-12s\n', 'Exp', 'DDI KM rms', 'DDI kin rms');
for a = 1:Nexp
    [km_ddi, kin_ddi] = constraint_residuals(N, dt_vec(a), c_sp, We_vec(a), R_traj(a,:), V_traj(a,:), S_traj(a,:), pb(a,:), dpb(a,:));
    fprintf('  %-5d  %-12.4e  %-12.4e\n', a, rms(km_ddi), rms(kin_ddi));
end
end

function [km_res,kin_res] = constraint_residuals(N,dt,c,We,R,V,S,pb,dpb)
km_res  = zeros(1, N);
kin_res = zeros(1, N);

kin_res(1) = V(1);
for n = 2:N-1
    dV = (V(n+1)-V(n-1)) / (2*dt);
    dS = (S(n+1)-S(n-1)) / (2*dt);
    km_res(n) = (1-V(n)/c)*R(n)*dV + 1.5*(1-V(n)/(3*c))*V(n)^2 ...
              - (1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
              - R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
    kin_res(n) = V(n) - (R(n+1)-R(n-1)) / (2*dt);
end
kin_res(N) = V(N) - (R(N)-R(N-1)) / dt;
end

%% =========================================================================
function res = kkt_residual_only( ...
    N,dt,c,We,R_obs,pb,dpb,wn, ...
    w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,sol)

iR = @(n) n;
iV = @(n) N+n;
iS = @(n) 2*N+n;
iL = @(n) 3*N+n;
iM = @(n) 4*N+n;

R   = max(sol(iR(1):iR(N)),0.005);
V   = sol(iV(1):iV(N));
S   = sol(iS(1):iS(N));
lam = sol(iL(1):iL(N));
mu  = sol(iM(1):iM(N));

res = zeros(5*N,1);

res(iR(1)) = R(1)-1;
res(iV(1)) = V(1);
res(iL(1)) = lam(1);
res(iL(N)) = lam(N);
res(iM(1)) = mu(1);
res(iM(N)) = V(N) - (R(N)-R(N-1)) / dt;

for n = 2:N-1
    kk = assign(n);
    dV = (V(n+1)-V(n-1)) / (2*dt);
    dS = (S(n+1)-S(n-1)) / (2*dt);

    A = (1-V(n)/c)*dV - 1/(We*R(n)^2) - (dpb(n)+dS)/c;
    B = -R(n)/c*dV + 3*V(n) - 3*V(n)^2/(2*c) - (pb(n)+S(n)-1)/c;

    Fc = (1-V(n)/c)*R(n)*dV + 1.5*(1-V(n)/(3*c))*V(n)^2 ...
       - (1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
       - R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
    res(iL(n)) = Fc;

    Gc = V(n) - (R(n+1)-R(n-1)) / (2*dt);
    res(iM(n)) = Gc;

    if n < N-1
        mu_term = (mu(n+1)-mu(n-1)) / (2*dt);
    else
        mu_term = mu(N)/dt - mu(N-2)/(2*dt);
    end

    res(iR(n)) = 2*wn(n)*(w_obs_R*(R(n)-R_obs(n)) + wR*(R(n)-Rbar(kk))) ...
               + A*lam(n) + mu_term;

    cn1 = -(1-V(n+1)/c)*R(n+1)/(2*dt);
    cm1 =  (1-V(n-1)/c)*R(n-1)/(2*dt);
    res(iV(n)) = 2*wn(n)*wV*(V(n)-Vbar(kk)) ...
               + B*lam(n) + cn1*lam(n+1) + cm1*lam(n-1) + mu(n);

    lS = -(1+V(n)/c)*lam(n);
    if n+1 <= N-1, lS = lS + R(n+1)*lam(n+1)/(2*c*dt); end
    if n-1 >= 2,   lS = lS - R(n-1)*lam(n-1)/(2*c*dt); end

    sm = 2*beta_S*dt*((S(n)-S(n-1)) - (S(n+1)-S(n)));
    res(iS(n)) = 2*wn(n)*wS*(S(n)-Sbar(kk)) + sm + lS;
end

kN = assign(N);
k1 = assign(1);

res(iR(N)) = 2*wn(N)*(w_obs_R*(R(N)-R_obs(N)) + wR*(R(N)-Rbar(kN))) ...
           - mu(N-1)/(2*dt) - mu(N)/dt;

cN = (1-V(N-1)/c)*R(N-1)/(2*dt);
res(iV(N)) = 2*wn(N)*wV*(V(N)-Vbar(kN)) + cN*lam(N-1) + mu(N);

sm1 = 2*beta_S*dt*(S(1)-S(2));
res(iS(1)) = 2*wn(1)*wS*(S(1)-Sbar(k1)) + sm1 + R(2)*lam(2)/(2*c*dt);

smN = 2*beta_S*dt*(S(N)-S(N-1));
res(iS(N)) = 2*wn(N)*wS*(S(N)-Sbar(kN)) + smN - R(N-1)*lam(N-1)/(2*c*dt);
end

%% =========================================================================
function [R_out,V_out,S_out] = solve_KKT_fixed( ...
    N,dt,c,We,R_obs,pb,dpb,wn, ...
    w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,max_nr, ...
    R_w,V_w,S_w)

iR = @(n) n;
iV = @(n) N+n;
iS = @(n) 2*N+n;
iL = @(n) 3*N+n;
iM = @(n) 4*N+n;
feas_idx = [(3*N+2):(4*N-1), (4*N+1):(5*N)];

sol = [max(R_w,0.02); V_w; S_w; zeros(2*N,1)];
sol(iR(1)) = 1;
sol(iV(1)) = 0;

for nr = 1:max_nr
    R   = max(sol(iR(1):iR(N)),0.005);
    V   = sol(iV(1):iV(N));
    S   = sol(iS(1):iS(N));
    lam = sol(iL(1):iL(N));
    mu  = sol(iM(1):iM(N));

    res = zeros(5*N,1);
    nz  = 0;
    TI  = zeros(72*N,1);
    TJ  = zeros(72*N,1);
    TS  = zeros(72*N,1);

    res(iR(1)) = R(1)-1;
    nz=nz+1; TI(nz)=iR(1); TJ(nz)=iR(1); TS(nz)=1;

    res(iV(1)) = V(1);
    nz=nz+1; TI(nz)=iV(1); TJ(nz)=iV(1); TS(nz)=1;

    res(iL(1)) = lam(1);
    nz=nz+1; TI(nz)=iL(1); TJ(nz)=iL(1); TS(nz)=1;

    res(iL(N)) = lam(N);
    nz=nz+1; TI(nz)=iL(N); TJ(nz)=iL(N); TS(nz)=1;

    res(iM(1)) = mu(1);
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iM(1); TS(nz)=1;

    res(iM(N)) = V(N) - (R(N)-R(N-1)) / dt;
    nz=nz+1; TI(nz)=iM(N); TJ(nz)=iV(N);   TS(nz)=1;
    nz=nz+1; TI(nz)=iM(N); TJ(nz)=iR(N);   TS(nz)=-1/dt;
    nz=nz+1; TI(nz)=iM(N); TJ(nz)=iR(N-1); TS(nz)=1/dt;

    for n = 2:N-1
        kk = assign(n);
        dV = (V(n+1)-V(n-1)) / (2*dt);
        dS = (S(n+1)-S(n-1)) / (2*dt);

        A = (1-V(n)/c)*dV - 1/(We*R(n)^2) - (dpb(n)+dS)/c;
        B = -R(n)/c*dV + 3*V(n) - 3*V(n)^2/(2*c) - (pb(n)+S(n)-1)/c;

        Fc = (1-V(n)/c)*R(n)*dV + 1.5*(1-V(n)/(3*c))*V(n)^2 ...
           - (1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
           - R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
        res(iL(n)) = Fc;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iR(n);   TS(nz)=A;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n);   TS(nz)=B;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n+1); TS(nz)=(1-V(n)/c)*R(n)/(2*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n-1); TS(nz)=-(1-V(n)/c)*R(n)/(2*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n);   TS(nz)=-(1+V(n)/c);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n+1); TS(nz)=-R(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n-1); TS(nz)=R(n)/(2*c*dt);

        Gc = V(n) - (R(n+1)-R(n-1)) / (2*dt);
        res(iM(n)) = Gc;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n);   TS(nz)=1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n+1); TS(nz)=-1/(2*dt);
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n-1); TS(nz)=1/(2*dt);

        if n < N-1
            mu_term = (mu(n+1)-mu(n-1)) / (2*dt);
            mu_plus_coeff = 1 / (2*dt);
        else
            mu_term = mu(N)/dt - mu(N-2)/(2*dt);
            mu_plus_coeff = 1 / dt;
        end

        res(iR(n)) = 2*wn(n)*(w_obs_R*(R(n)-R_obs(n)) + wR*(R(n)-Rbar(kk))) ...
                   + A*lam(n) + mu_term;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iR(n);   TS(nz)=2*wn(n)*(w_obs_R+wR)+2*lam(n)/(We*R(n)^3);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n);   TS(nz)=-lam(n)*dV/c;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n+1); TS(nz)=lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n-1); TS(nz)=-lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n+1); TS(nz)=-lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n-1); TS(nz)=lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iL(n);   TS(nz)=A;
        if n < N-1
            nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n+1); TS(nz)=mu_plus_coeff;
        else
            nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(N); TS(nz)=mu_plus_coeff;
        end
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n-1); TS(nz)=-1/(2*dt);

        cn1 = -(1-V(n+1)/c)*R(n+1)/(2*dt);
        cm1 =  (1-V(n-1)/c)*R(n-1)/(2*dt);
        res(iV(n)) = 2*wn(n)*wV*(V(n)-Vbar(kk)) ...
                   + B*lam(n) + cn1*lam(n+1) + cm1*lam(n-1) + mu(n);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n);   TS(nz)=-lam(n)*dV/c;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n+1); TS(nz)=-(1-V(n+1)/c)*lam(n+1)/(2*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n-1); TS(nz)=(1-V(n-1)/c)*lam(n-1)/(2*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n);   TS(nz)=2*wn(n)*wV + (3-3*V(n)/c)*lam(n);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n+1); TS(nz)=-R(n)*lam(n)/(2*c*dt)+R(n+1)*lam(n+1)/(2*c*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n-1); TS(nz)=R(n)*lam(n)/(2*c*dt)-R(n-1)*lam(n-1)/(2*c*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iS(n);   TS(nz)=-lam(n)/c;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n);   TS(nz)=B;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n+1); TS(nz)=cn1;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n-1); TS(nz)=cm1;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n);   TS(nz)=1;

        lS = -(1+V(n)/c)*lam(n);
        if n+1 <= N-1, lS = lS + R(n+1)*lam(n+1)/(2*c*dt); end
        if n-1 >= 2,   lS = lS - R(n-1)*lam(n-1)/(2*c*dt); end

        sm = 2*beta_S*dt*((S(n)-S(n-1)) - (S(n+1)-S(n)));
        res(iS(n)) = 2*wn(n)*wS*(S(n)-Sbar(kk)) + sm + lS;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n);   TS(nz)=2*wn(n)*wS + 4*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n+1); TS(nz)=-2*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n-1); TS(nz)=-2*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iV(n);   TS(nz)=-lam(n)/c;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n);   TS(nz)=-(1+V(n)/c);
        if n+1 <= N-1
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iR(n+1); TS(nz)=lam(n+1)/(2*c*dt);
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n+1); TS(nz)=R(n+1)/(2*c*dt);
        end
        if n-1 >= 2
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iR(n-1); TS(nz)=-lam(n-1)/(2*c*dt);
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n-1); TS(nz)=-R(n-1)/(2*c*dt);
        end
    end

    kN = assign(N);
    k1 = assign(1);

    res(iR(N)) = 2*wn(N)*(w_obs_R*(R(N)-R_obs(N)) + wR*(R(N)-Rbar(kN))) ...
               - mu(N-1)/(2*dt) - mu(N)/dt;
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iR(N);   TS(nz)=2*wn(N)*(w_obs_R+wR);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iM(N-1); TS(nz)=-1/(2*dt);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iM(N);   TS(nz)=-1/dt;

    cN = (1-V(N-1)/c)*R(N-1)/(2*dt);
    res(iV(N)) = 2*wn(N)*wV*(V(N)-Vbar(kN)) + cN*lam(N-1) + mu(N);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N);   TS(nz)=2*wn(N)*wV + 1e-10;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iR(N-1); TS(nz)=(1-V(N-1)/c)*lam(N-1)/(2*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N-1); TS(nz)=-R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iL(N-1); TS(nz)=cN;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iM(N);   TS(nz)=1;

    sm1 = 2*beta_S*dt*(S(1)-S(2));
    res(iS(1)) = 2*wn(1)*wS*(S(1)-Sbar(k1)) + sm1 + R(2)*lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(1); TS(nz)=2*wn(1)*wS + 2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(2); TS(nz)=-2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iR(2); TS(nz)=lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iL(2); TS(nz)=R(2)/(2*c*dt);

    smN = 2*beta_S*dt*(S(N)-S(N-1));
    res(iS(N)) = 2*wn(N)*wS*(S(N)-Sbar(kN)) + smN - R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N);   TS(nz)=2*wn(N)*wS + 2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N-1); TS(nz)=-2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iR(N-1); TS(nz)=-lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iL(N-1); TS(nz)=-R(N-1)/(2*c*dt);

    J = sparse(TI(1:nz), TJ(1:nz), TS(1:nz), 5*N, 5*N);
    reg = 1e-6 * max(max(abs(diag(J))), 1);
    delta = -(J + reg*speye(5*N)) \ res;

    rn       = norm(res, inf);
    full_nrm = norm(res, 2);
    feas_nrm = norm(res(feas_idx), 2);
    alpha    = 1;
    accepted = false;
    for ls = 1:20
        trial = sol + alpha*delta;
        if all(trial(1:N) > 0.02) && all(trial(1:N) < 2)
            trial(iR(1)) = 1;
            trial(iV(1)) = 0;
            trial(iL(1)) = 0;
            trial(iL(N)) = 0;
            trial(iM(1)) = 0;
            trial(1:N)   = min(max(trial(1:N),0.005),2);

            trial_res = kkt_residual_only( ...
                N,dt,c,We,R_obs,pb,dpb,wn, ...
                w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,trial);
            trial_full = norm(trial_res,2);
            trial_feas = norm(trial_res(feas_idx),2);
            if (trial_feas < feas_nrm) || ...
               (trial_feas <= feas_nrm && trial_full < full_nrm)
                sol = trial;
                rn = norm(trial_res,inf);
                accepted = true;
                break;
            end
        end
        alpha = alpha * 0.5;
    end
    if ~accepted, break; end
    if rn < 1e-8, break; end
end

R_out = sol(iR(1):iR(N))';
V_out = sol(iV(1):iV(N))';
S_out = sol(iS(1):iS(N))';
end

%% =========================================================================
function dy = differentiate_fd(y, dt)
%DIFFERENTIATE_FD  Fixed finite-difference derivative on a uniform grid.
%
% Uses fourth-order accurate one-sided/centered formulas when enough
% samples are available, and falls back to lower-order formulas otherwise.

y = y(:).';
N = numel(y);
dy = zeros(1, N);

if N == 1
    return;
elseif N == 2
    dy(:) = (y(2) - y(1)) / dt;
    return;
elseif N == 3
    dy(1) = (-3*y(1) + 4*y(2) - y(3)) / (2*dt);
    dy(2) = (y(3) - y(1)) / (2*dt);
    dy(3) = (3*y(3) - 4*y(2) + y(1)) / (2*dt);
    return;
elseif N == 4
    dy(1) = (-11*y(1) + 18*y(2) - 9*y(3) + 2*y(4)) / (6*dt);
    dy(2) = (-2*y(1) - 3*y(2) + 6*y(3) - y(4)) / (6*dt);
    dy(3) = (y(1) - 6*y(2) + 3*y(3) + 2*y(4)) / (6*dt);
    dy(4) = (11*y(4) - 18*y(3) + 9*y(2) - 2*y(1)) / (6*dt);
    return;
end

dy(1)   = (-25*y(1) + 48*y(2) - 36*y(3) + 16*y(4) - 3*y(5)) / (12*dt);
dy(2)   = (-3*y(1) - 10*y(2) + 18*y(3) - 6*y(4) + y(5)) / (12*dt);
dy(N-1) = (3*y(N) + 10*y(N-1) - 18*y(N-2) + 6*y(N-3) - y(N-4)) / (12*dt);
dy(N)   = (25*y(N) - 48*y(N-1) + 36*y(N-2) - 16*y(N-3) + 3*y(N-4)) / (12*dt);

for n = 3:N-2
    dy(n) = (y(n-2) - 8*y(n-1) + 8*y(n+1) - y(n+2)) / (12*dt);
end
end

%% =========================================================================
function ys = prefilter_trace(y)
%PREFILTER_TRACE  Fixed local smoothing for warm-start differentiation.
%
% Uses a single five-point binomial pass. This is a numerical warm-start
% stabilizer only, not a model parameter or constitutive prior.

y = y(:).';
N = numel(y);
ys = y;

if N <= 4
    return;
end

w = [1 4 6 4 1] / 16;

ys(1)   = y(1);
ys(2)   = (5*y(1) + 10*y(2) + y(3)) / 16;
ys(N-1) = (y(N-2) + 10*y(N-1) + 5*y(N)) / 16;
ys(N)   = y(N);

for n = 3:N-2
    ys(n) = w(1)*y(n-2) + w(2)*y(n-1) + w(3)*y(n) + w(4)*y(n+1) + w(5)*y(n+2);
end
end

%% =========================================================================
function S_init = infer_S_warm_start(N, dt, c, We, R, V, pb, dpb)
%INFER_S_WARM_START  Constitutive-free S warm start from the KM closure.
%
% This uses only observed R, numerically inferred V, and fixed pb(R) to
% solve the linear part of the native-grid KM equation for S.

R = max(R(:).', 0.02);
V = V(:).';
pb = pb(:).';
dpb = dpb(:).';

dV = differentiate_fd(V, dt);

q = (1 - V/c) .* R .* dV ...
  + 1.5 * (1 - V/(3*c)) .* V.^2 ...
  - (1 + V/c) .* (pb - 1 ./ (We * R) - 1) ...
  - (R / c) .* (dpb + V ./ (We * R.^2));

A = zeros(N, N);

if N == 1
    A(1,1) = 1 + V(1)/c;
elseif N == 2
    A(1,1) = 1 + V(1)/c - R(1)/(c*dt);
    A(1,2) = R(1)/(c*dt);
    A(2,1) = -R(2)/(c*dt);
    A(2,2) = 1 + V(2)/c + R(2)/(c*dt);
else
    A(1,1) = 1 + V(1)/c - R(1)/(c*dt);
    A(1,2) = R(1)/(c*dt);

    for n = 2:N-1
        A(n,n-1) = -R(n)/(2*c*dt);
        A(n,n)   = 1 + V(n)/c;
        A(n,n+1) = R(n)/(2*c*dt);
    end

    A(N,N-1) = -R(N)/(c*dt);
    A(N,N)   = 1 + V(N)/c + R(N)/(c*dt);
end

S_init = (A \ q(:)).';
end
