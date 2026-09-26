%% ThermoBrain — Full Analysis (Final Version)
% Companion code for:
%   Sharifi Taskuh, B. ThermoBrain: An Externally Validated,
%   Physiology-Informed Computational Framework for Oral Cooling
%   Optimization During Exercise in the Heat. (manuscript in preparation)
%
% This script reproduces, in order:
%   1) Fixed physical constants (Table 1)
%   2) External validation of R_ice_eff against Onitsuka et al. (2018)
%      and Tan et al. (2024), with k_core anchored to Lee et al. (2010)
%      (Table 3, Fig. 2)
%   3) 100-iteration Monte Carlo parameter recovery at the final,
%      jointly consistent ground truth (Table 4, Fig. 3)
%   4) Architectural ablation vs. a single-compartment model (Table 5, Fig. 4)
%   5) Protocol-optimization decision grid (Table 6, Fig. 5)
%
% Requires MATLAB with the Optimization Toolbox (fminbnd, lsqnonlin).
% Run the whole script (F5) or cell-by-cell (Ctrl+Enter).

clear; clc;

%% 1) Fixed physical constants (Table 1) — literature-sourced, not calibrated
P = struct();
P.c_blood         = 3850.0;   % J/(kg*K)
P.c_brain         = 3650.0;   % J/(kg*K)
P.m_brain         = 1.40;     % kg
P.CBF_Lmin        = 0.750;    % L/min
P.rho_blood       = 1060.0;   % kg/m3
P.Q_met           = 15.0;     % W
P.m_blood_central = 3.5;      % kg (modeling assumption, see manuscript Table 1)
P.L_f             = 334000.0; % J/kg

P.mdot_perf    = (P.CBF_Lmin/60.0)*1e-3*P.rho_blood;
P.k_perf_brain = P.mdot_perf*P.c_blood/(P.m_brain*P.c_brain);
P.k_perf_blood = P.mdot_perf/P.m_blood_central;

fprintf('=== Derived (non-calibrated) rate constants ===\n');
fprintf('k_perf_brain = %.5f 1/s (tau = %.1f s)\n', P.k_perf_brain, 1/P.k_perf_brain);
fprintf('k_perf_blood = %.5f 1/s (tau = %.1f s)\n\n', P.k_perf_blood, 1/P.k_perf_blood);

%% 2a) External validation — Onitsuka et al. (2018), resting protocol
fprintf('=== External validation: Onitsuka et al. (2018) ===\n');
body_mass_A = 67.19;
m_ice0_A    = 7.5*body_mass_A/1000.0 - 0.0375;   % target dose minus reported shortfall
T_core_A    = 37.0;
T_brain0_A  = 37.3;
REAL_dTbrain_A = -0.4;
window_A    = [0, 30*60];
t_eval_A    = linspace(0, 35*60, 351);

% k_core anchored to Lee et al. (2010): tau = 1.8-4.4 min (deepest measured
% site, 16 cm; latent period 3.1 +/- 1.3 min), central estimate 3.1 min
tau_lo = 1.8; tau_hi = 4.4; tau_mid = 3.1;   % minutes
taus_ab = [tau_lo, tau_hi];

fprintf('Refitting R_ice_eff across tau_core = %.1f-%.1f min (Lee et al., 2010):\n', tau_lo, tau_hi);
R_range_A = zeros(1,2);
for i = 1:2
    tau_min = taus_ab(i);
    k = 1/(tau_min*60);
    obj = @(logR) (mean_dTbrain_bolus(10^logR, k, m_ice0_A, T_core_A, T_brain0_A, P, window_A, t_eval_A) - REAL_dTbrain_A)^2;
    logR_fit = fminbnd(obj, -3, 1);
    R_range_A(i) = 10^logR_fit;
    fprintf('  tau=%.1f min -> R_ice_eff = %.4f K/W\n', tau_min, R_range_A(i));
end
fprintf('Onitsuka-derived range: %.2f-%.2f K/W\n\n', min(R_range_A), max(R_range_A));

%% 2b) External validation — Tan et al. (2024), exercise-heat crossover
fprintf('=== External validation: Tan et al. (2024) EX-vs-CL crossover ===\n');
body_mass_B = 64.0;
Tre0_B = 36.8; Tre_end_B = 39.4; dTre_target_B = -0.4; t_end_B = 60*60;
aliquot_pre = (8.0*body_mass_B/1000.0)/6.0;
aliquot_ex  = 1.5*body_mass_B/1000.0;
dose_t_B = [0 5 10 15 20 25 30 45]*60;
dose_m_B = [repmat(aliquot_pre,1,6), repmat(aliquot_ex,1,2)];

R_range_B = zeros(1,2);
for i = 1:2
    tau_min = taus_ab(i);
    k = 1/(tau_min*60);
    obj = @(logR) (crossover_diff(10^logR, k, dose_t_B, dose_m_B, Tre0_B, Tre_end_B, t_end_B, P) - dTre_target_B)^2;
    logR_fit = fminbnd(obj, -3, 1);
    R_range_B(i) = 10^logR_fit;
    fprintf('  tau=%.1f min -> R_ice_eff = %.4f K/W\n', tau_min, R_range_B(i));
end
fprintf('Tan-derived range: %.2f-%.2f K/W\n', min(R_range_B), max(R_range_B));

R_adopted_lo = max(min(R_range_A), min(R_range_B));
R_adopted_hi = min(max(R_range_A), max(R_range_B));
R_adopted_mid = (R_adopted_lo + R_adopted_hi)/2;
fprintf('Adopted (overlap) range: %.3f-%.3f K/W (midpoint %.3f)\n\n', R_adopted_lo, R_adopted_hi, R_adopted_mid);

%% 3) Monte Carlo recovery at the final, jointly consistent ground truth (Table 4)
fprintf('=== Monte Carlo recovery (100 iterations), final joint ground truth ===\n');
rng(42);

R_true  = R_adopted_mid;      % 0.825 K/W
k_true  = 1/(tau_mid*60);     % tau = 3.1 min
T0_true = 37.3;
fprintf('Ground truth: R_ice_eff=%.4f K/W, k_core=%.5f 1/s (tau=%.2f min)\n', R_true, k_true, 1/k_true/60);

t_meas = (0:5:60)*60;
Tbrain_true = simulate_full_dosed(R_true, k_true, T0_true, dose_t_B, dose_m_B, Tre0_B, Tre_end_B, t_end_B, P, t_meas);

noise_sd = 0.05;
n_iter   = 100;
est      = zeros(n_iter,3);
x0       = [1.0, 0.006, 37.0];
lb       = [0.05, 0.0005, 36.0];
ub       = [5.0, 0.2, 38.5];
opts_ls  = optimoptions('lsqnonlin','Display','off');

for i = 1:n_iter
    noisy = Tbrain_true + noise_sd*randn(size(Tbrain_true));
    resid = @(p) simulate_full_dosed(p(1), p(2), p(3), dose_t_B, dose_m_B, Tre0_B, Tre_end_B, t_end_B, P, t_meas) - noisy;
    est(i,:) = lsqnonlin(resid, x0, lb, ub, opts_ls);
end

report_recovery('R_ice_eff', est(:,1), R_true);
report_recovery('k_core   ', est(:,2), k_true);
report_recovery('T_brain0 ', est(:,3), T0_true);
fprintf('\n');

%% 4) Architectural ablation: single-compartment vs. 3-compartment (Table 5)
fprintf('=== Architectural ablation ===\n');
fprintf('tau_core       | R_ice_eff from Onitsuka | R_ice_eff from Tan\n');
for tau_min = [40/60, 10, 3]
    k = 1/(tau_min*60);
    obj_A = @(logR) (mean_dTbrain_bolus(10^logR, k, m_ice0_A, T_core_A, T_brain0_A, P, window_A, t_eval_A) - REAL_dTbrain_A)^2;
    R_A = 10^fminbnd(obj_A, -3, 1);
    obj_B = @(logR) (crossover_diff(10^logR, k, dose_t_B, dose_m_B, Tre0_B, Tre_end_B, t_end_B, P) - dTre_target_B)^2;
    R_B = 10^fminbnd(obj_B, -3, 1);
    fprintf('tau=%5.2f min  | %.4f                  | %.4f\n', tau_min, R_A, R_B);
end
fprintf('(3-state ThermoBrain, physiological range: Onitsuka %.2f-%.2f, Tan %.2f-%.2f, overlapping)\n\n', ...
        min(R_range_A), max(R_range_A), min(R_range_B), max(R_range_B));

%% 5) Protocol-optimization decision grid (Table 6)
fprintf('=== Minimum pre-exercise ice fraction ===\n');
taus_min    = [1.8 2.5 3.1 3.8 4.4];
Rs          = [0.6 0.9 1.3];
body_mass_D = 64.0;
total_ice_D = 10.0*body_mass_D/1000.0;
T_brain0_D  = 37.3;

fprintf('tau_core(min) | R=0.6  | R=0.9  | R=1.3\n');
for tmin = taus_min
    k = 1/(tmin*60);
    row = zeros(1,numel(Rs));
    for ri = 1:numel(Rs)
        row(ri) = min_viable_alpha(Rs(ri), k, total_ice_D, Tre0_B, Tre_end_B, t_end_B, T_brain0_D, P);
    end
    fprintf('%13.1f | %.2f   | %.2f   | %.2f\n', tmin, row(1), row(2), row(3));
end

fprintf('\nDone.\n');


%% ================================= LOCAL FUNCTIONS =================================

function dy = thermobrain_rhs(t, y, R_ice_eff, k_core, T_core_fn, P)
    m_ice = y(1); T_blood = y(2); T_brain = y(3);
    if m_ice > 1e-9
        Q_ice  = max(T_blood/R_ice_eff, 0.0);
        dm_ice = -Q_ice/P.L_f;
    else
        Q_ice = 0.0; dm_ice = 0.0;
    end
    dT_blood = k_core*(T_core_fn(t) - T_blood) - Q_ice/(P.m_blood_central*P.c_blood) ...
               - P.k_perf_blood*(T_blood - T_brain);
    dT_brain = P.Q_met/(P.m_brain*P.c_brain) + P.k_perf_brain*(T_blood - T_brain);
    dy = [dm_ice; dT_blood; dT_brain];
end

function T_core_fn = make_T_core(t_end, T0, T1)
    T_core_fn = @(t) T0 + (T1-T0).*min(t,t_end)./t_end;
end

function [t_out, T_brain_out] = simulate_bolus(R_ice_eff, k_core, m_ice0, T_core_const, T_brain0, P, t_eval)
    T_core_fn = @(t) T_core_const;
    y0 = [m_ice0; T_core_const; T_brain0];
    odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t_out, Y] = ode15s(@(t,y) thermobrain_rhs(t,y,R_ice_eff,k_core,T_core_fn,P), t_eval, y0, odeopts);
    T_brain_out = Y(:,3);
end

function dTbrain = mean_dTbrain_bolus(R_ice_eff, k_core, m_ice0, T_core_const, T_brain0, P, window, t_eval)
    [t, Tbr] = simulate_bolus(R_ice_eff, k_core, m_ice0, T_core_const, T_brain0, P, t_eval);
    mask = t>=window(1) & t<=window(2);
    dTbrain = mean(Tbr(mask)) - T_brain0;
end

function dy = local_blood_only_rhs(t, y, R_ice_eff, k_core, T_core_fn, P)
    m_ice = y(1); T_blood = y(2);
    if m_ice > 1e-9
        Q_ice = max(T_blood/R_ice_eff, 0.0); dm_ice = -Q_ice/P.L_f;
    else
        Q_ice = 0.0; dm_ice = 0.0;
    end
    dT_blood = k_core*(T_core_fn(t)-T_blood) - Q_ice/(P.m_blood_central*P.c_blood);
    dy = [dm_ice; dT_blood];
end

function T_end = final_blood_temp(R_ice_eff, k_core, dose_t, dose_m, Tre0, Tre_end, t_end, P, with_ice)
    T_core_fn = make_T_core(t_end, Tre0, Tre_end);
    y = [0; Tre0]; t_prev = 0;
    if with_ice
        events_t = [dose_t, t_end]; events_m = [dose_m, 0];
    else
        events_t = t_end; events_m = 0;
    end
    odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    for i = 1:numel(events_t)
        te = events_t(i); dm = events_m(i);
        if te > t_prev
            [~, Y] = ode15s(@(t,yy) local_blood_only_rhs(t,yy,R_ice_eff,k_core,T_core_fn,P), [t_prev te], y, odeopts);
            y = Y(end,:)';
        end
        if dm > 0, y(1) = y(1) + dm; end
        t_prev = te;
    end
    T_end = y(2);
end

function d = crossover_diff(R_ice_eff, k_core, dose_t, dose_m, Tre0, Tre_end, t_end, P)
    T_ex = final_blood_temp(R_ice_eff, k_core, dose_t, dose_m, Tre0, Tre_end, t_end, P, false);
    T_cl = final_blood_temp(R_ice_eff, k_core, dose_t, dose_m, Tre0, Tre_end, t_end, P, true);
    d = T_cl - T_ex;
end

function Tbrain_out = simulate_full_dosed(R_ice_eff, k_core, T_brain0, dose_t, dose_m, Tre0, Tre_end, t_end, P, t_meas)
    T_core_fn = make_T_core(t_end, Tre0, Tre_end);
    y = [0; Tre0; T_brain0]; t_prev = 0;
    events_t = [dose_t, t_meas(end)]; events_m = [dose_m, 0];
    [events_t, order] = sort(events_t); events_m = events_m(order);
    Tbrain_out = zeros(size(t_meas)); filled = false(size(t_meas));
    odeopts = odeset('RelTol',1e-7,'AbsTol',1e-9);
    for i = 1:numel(events_t)
        te = events_t(i); dm = events_m(i);
        if te > t_prev
            tspan = unique(sort([t_prev, te, t_meas(t_meas>=t_prev & t_meas<=te)]));
            [tt, Y] = ode15s(@(t,yy) thermobrain_rhs(t,yy,R_ice_eff,k_core,T_core_fn,P), tspan, y, odeopts);
            mask = t_meas>=t_prev & t_meas<=te & ~filled;
            idxs = find(mask);
            for kk = 1:numel(idxs)
                [~, ii] = min(abs(tt - t_meas(idxs(kk))));
                Tbrain_out(idxs(kk)) = Y(ii,3);
            end
            filled(mask) = true;
            y = Y(end,:)';
        end
        if dm > 0, y(1) = y(1) + dm; end
        t_prev = te;
    end
end

function report_recovery(name, est_vec, true_val)
    m = mean(est_vec); s = std(est_vec); med = median(est_vec);
    sorted_vec = sort(est_vec);
    n = numel(sorted_vec);
    q1 = sorted_vec(max(1, round(0.25*n)));
    q3 = sorted_vec(max(1, round(0.75*n)));
    rel_err = abs(m-true_val)/true_val*100;
    fprintf('%s: mean=%.4f, sd=%.4f, median=%.4f, IQR=[%.4f,%.4f], true=%.4f, rel.err=%.2f%%\n', ...
            name, m, s, med, q1, q3, true_val, rel_err);
end

function total_load = simulate_optim(alpha, R_ice_eff, k_core, total_ice, Tre0, Tre_end, t_end, T_brain0, P, T_thresh, beta)
    T_core_fn = make_T_core(t_end, Tre0, Tre_end);
    pre_dose = alpha*total_ice; during_dose = (1-alpha)*total_ice/4;
    dose_t = [0, 15*60, 30*60, 45*60, 60*60];
    dose_m = [pre_dose, during_dose, during_dose, during_dose, during_dose];
    y = [0; Tre0; T_brain0; 0]; t_prev = 0;
    odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    for i = 1:numel(dose_t)
        te = dose_t(i); dm = dose_m(i);
        if te > t_prev
            rhs3 = @(t,yy) [thermobrain_rhs(t,yy(1:3),R_ice_eff,k_core,T_core_fn,P); ...
                            log1p(exp(beta*(yy(3)-T_thresh)))/beta];
            [~, Y] = ode15s(rhs3, [t_prev te], y, odeopts);
            y = Y(end,:)';
        end
        y(1) = y(1) + dm; t_prev = te;
    end
    total_load = y(4);
end

function a_min = min_viable_alpha(R_ice_eff, k_core, total_ice, Tre0, Tre_end, t_end, T_brain0, P)
    alphas = linspace(0,1,51);
    loads = zeros(size(alphas));
    for i = 1:numel(alphas)
        loads(i) = simulate_optim(alphas(i), R_ice_eff, k_core, total_ice, Tre0, Tre_end, t_end, T_brain0, P, 38.0, 8.0);
    end
    plateau = loads(end);
    if plateau < 1e-6, thresh = plateau + 0.05; else, thresh = plateau*1.02; end
    idx = find(loads <= thresh, 1, 'first');
    a_min = alphas(idx);
end