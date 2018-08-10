% Author: Hua Feng

clear all

syms R M M_dot alpha f eta % independent variables
syms H rho Sigma p c_s nu T_c % dependent variables
syms G c sigma k_es % constants

eq1 = rho == Sigma / 2 / H;
eq2 = H == c_s * R^(3/2) / (G * M)^(1/2);
eq3 = c_s^2 == p / rho;
eq4 = p == 4 * sigma / (3 * c) * T_c^4;
eq5 = 3 * G * M * M_dot / (8 * pi * R^3) * f == M_dot / (2 * pi * R^2) * p / rho;
eq6 = nu * Sigma == M_dot / 3 / pi * f;
eq7 = nu == 2 / 3 * alpha * c_s * H;
i = 8; % which solution


S = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], [H, rho, Sigma, p, c_s, nu, T_c]);

% disp(simplify(S.H(i), 'IgnoreAnalyticConstraints', true, 'Steps', 100))
% disp(simplify(S.T_c(i), 'IgnoreAnalyticConstraints', true, 'Steps', 100))

% find the trapping radius where Q_rad = Q_adv
eq_t = 4 * sigma * S.T_c(i)^4 / (k_es * S.H(i)) == M_dot / (2 * pi * R^2) * S.p(i);
R_t = solve(eq_t, R);

% constants in cgs
M_sun = 1.98855e33;
G = 6.6726e-8;
k_es = 0.342;
eta = 0.1;
c = 2.99792458e10;
sigma = 5.6705e-5;
R_isco_sun = 6 * G * M_sun / c^2;
L_Edd_sun = 4 * pi * G * M_sun * c / k_es;
M_Edd_sun = L_Edd_sun / (eta * c^2);

disp('R_tr')
disp(simplify(R_t(end), 'IgnoreAnalyticConstraints', true, 'Steps', 100))
disp(eval(subs(R_t(end), [M_dot, f], [M_Edd_sun, 1.0])))

disp('rho')
disp(eval(subs(S.rho(i), [M, M_dot, R, alpha, f], [M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0])))

disp('p')
disp(eval(subs(S.p(i), [M, M_dot, R, alpha, f], [M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0])))

disp('H')
disp(eval(subs(S.H(i), [M, M_dot, R, alpha, f], [M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0])))

disp('T_c')
disp(eval(subs(S.T_c(i), [M, M_dot, R, alpha, f], [M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0])))


disp('L_bb')
R_bb = eval(subs(R_t(end), [M_dot, f], [M_Edd_sun, 1.0]));
T_bb = eval(subs(S.T_c(i), [R, M, M_dot, R, alpha, f], [R_bb, M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0]));
rho_bb = eval(subs(S.rho(i), [R, M, M_dot, R, alpha, f], [R_bb, M_sun, M_Edd_sun, R_isco_sun, 1.0, 1.0]));
tau_es = k_es * rho_bb * R_bb;
L_bb = 4 * pi * sigma * R_bb^2 * T_bb^4 / tau_es;
disp(L_bb / L_Edd_sun)
