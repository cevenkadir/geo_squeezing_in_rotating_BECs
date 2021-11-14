clear;
clc;

global M_x; % division number of spatial mesh along the x-axis
global M_y; % division number of spatial mesh along the y-axis
global dx; % unitless spatial mesh size for the x-axis
global dy; % unitless spatial mesh size for the y-axis
global x; % unitless spatial array for the x-axis
global y; % unitless spatial array for the y-axis
global xx; % unitless spatial mesh for the x-axis
global yy; % unitless spatial mesh for the y-axis
global xi; % unitless Fourier frequency array with respect to the x-axis
global nu; % unitless Fourier frequency array with respect to the y-axis
global xxi; % unitless Fourier frequency mesh with respect to the x-axis
global nnu; % unitless Fourier frequency mesh with respect to the y-axis
global V; % unitless potential
global E; % unitless energy of state
global phi_n; % unitless wave function
global L; L = 20; % unitless half length of bounded domain
global N; N = (8.1e5)^(4/5); % number of atoms
global g; g = 6.552e-2; % unitless 2D unitless coupling constant
global gamma_y; gamma_y = 1; % omega_y/omega_x
global Omega; Omega = 0.8; % unitless rotation frequency
global initial_wf; initial_wf = 'phi_a'; % initial wave function
global epsilon_energy; epsilon_energy = 1e-12; % stopping criteria value
% for the energy minimization

% for the multi-grid method
p_initial = 6;
p_end = 9;
Ms = 2.^(p_initial:p_end);

CPU_t_initial = cputime; % storing the time value at the beginning

% starting the multigrid method
for M_i = Ms
    M_x = M_i;
    M_y = M_i;

    dx = 2 * L / (M_x);
    dy = 2 * L / (M_y);

    x = (-M_x / 2:M_x / 2 - 1) * dx;
    y = (-M_y / 2:M_y / 2 - 1) * dy;

    % if the current division number is not equal to 2^p_initial,
    if M_i ~= Ms(1)
        % store the last spatial meshes for the interpolation
        xx_last = xx;
        yy_last = yy;

        % calculate new spatial meshes
        [xx, yy] = meshgrid(x, y);

        % interpolate
        phi_n = interp2(xx_last, yy_last, phi_n, xx, yy, 'cubic');
        phi_n(isnan(phi_n)) = 0;
    else
        % calculate new spatial meshes
        [xx, yy] = meshgrid(x, y);

        % defining the functions of possible initial wave functions
        phi_0s.phi_a = calc_phi_a;
        phi_0s.phi_b = calc_phi_b;
        phi_0s.phi_b_bar = calc_phi_b_bar;
        phi_0s.phi_c = calc_phi_c;
        phi_0s.phi_c_bar = calc_phi_c_bar;
        phi_0s.phi_d = calc_phi_d;
        phi_0s.phi_d_bar = calc_phi_d_bar;
        phi_0s.phi_e = calc_phi_e;
        phi_0s.phi_e_bar = calc_phi_e_bar;

        % initialize the wave function with the chosen state
        phi_n = phi_0s.(initial_wf)();
    end

    xi = (-M_x / 2:M_x / 2 - 1) .* 2 * pi / (2 * L);
    nu = (-M_y / 2:M_y / 2 - 1) .* 2 * pi / (2 * L);
    [xxi, nnu] = meshgrid(xi, nu);

    V = 0.5 * (xx.^2 + gamma_y^2 * yy.^2);

    phi_n = normalize(phi_n);

    E = calc_E(phi_n);

    % running the PG_C method once
    [mu_n, r_n, p_n, phi_n] = PG_C();

    E = calc_E(phi_n);

    count = 0;

    % run the PCG_C method until energy_err_n is less than epsilon_energy
    while true
        r_n_minus_1 = r_n;
        p_n_minus_1 = p_n;

        [mu_n, r_n, p_n, phi_n] = PCG_C(r_n_minus_1, p_n_minus_1);

        E_new = calc_E(phi_n);

        count = count + 1;

        energy_err_n = abs(E-E_new);

        E = E_new;

        if energy_err_n < epsilon_energy
            break
        end

    end

end

CPU_t_final = cputime; % storing the time value at the end
CPU_t = CPU_t_final - CPU_t_initial; % calculating the CPU time

prob_density = abs(phi_n).^2; % calculating the probability density

hold on

% setting the interpreter of plot to LaTeX
set(gca, 'DefaultTextInterpreter', 'latex');
set(gca, 'defaultAxesTickLabelInterpreter', 'latex');
set(gca, 'defaulttextinterpreter', 'latex');
set(gca, 'defaultLegendInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

% plotting the density
surf(xx, yy, prob_density, 'EdgeColor', 'None', 'facecolor', 'interp')

view(2) % displaying the plot in a 2D view
daspect([1, 1, 1]) % using equal unit lengths along each axis
colormap jet % setting current colormap to jet

% specifying the axis limits
xlim([-L, L - dx])
ylim([-L, L - dy])

% labeling each axis
xlabel('$x$')
ylabel('$y$')

% setting the plot title
title(['$N = $ ', num2str(N),'$, g = $ ', num2str(g), ', $\Omega = $ ', num2str(Omega), ', $E =$ ', num2str(real(E), '%.5f'), ', CPU time: ', num2str(CPU_t), ' s'])

% creating a colorbar
c = colorbar;
c.Label.String = '$|\psi(x,y)|^2$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

hold off

% preparing the filname
Omega_str_wo_dot = num2str(Omega, '%.2f');
Omega_str_wo_dot(Omega_str_wo_dot == '.') = [];
g_str_wo_dot = num2str(g, '%.2f');
g_str_wo_dot(g_str_wo_dot == '.') = [];
filename = ['N_', num2str(N, '%.0f'), '_g_', g_str_wo_dot, '_Omega_', Omega_str_wo_dot];

% saving several variables in .MAT format
save([filename, '.mat'], 'phi_n', 'g', 'Omega', 'gamma_y', 'N', 'L', 'M_x', 'M_y', 'x', 'y', 'mu_n', 'E')

% saving the plot in .EPS format
exportgraphics(gcf, [filename, '.eps'], 'Resolution', 300)

%-------------------------------------------------------------------------%
% FUNCTIONS

function output = allsum(x)
output = sum(x, 'all');
end

function phi_a = calc_phi_a()
% calculating (a) initial condition

global xx;
global yy;

phi_a = exp(-(xx.^2 + yy.^2)/2) / sqrt(pi);
end

function phi_b = calc_phi_b()
% calculating (b) initial condition

global xx;
global yy;

phi_a = calc_phi_a();
phi_b = phi_a .* (xx + 1j * yy);
end

function phi_b_bar = calc_phi_b_bar()
% calculating (b-bar) initial condition

phi_b = calc_phi_b();
phi_b_bar = conj(phi_b);
end

function phi_c = calc_phi_c()
% calculating (c) initial condition

phi_a = calc_phi_a();
phi_b = calc_phi_b();
phi_c = (phi_a + phi_b) ./ 2;
end

function phi_c_bar = calc_phi_c_bar()
% calculating (c-bar) initial condition

phi_c = calc_phi_c();
phi_c_bar = conj(phi_c);
end

function phi_d = calc_phi_d()
% calculating (d) initial condition

global Omega;

phi_a = calc_phi_a();
phi_b = calc_phi_b();
phi_d = (1 - Omega) * phi_a + Omega * phi_b;
end

function phi_d_bar = calc_phi_d_bar()
% calculating (d-bar) initial condition

phi_d = calc_phi_d();
phi_d_bar = conj(phi_d);
end

function phi_e = calc_phi_e()
% calculating (e) initial condition
global Omega;

phi_a = calc_phi_a();
phi_b = calc_phi_b();
phi_e = (1 - Omega) * phi_b + Omega * phi_a;
end

function phi_e_bar = calc_phi_e_bar()
% calculating (d-bar) initial condition

phi_e = calc_phi_e();
phi_e_bar = conj(phi_e);
end

function normalized_state = normalize(state)
% normalizing the given state

global dx;
global dy;
global N;

integral = allsum(abs(state).^2) * dx * dy;
norm_factor = sqrt(N/integral);

normalized_state = norm_factor * state;
end

function H_times_state = H_times_state(state)
% calculating H times the given state

global xx;
global yy;
global xxi
global nnu;
global Omega;
global g;
global V;
global phi_n;

H_1 = 0.5 * xxi.^2 + Omega * yy .* xxi;
state_kp = fftshift(fft(state, [], 2), 2); % FFT along the x-axis
state_kp = H_1 .* state_kp;
H_1_times_state = ifft(ifftshift(state_kp, 2), [], 2); % IFFT along the
% x-axis

H_2 = 0.5 * nnu.^2 - Omega * xx .* nnu;
state_qj = fftshift(fft(state, [], 1), 1); % FFT along the y-axis
state_qj = H_2 .* state_qj;
H_2_times_state = ifft(ifftshift(state_qj, 1), [], 1); % IFFT along the
% y-axis

H_times_state = H_1_times_state + H_2_times_state + (V + g * abs(phi_n).^2) .* state;
end

function E = calc_E(state)
% calculating the energy of given state

global xx;
global yy;
global xxi
global nnu;
global Omega;
global g;
global V;
global dx;
global dy;
global N;

H_1_op = 0.5 * xxi.^2 + Omega * yy .* xxi;
state_kp = fftshift(fft(state, [], 2), 2); % FFT along the x-axis
state_kp = H_1_op .* state_kp;
H_1_times_state = ifft(ifftshift(state_kp, 2), [], 2); % IFFT along the
% x-axis

H_2_op = 0.5 * nnu.^2 - Omega * xx .* nnu;
state_qj = fftshift(fft(state, [], 1), 1); % FFT along the y-axis
state_qj = H_2_op .* state_qj;
H_2_times_state = ifft(ifftshift(state_qj, 1), [], 1); % IFFT along the
% y-axis

H_state = H_1_times_state + H_2_times_state + (V + 0.5 * g * abs(state).^2) .* state;

E = allsum((conj(state).*H_state)) * dx * dy / N;
E = real(E);
end

function mu = calc_mu()
% calculating the chemical energy of phi_n

global phi_n;
global N;
global dx;
global dy;

mu = allsum((conj(phi_n).*H_times_state(phi_n))) * dx * dy / N;
mu = real(mu);
end

function P_C_times_state = calc_P_C_times_state(state)
% calculating the preconditioner P_C times the given state

global xxi
global nnu;
global g;
global V;
global phi_n;
global dx;
global dy;
global N;

phi_n_qp = fftshift(fft2(phi_n));
phi_n_qp = -(xxi.^2 + nnu.^2) .* phi_n_qp;
lap_times_phi_n = ifft2(ifftshift(phi_n_qp));

alpha_Delta = -0.5 * conj(phi_n) .* lap_times_phi_n + V .* abs(phi_n).^2 + g * abs(phi_n).^4;
alpha_Delta = real(sum(sum(alpha_Delta))*dx*dy/N);

alpha_V = alpha_Delta;

P_V = sqrt(alpha_V+V+g*abs(phi_n).^2).^(-1);
P_Delta = (alpha_Delta + 0.5 * (xxi.^2 + nnu.^2)).^(-1);

P_C_times_state = P_V .* state;

P_C_times_state_qp = fftshift(fft2(P_C_times_state));
P_C_times_state_qp = P_Delta .* P_C_times_state_qp;
P_C_times_state = ifft2(ifftshift(P_C_times_state_qp));

P_C_times_state = P_V .* P_C_times_state;
end

function p_n_norm = calc_p_n_norm(p_n)
% calculating the norm of p_n

global dx;
global dy;

p_n_norm = sqrt(allsum(abs(p_n).^2)*dx*dy);
end

function [mu_n, r_n, p_n, phi_n_plus_1] = PG_C()
% finding phi_n_plus_1 by running the PG_C method once

global phi_n;
global dx;
global dy;
global N;
global g;
global E;

mu_n = calc_mu();

r_n = H_times_state(phi_n) - mu_n * phi_n;

d_n = -calc_P_C_times_state(r_n);

p_n = d_n - real(allsum((d_n) .* conj(phi_n))*dx*dy/N) .* phi_n;

p_n_norm = calc_p_n_norm(p_n);

theta_n_numerator = -real(allsum((2*H_times_state(phi_n)) .* conj(p_n))*dx*dy/N) * p_n_norm;
rho_p_phi = real(phi_n.*conj(p_n));
g_n = 2 * g * rho_p_phi .* phi_n;
theta_n_denominator = 2 * sqrt(N) * (allsum(conj(H_times_state(p_n)) .* (p_n) + real(conj(p_n) .* g_n)) * dx * dy / N - mu_n * p_n_norm^2 / N);

% if the denominator of theta_n is positive
if theta_n_denominator > 0
    % choose theta_n as theta_opt_n
    theta_n = real(theta_n_numerator/theta_n_denominator);
else
    % choose theta_n as a small positive value
    theta_n = real(theta_n_numerator);
end

% decrease theta_n until the energy starts to decrease or theta_n is less
% than eps (2.2204e-16)
while true
    phi_n_plus_1_temp = cos(theta_n) .* phi_n + sin(theta_n) .* p_n * sqrt(N) / p_n_norm;

    E_temp = calc_E(phi_n_plus_1_temp);

    if E_temp < E || theta_n < eps
        break
    else
        theta_n = theta_n / 2;
    end
end

phi_n_plus_1 = cos(theta_n) .* phi_n + sin(theta_n) .* p_n * sqrt(N) / p_n_norm;
phi_n_plus_1 = normalize(phi_n_plus_1);
end

function [mu_n, r_n, p_n, phi_n_plus_1] = PCG_C(r_n_minus_1, p_n_minus_1)
% finding phi_n_plus_1 by running the PCG_C method once

global phi_n;
global dx;
global dy;
global N;
global g;
global E;

mu_n = calc_mu();

r_n = H_times_state(phi_n) - mu_n * phi_n;

% calculate beta_n_PR
beta_n_PR = real((allsum((r_n - r_n_minus_1) .* conj(calc_P_C_times_state(r_n)))*dx*dy/N)/(allsum((r_n_minus_1) .* conj(calc_P_C_times_state(r_n_minus_1))) * dx * dy / N));
% if beta_n_PR is negative, take beta_n as 0
beta_n = max([beta_n_PR, 0]);

d_n = -calc_P_C_times_state(r_n) + beta_n * p_n_minus_1;

p_n = d_n - real(allsum((d_n) .* conj(phi_n))*dx*dy/N) .* phi_n;

p_n_norm = calc_p_n_norm(p_n);

theta_n_numerator = -real(allsum((2*H_times_state(phi_n)) .* conj(p_n))*dx*dy/N) * p_n_norm;
rho_p_phi = real(phi_n.*conj(p_n));
g_n = 2 * g * rho_p_phi .* phi_n;
theta_n_denominator = 2 * sqrt(N) * (allsum(conj(H_times_state(p_n)) .* (p_n) + real(conj(p_n) .* g_n)) * dx * dy / N - mu_n * p_n_norm^2 / N);

% if the numerator of theta_n is negative
if theta_n_numerator < 0
    % run the PG_C method once
    [mu_n, r_n, p_n, phi_n_plus_1] = PG_C();
else
    % choose theta_n as theta_opt_n
    
    theta_n = real(theta_n_numerator/theta_n_denominator);
    
    % decrease theta_n until the energy starts to decrease or theta_n is
    % less than eps (2.2204e-16)
    while true
        phi_n_plus_1_temp = cos(theta_n) .* phi_n + sin(theta_n) .* p_n * sqrt(N) / p_n_norm;

        E_temp = calc_E(phi_n_plus_1_temp);

        if E_temp < E || theta_n < eps
            break
        else
            theta_n = theta_n / 2;
        end
    end

    phi_n_plus_1 = cos(theta_n) .* phi_n + sin(theta_n) .* p_n * sqrt(N) / p_n_norm;

    phi_n_plus_1 = normalize(phi_n_plus_1);
end
end