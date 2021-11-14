clear;
clc;

% uncomment the following line to load the MAT-File (exported by
% multigrid_PCG_C_method.m).
% load('N_[]_g_[]_Omega_[].mat')

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
global phi_n; % unitless wave function
global L; % unitless half length of bounded domain
global N; % number of atoms
global g; % unitless 2D unitless coupling constant
global gamma_y; % omega_y/omega_x
global Omega; Omega = 1; % unitless rotation frequency
global epsilon_el; epsilon_el = 0.125; % trap ellipticity
global dt; dt = 1e-3; % unitless time step size
global t_last; t_last = 6.4;

zeta = epsilon_el / 2;

dx = 2 * L / (M_x);
dy = 2 * L / (M_y);

x = (-M_x / 2:M_x / 2 - 1) * dx;
y = (-M_y / 2:M_y / 2 - 1) * dy;
[xx, yy] = meshgrid(x, y);

xi = (-M_x / 2:M_x / 2 - 1) .* 2 * pi / (2 * L);
nu = (-M_y / 2:M_y / 2 - 1) .* 2 * pi / (2 * L);
[xxi, nnu] = meshgrid(xi, nu);

% harmonic oscillator + saddle
V = 0.5 * (xx.^2 + gamma_y^2 * yy.^2) + 0.5 * epsilon_el * (xx.^2 - yy.^2);

% opening a video file
Omega_str_wo_dot = num2str(Omega, '%.2f');
Omega_str_wo_dot(Omega_str_wo_dot == '.') = [];
g_str_wo_dot = num2str(g, '%.2f');
g_str_wo_dot(g_str_wo_dot == '.') = [];
filename = ['t_evol_with_V_S_', 'N_', num2str(N, '%.0f'), '_g_', g_str_wo_dot, '_Omega_', Omega_str_wo_dot];
v = VideoWriter(filename);
open(v);

majors = []; % initializing an array for major radius
minors = []; % initializing an array for minor radius
N_t_i = 0;
while (N_t_i <= t_last / dt)
    H_1x = 0.5 * xxi.^2 + Omega * yy .* xxi;
    state_kp = fftshift(fft(phi_n, [], 2), 2); % FFT along the x-axis
    state_kp = exp(-0.5j*(N_t_i > 0)*dt*H_1x) .* state_kp;
    phi_n = ifft(ifftshift(state_kp, 2), [], 2); % IFFT along the x-axis
    phi_n = normalize(phi_n);

    H_1y = 0.5 * nnu.^2 - Omega * xx .* nnu;
    state_qj = fftshift(fft(phi_n, [], 1), 1); % FFT along the y-axis
    state_qj = exp(-0.5j*(N_t_i > 0)*dt*H_1y) .* state_qj;
    phi_n = ifft(ifftshift(state_qj, 1), [], 1); % IFFT along the y-axis
    phi_n = normalize(phi_n);

    H_2 = V + g * abs(phi_n).^2;
    phi_n = exp(-1j*(N_t_i > 0)*dt*H_2) .* phi_n;
    phi_n = normalize(phi_n);

    H_1y = 0.5 * nnu.^2 - Omega * xx .* nnu;
    state_qj = fftshift(fft(phi_n, [], 1), 1); % FFT along the y-axis
    state_qj = exp(-0.5j*(N_t_i > 0)*dt*H_1y) .* state_qj;
    phi_n = ifft(ifftshift(state_qj, 1), [], 1); % IFFT along the y-axis
    phi_n = normalize(phi_n);

    H_1x = 0.5 * xxi.^2 + Omega * yy .* xxi;
    state_kp = fftshift(fft(phi_n, [], 2), 2); % FFT along the x-axis
    state_kp = exp(-0.5j*(N_t_i > 0)*dt*H_1x) .* state_kp;
    phi_n = ifft(ifftshift(state_kp, 2), [], 2); % IFFT along the x-axis
    phi_n = normalize(phi_n);

    % finding the major and minor radii
    [major_n, minor_n] = find_major_minor_radii();
    majors = [majors, major_n];
    minors = [minors, minor_n];

    % plotting the density distribution
    surf(xx, yy, abs(phi_n).^2, 'EdgeColor', 'None', 'facecolor', 'interp');
    title(['zeta t = ', num2str(zeta * N_t_i * dt, '%.3f')]);

    view(2); % displaying the plot in a 2D view
    daspect([1, 1, 1]); % using equal unit lengths along each axis
    colormap jet; % setting current colormap to jet

    % specifying the axis limits
    xlim([-L, L - dx]);
    ylim([-L, L - dy]);

    % labeling each axis
    xlabel('x')
    ylabel('y')

    % get a frame from the figure and write it into the video file
    frame = getframe(gcf);
    drawnow;
    writeVideo(v, frame);

    N_t_i = N_t_i + 1;
end

% closing the video file
close(v);

% plotting time evolution of radii
figure(2)
hold on;

% setting the interpreter of plot to LaTeX
set(gca, 'DefaultTextInterpreter', 'latex');
set(gca, 'defaultAxesTickLabelInterpreter', 'latex');
set(gca, 'defaulttextinterpreter', 'latex');
set(gca, 'defaultLegendInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

t = 0:dt:t_last; % time array

plot(zeta*t, majors, 'DisplayName', '$R_+$')
plot(zeta*t, minors, 'DisplayName', '$R_-$')

A = (majors(1) + minors(1)) / 2; % free parameter of exp. fit

plot(zeta*t, A*exp(zeta * t), '--', 'DisplayName', '$A \exp (\zeta t)$')
plot(zeta*t, A*exp(-zeta * t), '--', 'DisplayName', '$A \exp (\zeta t)$')

legend('Location', 'northwest', 'Interpreter', 'latex')

xlabel('$\zeta t$')
ylabel('$R$')

hold off;

%-------------------------------------------------------------------------%
% FUNCTIONS

function output = allsum(x)
output = sum(x, 'all');
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

function [major_n, minor_n] = find_major_minor_radii()
% calculating the major and minor radii of cloud scaled by ell_B

global phi_n;
global dx;

% converting the density matrix to a grayscale image
img = mat2gray(abs(phi_n).^2);

% binarizing the grayscale image
binarized_img = imbinarize(img, 'adaptive');

% extracting the extra area from the image
binarized_img = imclearborder(binarized_img);
binarized_img = bwareafilt(binarized_img, 1);

% calculating the major and minor axis length of the ellipse
s = regionprops(binarized_img, {'MajorAxisLength', 'MinorAxisLength'});

major_n = s.MajorAxisLength * dx * sqrt(2) / 2;
minor_n = s.MinorAxisLength * dx * sqrt(2) / 2;
end