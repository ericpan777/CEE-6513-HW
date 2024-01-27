%================================================================================
% 
%   Due 2023/10/25 23:59
%       Georgia Institute of Technology
%       CEE_6513
%       HW_4
%       Student: Yu-Chen Pan
%       GTID: 903918558
%
%================================================================================
%  
%   Problem_1:
%       -d^2(u)/dx^2 = 1
%       u(0) = 0
%       u'(L) = 0
%
%================================================================================
%   
%   Exact solution:
%       u(x) = (-x^2)/2 + L*x
%
%================================================================================
%   
%   Second-order Finite Difference Method:
%
%       u_0  u_1  u_2  ...  ...  ...  u_m  u_m+1
%       |----|----|----|----|----|----|----|
%       0    h    2h   ...  ...  ...       L
%
%       h^(-2)*(u_i-1 - 2*u_i + u_i+1) = f(x_i)
%
%       Au = B
%
%       A = [1    0
%            1   -2    1
%                 1   -2    1
%                      1   -2    1
%                      .    .    .
%                           .    .    .
%                                .    .    .
%                                     1   -2    1
%                                          1/h  -1/h]
%
%       u = [u_0    u_1    u_2    ...    u_m    u_m+1]^T
%
%       B = [0      f(x_1) f(x_2) ...    f(x_m) (h/2)*f(x_m+1)]^T
%
%================================================================================

% Define problem parameters
L = 1;          % Length of the domain
n = 15;         % mesh refinements
f = -1;         % Source term

% Initialize arrays for mesh sizes, matrices, and solutions
m = zeros(n, 1);
h = zeros(n, 1);
As = cell(1, n);
Bs = cell(1, n);
us = cell(1, n);

% Loop over mesh refinements to compute numerical solutions
for i = 1:n
    % Calculate number of internal nodes and mesh size
    m(i) = (2^(i) - 1);
    h(i) = L/(m(i) + 1);
    
    % Initialize matrices and vectors
    A = zeros(m(i) + 2);
    B = zeros(m(i) + 2, 1);

    % Set Dirichlet boundary condition at the left end
    A(1, 1) = 1;
    B(1) = 0;
    
    % Populate A matrix and B vector for internal nodes
    for j = 2:m(i)+1
        A(j, j - 1) = 1/(h(i)^2);
        A(j, j) = -2/(h(i)^2);
        A(j, j + 1) = 1/(h(i)^2);
        B(j) = f;
    end
    
    % Set Neumann boundary condition at the right end
    A(m(i) + 2, m(i) + 1) = 1/h(i);
    A(m(i) + 2, m(i) + 2) = -1/h(i);
    B(m(i) + 2) = (h(i)/2)*f;
    
    % Solve the system using sparse matrix
    A_sparse = sparse(A);
    u = A_sparse\B;

    % Store matrices and solution
    As{i} = A;
    Bs{i} = B;
    us{i} = u;
end

% Plot numerical solutions for different mesh sizes
figure;
for i = 1:n
    x = (0:h(i):L)';
    plot(x, us{i}, '-', 'DisplayName', sprintf('mesh size = %.3e', h(i)));
    hold on;
end

% Plot the exact solution
x_exact = linspace(0, L, 5000);
y_exact = (-x_exact.^2)/2 + x_exact;
plot(x_exact, y_exact, 'k:', 'LineWidth', 1, 'DisplayName', 'Exact Solution');

% Set plot properties for convergence plot
xlabel('x', 'FontSize', 18); 
ylabel('u', 'FontSize', 18); 
title('Convergence', 'FontSize', 24); 
lgd = legend('show', 'Location','southeast');
fontsize(lgd,16,'points')
grid on;

% Compute L2-norm errors for each mesh size using a fine grid
num_fine_points = 5000;
x_fine = linspace(0, L, num_fine_points);
L2_n_errors = zeros(1, n);

for i = 1:n
    u_interpolated = interp1((0:h(i):L)', us{i}, x_fine, 'linear');
    y_exact_fine = (-x_fine.^2)/2 + x_fine;
    L2_n_errors(i) = norm(u_interpolated - y_exact_fine, 2) / norm(y_exact_fine, 2);
end

% Display L2-norm errors
disp('Errors for each mesh size:');
disp(L2_n_errors);

% Compute RMS errors for each mesh size using a fine grid
rms_errors = zeros(1, n);

for i = 1:n
    u_interpolated = interp1((0:h(i):L)', us{i}, x_fine, 'linear');
    y_exact_fine = (-x_fine.^2)/2 + x_fine;
    
    % RMS Error computation
    differences = u_interpolated - y_exact_fine;
    rms_errors(i) = sqrt(mean(differences.^2));
end

% Display RMS errors
disp('RMS Errors for each mesh size:');
disp(rms_errors);

% Compute and display convergence rates using RMS error
convergence_rates = zeros(1, n-1);
for i = 1:n-1
    convergence_rates(i) = log(rms_errors(i)/rms_errors(i+1)) / log(2);
end
disp('Convergence rates:');
disp(convergence_rates);

% Plot L2-norm error vs. mesh size (log-log)
figure;
loglog(h, L2_n_errors, '-o', 'LineWidth', 2, 'MarkerSize', 3);
xlabel('Mesh Size (h)', 'FontSize', 18);
ylabel('Error', 'FontSize', 18);
title('L2-norm Error vs. Mesh Size (log-log)', 'FontSize', 24);
grid on;

% Plot L2-norm error vs. mesh size (linear)
figure;
plot(h, L2_n_errors, '-o', 'LineWidth', 2, 'MarkerSize', 3);
xlabel('Mesh Size (h)', 'FontSize', 18);
ylabel('Error', 'FontSize', 18);
title('L2-norm Error vs. Mesh Size (linear)', 'FontSize', 24);
grid on;

% Plot RMS error vs. mesh size (log - log)
figure;
loglog(h, rms_errors, '-o', 'LineWidth', 2, 'MarkerSize', 3);
xlabel('Mesh Size (h)', 'FontSize', 18);
ylabel('Error', 'FontSize', 18);
title('RMS Error vs. Mesh Size (log-log)', 'FontSize', 24);
grid on;

% Plot RMS error vs. mesh size (linear)
figure;
plot(h, rms_errors, '-o', 'LineWidth', 2, 'MarkerSize', 3);
xlabel('Mesh Size (h)', 'FontSize', 18);
ylabel('Error', 'FontSize', 18);
title('RMS Error vs. Mesh Size (linear)', 'FontSize', 24);
grid on;

% Plot convergence rate vs. mesh size
figure;
semilogx(h(1:end-1), convergence_rates, '-o', 'LineWidth', 2, 'MarkerSize', 3);
xlabel('Mesh Size (h)', 'FontSize', 18);
ylabel('Convergence Rate', 'FontSize', 18);
title('Convergence Rate vs. Mesh Size', 'FontSize', 24);
ylim([1.9 2.1]);
grid on;

