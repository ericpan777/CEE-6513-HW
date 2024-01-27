% Clear the environment and close figures
clc;
close all;

% Set the parameters for spatial discretization
L = 1;
dx = 0.01;
Nx = L / dx + 1;
x = linspace(0, L, Nx);

% Set the parameters for temporal discretization
T = 100;
dt = 0.01;
Nt = T / dt + 1;
t = linspace(0, T, Nt);

% Wave speed and Courant number
v = 1;
c = v * (dt / dx);

% Initialize solution matrix
U = zeros(Nx, Nt);

% Define the time index for the initial condition
s1 = floor(pi / dt);

% Set the initial condition
U(1, 1:s1) = 1 - cos(2 * t(1:s1));

% Perform the finite difference scheme
for n = 3:Nt-1
    for i = 2:Nx-1
        U(i, n) = applyScheme(U, i, n, c^2);
    end
end

% Plot the results at specific time step
plot_time_step = [20, 40, 60, 80, 100];
figNum = 1;

for i = 1:length(plot_time_step)
    plotWave(U, x, t, plot_time_step(i), figNum);
    figNum = figNum + 1;
end

plotAllWaves(U, x, t, plot_time_step, figNum);
figNum = figNum + 1;

% Plot the results at specific time step
plot_time_step = [20, 40, 60, 80, 100];
plot_times = plot_time_step/dt;

for i = 1:length(plot_times)
    plotWave(U, x, t, plot_times(i), figNum);
    figNum = figNum + 1;
end

plotAllWaves(U, x, t, plot_times, figNum);

% Create a movie for the traveling wave
createWaveMovie(U, x, t, Nt);

function Uj = applyScheme(U, i, n, c2)
    % Applies the finite difference scheme for a single grid point
    U1 = 2 * U(i, n-1) - U(i, n-2);
    U2 = U(i-1, n-1) - 2 * U(i, n-1) + U(i+1, n-1);
    Uj = U1 + c2 * U2;
end

function plotWave(U, x, t, k, figNum)
    % Plot the wave at a given time step
    figure(figNum);
    plot(x, U(:, k), 'linewidth', 2);
    grid on;
    axis([min(x), max(x), -2, 2]);
    xlabel('X axis', 'FontSize', 14);
    ylabel('Wave Amplitude', 'FontSize', 14);
    title(sprintf('TIME STEP = %d, TIME = %.2f seconds', k, t(k)), 'FontSize', 14);
    set(gca, 'FontSize', 14);
    set(gcf, 'color', 'white');
end

function plotAllWaves(U, x, t, plot_times, figNum)
    % Plot all the selected time steps in one figure
    figure(figNum);
    hold on;
    colors = lines(numel(plot_times));
    for i = 1:length(plot_times)
        k = plot_times(i);
        plot(x, U(:, k), 'linewidth', 2, 'Color', colors(i, :));
        legendInfo{i} = sprintf('TIME STEP = %d, TIME = %.2f seconds', k, t(k));
    end
    hold off;
    
    grid on;
    axis([min(x), max(x), -2, 2]);
    xlabel('X axis', 'FontSize', 14);
    ylabel('Wave Amplitude', 'FontSize', 14);
    title('Wave States at Different Times', 'FontSize', 14);
    legend(legendInfo, 'Location', 'northeastoutside');
    set(gca, 'FontSize', 14);
    set(gcf, 'color', 'white');
end


function createWaveMovie(U, x, t, Nt)
    % Create a movie for the traveling wave
    Speed = 2;
    frames = struct('cdata', [], 'colormap', []);
    movieFigure = figure;
    set(movieFigure, 'color', 'white');

    pause(20);

    for j = 1:floor((Nt-1) / Speed)
        plotWave(U, x, t, j * Speed, movieFigure);
        frames(j) = getframe(movieFigure);
    end

    % Play back the movie
    movieFigurePlayback = figure; % Create a new figure for playback
    movie(movieFigurePlayback, frames, 1, 10);
end
