function [x, y] = odeEulerMod(ode, a, b, h, Y, opt)
% odeEulerMod uses the "modified" Euler's algorithm to solve a first order
% ODE 
%
% Inputs:
%   Ordinary differential equation, ODE
%   Left bound of "x", a
%   Right bound of "x", b
%   Step size, h
%   Value of the solution at y (initial value), Y
%   Plotting option, opt
%       If you want to the function to plot, set opt = plot
%
% Outputs:
%   Vector containing the "independent" variable solution, x
%   Vector containing the "dependent" variable solution, y
%
% Syntax:
%   [x, y] = odeEulerMod(ode, a, b, h, Y)

if nargin == 5
    opt = 'no';
end

% Initialize:
N = (b - a)/h; % Number of points (N-1 points)
x = zeros(N+1,1); % Initialize x-vector
y = zeros(N+1,1); % Initialize y-vector

% Apply initial conditions:
x(1) = a;
y(1) = Y;

for i = 1:N
    x(i+1) = x(i) + h;
    m1 = ode(x(i),y(i)); % Slope at beginning of interval
    yEu = y(i) + m1*h; 
    m2 = ode(x(i+1), yEu); % Slope at end of interval
    y(i+1) = y(i) + (m1 + m2)*h/2; % Calculate new point
end

% Plot:
if strcmpi(opt, 'plot')
    figure(1)
        plot(x,y,'k-','linewidth',1.5)
        
        grid on
        xlabel('X')
        ylabel('Y')
        title('Fuck this guy')
end
