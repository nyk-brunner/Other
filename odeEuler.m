function [x, y] = odeEuler(ode, a, b, Y, h, opt)
% odeEuler uses Euler's explicit method (also called "forward Euler method")
% for solving first order ordinary differential equations given an initial
% value (IVP-ODE)
%
% Inputs:
%   First order, ordinary differential equation that calculates dy/dx, ODE
%   Left limit of 'x', a
%   Right limit of 'x', b
%   Value of solution at initial point, Y
%   Plotting option, opt
%       If you want to plot solutions set 'opt' = plot
%
% Outputs:
%   Vector with x-coordinates of the solution, xs
%   Vector with y-coordinates of the solution, ys
%
% Syntax:
%   [x, y] = odeEuler(ode, a, b, Y, h, opt)

% Check input variables:
if nargin == 5 % No options given
    opt = 'no plot';
end

% Initialize:
N = (b - a)/h; % Calculate number of points
x = zeros(N+1,1); % Initialize x-vector
y = zeros(N+1,1); % Initialize y-vector

% Apply initial value:
x(1) = a;
y(1) = Y;

for i = 1:N
    x(i+1) = x(i) + h;
    y(i+1) = y(i) + ode(x(i),y(i))*h;
end

if strcmpi(opt,'plot')
    figure(1)
        plot(x,y,'k-','linewidth',1.5)
        
        grid on
        axis tight
        xlabel('x')
        ylabel('y')
        title({'Solved ODE',['N = ',num2str(N),', h = ',num2str(h)]})
end
