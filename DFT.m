function [kp, Ak, Bk, fk] = DFT(t, f, opt)
% DFT determines the real Discrete Fourier Transform of a function given by
% a finite number of points.
%
% Inputs:
%   Temporal (time) vector (independent variable), t
%   Values of the dependent variable, f
%
% Outputs:
%   Vector containing "k" values, k
%   Vector containing the "Ak" constants
%   Vector containing the "Bk" constants
%
% Syntax:
%   [kp, Ak, Bk, fk] = DFT(t, f, opt)

% Check inputs:
if strcmpi(class(f),'function_handle')
    % Actual function:
    t_plot = linspace(0,t(end)+(t(2)-t(1)),length(t)*3); % Time vector
    y_plot = f(t_plot); % Actual function
    f = f(t); % Replace function handle with sampled data
    T = true;
else
    T = false;
end

if nargin == 2
    opt = 'no';
end

N = length(f)/2; % Number of points
dt = (t(2*N) - t(1))/(2*N - 1); % Time step
kp = 0:N; 
Ak = zeros(N+1,1);
Bk = zeros(N+1,1);

% Apply "initial conditions":
Bk(1) = sum(f(1:2*N))/(2*N);

for k = 2:N
    for j = 1:2*N
        Ak(k) = Ak(k) + f(j)*sin(pi*(k-1)*t(j)/(dt*N));
        Bk(k) = Bk(k) + f(j)*cos(pi*(k-1)*t(j)/(dt*N));
    end
    Ak(k) = Ak(k)/N;
    Bk(k) = Bk(k)/N;
end

for j = 1:2*N
    Bk(N+1) = Bk(N+1) + f(j)*cos(pi*N*t(j)/(dt*N));
end
Bk(N+1) = Bk(N+1)/(2*N);
tau = t(end)+dt - t(1); % Time interval [s]
fk = kp/tau; % Vector of frequencies

% Determine power spectrum:
P = 0.25*(Ak.^2 + Bk.^2)/N;
    P = P/max(P); % Normalize vector

% Plot:
if strcmpi(opt,'plot')
    figure(1) % Sampling plot
        if T 
            plot(t_plot,y_plot,'k-','linewidth',1.5) % Actual function
            hold on
        end
        plot(t,f,'ro','linewidth',1.5) % Sampled points

        grid on
        xlabel('Time [s]')
        ylabel('f(t)')
        title('DFT Sampling')
        if T
            legend('Actula Function','Sampled Points','location','best')
        end
        
    figure(2) % Ak plots
        stem(fk, Ak,'k','linewidth',1.5,'markerfacecolor','k')

        grid on
        xlabel('Frequency [Hz]')
        ylabel('A_k')
        title('A_k vs Frequency')

    figure(3) %  Bk plots
        stem(fk,Bk,'k','linewidth',1.5,'markerfacecolor','k')

        grid on
        xlabel('Frequency [Hz]')
        ylabel('B_k')
        title('B_k vs Frequency')
    
    figure(4) % Power spectrum
        stem(fk, P,'k','linewidth',1.5) 
        
        grid on
        xlabel('Frequency [Hz]')
        ylabel('Power (normalized)')
        title('Power Spectrum')
end