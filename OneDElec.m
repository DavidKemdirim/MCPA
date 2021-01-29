% One Dimensional Electron Model 

clear
clc

a = 1; % simulates for 1 particle
% a = 5; % simulates for 5 particles at once

x = zeros(a,1);
xp = 0;
v = zeros(a,1);

t = 0;
dt = 1;

% F = 1; % default arbitrary force
F = 1e-24; % force for the electron to approach the speed of light but no exceed it

% m = 1; % default arbitrary mass
m = 9.11e-31; % electron mass in kg

steps = 1000; % time in "seconds" for the simulation

re = 0; % default rebound velocity factor
% re = -0.5; % Bounces back at half its initial speed

p = 0.05; % default probability of scattering
Tau = 20;

for i = 2:steps
    
    t(i) = t(i-1) + dt;
    
    v(:,i) = v(:,i-1) + F/m*dt;
    x(:,i) = x(:,i-1) + v(:,i-1)*dt + (F/m * dt^2)/2; 
    
    p = 1 - exp(-dt/Tau); % p is 0.0488 = 4.88% when T = 20 and dt = 1
    r = rand(a,1) < p;    
    v(r,i) = re*v(r,i);
    Vav(i,:) = mean(v,2);
    
    
    % Velocity vs Time
    subplot(3,1,1), plot(t,v,'r');
    grid on;
    hold on;
    subplot(3,1,1), plot(t,Vav,'k');
    hold off;
    xlabel('time')
    ylabel('v')
    title(['Average V: ', num2str(Vav(i,:)/1e7), 'e7 m/s'])
    
    % Velocity vs Position 
    subplot(3,1,2), plot(x(1,:),v,'b');
    grid on;
    hold on;
    subplot(3,1,2), plot(x(1,:),Vav,'k');
    hold off;    
    xlabel('x')
    ylabel('v')
        
    % Position vs Time
    subplot(3,1,3), plot(t,x,'g');
    grid on;
    xlabel('time')
    ylabel('x')     
      
    pause(0.01)
end 

display('Average V:')
Vav(steps,:)

