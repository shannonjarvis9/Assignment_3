clear all 
clear

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 
C.q = 1.60217662e-19;               % Charge of electron

e_conc = 1e15;                      % electron concentration /cm2
%Need to convert e_conc to /m2
e_conc = e_conc/0.001;


x_max = 200e-9;                     %maximum x dimension
y_max = 100e-9;                     %maximum y dimension

delta_t = .01e-12;                  % Time step 

numSteps = 150;                     % Number of time stpes 
numAtoms = 3000;                       % Number of particles 
numPlot = 25;                       % Number of particles to plot 

J_dens = zeros(numSteps,1);              %Vector for the current density
J_drift = zeros(numSteps,1);              %Vector for the current density
time = linspace(0,delta_t*numSteps, numSteps);

%Applied voltage in each direction
Applied_V_x = 0.1;
Applied_V_y = 0.1; 

E_x = Applied_V_x/x_max;
E_y = Applied_V_y/y_max;

F_x = E_x*C.q;
F_y = E_y*C.q;

ax = F_x/C.m_e;
ay = 0; %F_y/C.m_e;

e_field = sqrt(E_x^2 + E_y^2);
% Thermal Velocity = 1.870192676075498e+05
v_th = sqrt((2 * C.kb * C.temp) / C.m_e);

%mean free path = 3.740385352150996e-08
tau = 0.2e-12;
mfp= v_th*tau;


% Initialize the particle position
x = x_max*rand(numAtoms,1);     %Size of region in x-direction is 200nm 
y = y_max*rand(numAtoms,1);     %Size of region in y-direction is 100nm



%Assign a velocity from the Maxwell Boltzmann Distribution 
sigma = sqrt((C.kb * C.temp) / C.m_e);
MB = makedist('Normal', 'mu', 0, 'sigma', sigma);
Vx = icdf(MB,rand(numAtoms,1));
Vy = icdf(MB, rand(numAtoms,1));

% delta_t = (1e6)*(x_max/mean(Vx.^2 + Vy.^2));
time = linspace(0,delta_t*numSteps, numSteps);

figure(3);
hist(sqrt(Vx.^2 + Vy.^2));
title('Distribution of electron velocity');
xlabel('Velocity (m/s)');
ylabel('Frequency');

color = hsv(numPlot);

prob_scatter = 1- exp(-delta_t/0.2e-12);

%loop over time in steps, plot position of points at each time step 
%Don't loop over particles, just keep the current and prev posiiton 
for i = 1:numSteps

    %Add boundary conditions
    above_x_bounds = logical(x>=x_max);
    below_x_bounds = logical(x<=0);
    
    above_y_bounds = logical(y>=y_max);
    below_y_bounds = logical(y<=0);   

    x(above_x_bounds) = 0;
    
    x(below_x_bounds) = x_max;

    Vy(above_y_bounds) = -Vy(above_y_bounds);
    y(above_y_bounds) = y_max;
    

    Vy(below_y_bounds) = -Vy(below_y_bounds);
    y(below_y_bounds) = 0;
    
    %Rethermalize   
    for j=1:length(x)
        if prob_scatter > rand()
                Vx(j)= icdf(MB, rand());
                Vy(j)= icdf(MB, rand());
        end
    end
    

    y_prev = y;
    x_prev = x; 
    
    
    %Move electron    
    Vx = Vx + ax*delta_t;
    Vy = Vy + ay*delta_t;
    
    x = x + Vx*delta_t;
    y = y + Vy*delta_t;
    
   
    
    figure(1)
    xlabel('(m)');
    ylabel('(m)');
    title('Plot of the particle trajectory');
    axis ([0 x_max 0 y_max]);
    %pause(0.0001)
    hold on;
    
    for j=1:numPlot 
        plot([x_prev(j)';x(j)'], [y_prev(j)';y(j)'], 'color', color(j,:)); 
    end 
    
    Vavg = mean(Vx.^2 + Vy.^2); %it is already squared 
    %v_th^2 = (2 * C.kb * C.temp) / C.m_e);
    J_dens(i) = sqrt(Vavg)*C.q*e_conc;
%   J_drift(i) = sqrt(Vavg)*e_field*C.q*e_conc;
    
    mfp = tau*Vavg;
    %average time between colisions 
    tau_calc = (mfp*numAtoms)/Vavg;
    

    
end

figure(2)
plot(time,J_dens, '.');
title('Plot of the current over time');
xlabel('Time (s)');
ylabel('Drift Current (A/m)');

% figure(2)
% plot(time,Tavg);
% str = sprintf('Plot of the average temperature (average is %d K)', mean(Tavg));
% title(str);
% xlabel('Time');
% ylabel('Temperature (K)');


figure(3);
subplot(3,1,1);
hist(Vx);
title('Distribution of Vx');
xlabel('Velocity (m/s)');
ylabel('Frequency');

subplot(3,1,2); 
hist(Vy);
title('Distribution of Vy');
xlabel('Velocity (m/s)');
ylabel('Frequency');

subplot(3,1,3);
hist(sqrt(Vx.^2 + Vy.^2));
title('Distribution of electron velocity');
xlabel('Velocity (m/s)');
ylabel('Frequency');



% % Density Map
% figure(4);
% hist3([x y], 'CdataMode', 'auto');
% colorbar;
% grid on;
% title('Final particle position density map')
% xlabel('x position (m)');
% ylabel('y position (m)');
% zlabel('Particle count');
% hold off;

% Temperature Map

%Need to map average temperature to a box 
T = zeros(100,50);
NumPart = zeros(100,50);
x_bins = linspace(0,x_max, 100);
y_bins = linspace(0,y_max, 50);

for i = 1:numAtoms
    for x_idx = 1:99
        for y_idx = 1:49
            if (x(i,1) > x_bins(x_idx) && x(i,1) <= x_bins(x_idx+1) && y(i,1) > y_bins(y_idx) && y(i,1) <= y_bins(y_idx+1))
                T(x_idx,y_idx) = T(x_idx,y_idx) + ((Vx(i,1)^2 + Vy(i,1)^2)*C.m_e)/(2*C.kb);
                NumPart(x_idx,y_idx) = NumPart(x_idx,y_idx) + 1;
            end
        end
    end
end

[x_b, y_b] = meshgrid(1:2:100,1:2:200);

Temp = T./NumPart;
Temp(isnan(Temp))=0; %Nan if no particles in grid element
figure(4);
surf(x_b, y_b, Temp);
view(0,90);
title('Average Temperature Plot');
xlabel('x (nm)');
ylabel('y (nm)');
colorbar;


% Density Map
density = NumPart./sum(NumPart); 
density(isnan(density))=0; 
figure(5);
surf(x_b, y_b, density);
view(0,90);
colorbar;
grid on;
title('Final particle position density map (Normalized)')
xlabel('x position (nm)');
ylabel('y position (nm)');
colorbar;
hold off;



figure(6);
surf(x_b, y_b, NumPart);
view(0,90);
colorbar;
grid on;
title('Final particle position density map')
xlabel('x position (nm)');
ylabel('y position (nm)');
colorbar;
hold off;









































































