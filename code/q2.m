clear all 
clear

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 
C.q = 1.60217662e-19;               % Charge of electron





sigma_inside  = 1;          %Initialize resistivity inside the bottleneck
sigma_outside = 10e-2;      %Initialize resistivity outside the bottleneck


[Current, V_solution, Ex, Ey] = calculate_current(200,100, 40, 40, sigma_inside, sigma_outside);


[X,Y] = meshgrid(1:100,1:200);

%-------------------------------------------------------------------------- 
% Plot the node voltage:
%-------------------------------------------------------------------------- 
figure('Name','Voltage');  %Optional
surf(X,Y,V_solution);
grid;
title('Plot of V(x,y)', 'FontSize',20);
xlabel('y (nm)','FontSize',13);
ylabel('x (nm)','FontSize',13);
zlabel('Voltage (A)','FontSize',13);


%-------------------------------------------------------------------------- 
% To plot the electric field:
%-------------------------------------------------------------------------- 
figure('Name','Electric Field');  %Optional
quiver(Ex', Ey');
axis([0 200 0 100]);
grid;
title('Electric field (V/m) vector plot', 'FontSize',20);
xlabel('x (nm)','FontSize',13);
ylabel('y (nm)','FontSize',13);


%-------------------------------------------------------------------------- 
% Using electric field to calculate particle trajectory:
%-------------------------------------------------------------------------- 
e_conc = 1e15;                      % electron concentration /cm^2
e_conc = e_conc/0.001;              % electron concentration /m^2

x_max = 200e-9;                     %maximum x dimension
y_max = 100e-9;                     %maximum y dimension

delta_t = .01e-12;                  % Time step 
numSteps = 1000;                     % Number of time stpes 
numAtoms = 30000;                       % Number of particles 
numPlot = 50;                       % Number of particles to plot 

specular = 1; %BC for bottleneck boxes

T_avg = zeros(numSteps,1);
time = linspace(0,delta_t*numSteps, numSteps);


Fx = Ex*C.q;
Fy = Ey*C.q;

ax = Fx/C.m_e;
ay = Fy/C.m_e;



e_field = mean(sqrt(Ex.^2 + Ey.^2));


% Thermal Velocity = 1.870192676075498e+05
v_th = sqrt((2 * C.kb * C.temp) / C.m_e);

%Assign a velocity from the Maxwell Boltzmann Distribution 
sigma = sqrt((C.kb * C.temp) / C.m_e);
MB = makedist('Normal', 'mu', 0, 'sigma', sigma);
Vx = icdf(MB,rand(numAtoms,1));
Vy = icdf(MB, rand(numAtoms,1));


%mean free path = 3.740385352150996e-08
tau = 0.2e-12;
mfp= v_th*tau;


% Initialize the particle position
x = x_max*rand(numAtoms,1);     %Size of region in x-direction is 200nm 
y = y_max*rand(numAtoms,1);     %Size of region in y-direction is 100nm


%Check if pqarticle is in the box/not allowed region 
for i=1:numAtoms
    while (y(i)<= 40e-9 || y(i)>= 60e-9) && (x(i)>= 80e-9 && x(i) <= 120e-9)
        %Particle is in the box, generate new position 
        x(i) = x_max*rand(); 
        y(i) = y_max*rand(); 
    end
end 



color = hsv(numPlot);

prob_scatter = 1- exp(-delta_t/0.2e-12);

%loop over time in steps, plot position of points at each time step 
figure1 = figure;
% figure(10);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
hold on;
rectangle('position', [80e-9 0 40e-9 40e-9]);
rectangle('position', [80e-9 60e-9 40e-9 40e-9]);
hold on;
xlim(axes1,[0 x_max]);
ylim(axes1,[0 y_max]);
xlabel('nm');
ylabel('nm');

for i = 1:numSteps
    
    %Add boundary conditions
    above_x_bounds = logical(x>=x_max);
    below_x_bounds = logical(x<=0);
    
    above_y_bounds = logical(y>=y_max);
    below_y_bounds = logical(y<=0);
    
%     x(above_x_bounds) = x(above_x_bounds) - x_max;
%     x_prev(above_x_bounds) = 0;
    x(above_x_bounds) = 0;
    
%     x(below_x_bounds) = x(below_x_bounds) + x_max;
%     x_prev(below_x_bounds) = x_max;
    x(below_x_bounds) = x_max;

%     y(above_y_bounds) = -y(above_y_bounds) + 2*y_max;
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
    
    
        %Add bottleneck conditions 
     if specular == 1 
         for j=1:numAtoms
            if (y(j)<= 40e-9 || y(j) >= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                Vx(j) = - Vx(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
            if (y(j) <= 40e-9 && y(j) >= 60e-9) && (x(j) >= 80e-9 && x(j) <= 120e-9)
                Vy(j) = - Vy(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
         end
     else  %diffuse/rethermalize 
            for j=1:numAtoms
                while (y(j)<= 40e-9 || y(j)>= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                    %Particle is in the box, generate new position 
                    x(j) = x_max*rand(); 
                    y(j) = y_max*rand();
                    x(j) = x_prev(j);
                    y(j) = y_prev(j);
                end
            end    
     end    
         
     
    y_prev = y;
    x_prev = x; 
    
    
    
    %Find where the electron is located an apply the acceleration at that
    %position 
    x_bins = (0:1:200)*1e-9;
    y_bins = (0:1:200)*1e-9;
    
    for j = 1:numAtoms
        for x_idx = 1:200
            for y_idx = 1:100
                if (x(j,1) > x_bins(x_idx) && x(j,1) <= x_bins(x_idx+1) && y(j,1) > y_bins(y_idx) && y(j,1) <= y_bins(y_idx+1))
                        Vx(j,1) = Vx(j,1) + (ax(x_idx, y_idx)/1e-9)*delta_t;
                        Vy(j,1) = Vy(j,1) + (ay(x_idx, y_idx)/1e-9)*delta_t;
                       
                end
            end
        end
    end
    
%     %Move electron  
%     Vx = Vx + ax*delta_t;
%     Vy = Vy + ay*delta_t;
    
    x = x + Vx*delta_t;
    y = y + Vy*delta_t;

              
  
    xlabel('(m)');
    ylabel('(m)');
    title('Plot of the particle trajectory');
    axis ([0 x_max 0 y_max]);
    hold on;
    %pause(0.0001)
        
    for j=1:numPlot 
        
        plot([x_prev(j)';x(j)'], [y_prev(j)';y(j)'], 'color', color(j,:)); 
            
    end 
    
    Vavg = mean(Vx.^2 + Vy.^2); %it is already squared 
    %v_th^2 = (2 * C.kb * C.temp) / C.m_e);
    %Tavg(i) = ( Vavg*C.m_e)/(2*C.kb);
    
    J_dens(i) = sqrt(Vavg)*C.q*e_conc;
    
    mfp = tau*Vavg;
    %average time between colisions 
    tau_calc = (mfp*numAtoms)/Vavg;

    
end




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
surf(y_b, x_b, Temp);
title('Average Temperature Plot');
xlabel('x (nm)');
ylabel('y (nm)');


% Density Map
density = NumPart./sum(NumPart); 
density(isnan(density))=0; 
figure(5);
surf(y_b, x_b, NumPart);
grid on;
title('Final particle position density map')
xlabel('x position (nm)');
ylabel('y position (nm)');
hold off;






