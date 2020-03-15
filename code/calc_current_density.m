function [J_dens, calc_fc_current] = calc_current_density(Lb, Wb)
%     Lb: location of the bottleneck center on length axis
%     Wb: location of the bottleneck center on width axis
% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 
C.q = 1.60217662e-19;               % Charge of electron




sigma_inside  = 1;          %Initialize resistivity inside the bottleneck
sigma_outside = 10e-2;      %Initialize resistivity outside the bottleneck


[Current, V_solution, Ex, Ey] = calculate_current(200,100, Lb, Wb, sigma_inside, sigma_outside);

calc_fc_current = Current; 

[X,Y] = meshgrid(1:100,1:200);


%-------------------------------------------------------------------------- 
% Using electric field to calculate particle trajectory:
%-------------------------------------------------------------------------- 
e_conc = 1e15;                      % electron concentration /cm2
%Need to convert e_conc to /m2
e_conc = e_conc/0.001;


x_max = 200e-9;                     %maximum x dimension
y_max = 100e-9;                     %maximum y dimension

x1 = x_max/2 - Lb;
x2 = x_max/2 + Lb;

y1 = y_max/2 - Wb;
y2 = y_max/2 + Wb;

delta_t = .01e-12;                  % Time step 
numSteps = 50;                     % Number of time stpes 
numAtoms = 3000;                       % Number of particles 
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
    while (y(i)<= y1 || y(i)>= y2) && (x(i)>= x1 && x(i) <= x2)
        %Particle is in the box, generate new position 
        x(i) = x_max*rand(); 
        y(i) = y_max*rand(); 
    end
end 


prob_scatter = 1- exp(-delta_t/0.2e-12);


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
    
    
        %Add bottleneck conditions 
     if specular == 1 
         for j=1:numAtoms
            if (y(j)<= y1 || y(j) >= y2) && (x(j)>= x1 && x(j) <= x2)
                Vx(j) = - Vx(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
            if (y(j) <= y1 && y(j) >= y2) && (x(j) >= x1 && x(j) <= x2)
                Vy(j) = - Vy(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
         end
     else  %diffuse/rethermalize 
            for j=1:numAtoms
                while (y(j)<= y1 || y(j)>= y2) && (x(j)>= x1 && x(j) <= x2)
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
    
    %Move electron  
%     Vx = Vx + ax*delta_t;
%     Vy = Vy + ay*delta_t;
    
    
    x = x + Vx*delta_t;
    y = y + Vy*delta_t;

        
             
    
%     mfp = tau*Vavg;
%     %average time between colisions 
%     tau_calc = (mfp*numAtoms)/Vavg;
    Vavg = mean(Vx.^2 + Vy.^2); %it is already squared 

    J_dens_vect(i) = sqrt(Vavg)*C.q*e_conc;

    
end

J_dens = mean(J_dens_vect);


end

