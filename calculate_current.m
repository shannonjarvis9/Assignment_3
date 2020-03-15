function [Current, V_solution, Ex, Ey] = calculate_current(L,W, Lb, Wb, sigma_inside, sigma_outside)
% Function: calculate_current
% ---------------------------
% Inputs: 
%     L: length of region
%     W: width of region
%     Lb: location of the bottleneck center on length axis
%     Wb: location of the bottleneck center on width axis
%     sigma_inside: the value of sigma inside the bottleneck
%     sigma_outside: the value of sigma outside the bottleneck
% 
% Output: 
%     current: Calculated current 
%     



G = sparse(W*L,W*L);       %Initialize the G matrix 
F = sparse(W*L,1);         %Initialize the F matrix 
sigma = zeros(L, W);       % Initialize the sigma matrix  

delta_x = 1;     % Define parameters for FD method 
delta_y = 1;
V0 = 5;          % Define boundary voltage value 


%-------------------------------------------------------------------------- 
% Fill the sigma matrix:
%--------------------------------------------------------------------------


for i=1:L
    for j=1:W
        
        %In in the bottleneck region 
        in_x = logical( i >= (L- Lb)/2 && i <= (L + Lb)/2);
        in_y = logical( j <= Wb | j >= (W-Wb));

        if in_x && in_y
            sigma(i,j) = sigma_outside;
        else 
            sigma(i,j) = sigma_inside;
        end 
        
        
    end
end

%Indexing to plot the voltage 
[X,Y] = meshgrid(1:W,1:L);

% %-------------------------------------------------------------------------- 
% % Plot the sigma:
% %-------------------------------------------------------------------------- 
figure('Name','Question 2: Conductivity');  %Optional
surf(X,Y,sigma);
grid;
colorbar;
title('Conductivity Map', 'FontSize',20);
xlabel('y (distance)','FontSize',13);
ylabel('x (distance)','FontSize',13);
zlabel('Conductivity (\Omega^{-1})', 'FontSize',13);
view(0, 90); 



%-------------------------------------------------------------------------- 
% Fill the G matrix:
%--------------------------------------------------------------------------

%Following loop obtained from 4700Code/CondCode/GetCurrents.m
for i = 1:L
    for j = 1:W
        n = j + (i - 1) * W;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            F(n,1) = 1;
        elseif i == L
            G(n, :) = 0;
            G(n, n) = 1; %V0 = 1
        elseif j == 1
            nxm = j + (i - 2) * W;
            nxp = j + (i) * W;
            nyp = j + 1 + (i - 1) * W;

            rxm = (sigma(i, j) + sigma(i - 1, j)) / 2.0;
            rxp = (sigma(i, j) + sigma(i + 1, j)) / 2.0;
            ryp = (sigma(i, j) + sigma(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  W
            nxm = j + (i - 2) * W;
            nxp = j + (i) * W;
            nym = j - 1 + (i - 1) * W;

            rxm = (sigma(i, j) + sigma(i - 1, j)) / 2.0;
            rxp = (sigma(i, j) + sigma(i + 1, j)) / 2.0;
            rym = (sigma(i, j) + sigma(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*W;
            nxp = j + (i)*W;
            nym = j-1 + (i-1)*W;
            nyp = j+1 + (i-1)*W;

            rxm = (sigma(i,j) + sigma(i-1,j))/2.0;
            rxp = (sigma(i,j) + sigma(i+1,j))/2.0;
            rym = (sigma(i,j) + sigma(i,j-1))/2.0;
            ryp = (sigma(i,j) + sigma(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end



%-------------------------------------------------------------------------- 
% Solve for the node voltage and map back to spatial domain:
%--------------------------------------------------------------------------

V = G\F; 

V_solution = zeros(L,W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        V_solution(i,j) = V(n);
    end
end
    


%Indexing to plot the voltage 
[X,Y] = meshgrid(1:W,1:L);


% %-------------------------------------------------------------------------- 
% % Plot the node voltage:
% %-------------------------------------------------------------------------- 
% figure('Name','Voltage');  %Optional
% surf(X,Y,V_solution);
% grid;
% title('Plot of V(x,y)', 'FontSize',20);
% xlabel('y (distance)','FontSize',13);
% ylabel('x (distance)','FontSize',13);
% zlabel('Voltage (A)','FontSize',13);






%-------------------------------------------------------------------------- 
% Calculate the electric field from voltage:
%--------------------------------------------------------------------------

%Following loop obtained from 4700Code/CondCode/GetCurrents.m
for i = 1:L
    for j = 1:W
        if i == 1
            Ex(i, j) = (V_solution(i + 1, j) - V_solution(i, j));
        elseif i == L
            Ex(i, j) = (V_solution(i, j) - V_solution(i - 1, j));
        else
            Ex(i, j) = (V_solution(i + 1, j) - V_solution(i - 1, j)) * 0.5;
        end
        
        if j == 1
            Ey(i, j) = (V_solution(i, j + 1) - V_solution(i, j));
        elseif j == W
            Ey(i, j) = (V_solution(i, j) - V_solution(i, j - 1));
        else
            Ey(i, j) = (V_solution(i, j + 1) - V_solution(i, j - 1)) * 0.5;
        end
    end
end

%-------------------------------------------------------------------------- 
% Calculate the current density:
%--------------------------------------------------------------------------
Ex = -Ex;
Ey = -Ey; 

eFlowx = sigma .* Ex;
eFlowy = sigma .* Ey;


%-------------------------------------------------------------------------- 
% Calculate the current:
%--------------------------------------------------------------------------
C0 = sum(eFlowx(1, :));
Cnx = sum(eFlowx(L, :));
Current = (C0 + Cnx) * 0.5;


% %-------------------------------------------------------------------------- 
% % To plot the electric field:
% %-------------------------------------------------------------------------- 
% figure('Name','Electric Field');  %Optional
% quiver(Ex', Ey');
% axis([0 L 0 W]);
% grid;
% title('Quiver plot of the electric field (V/m)', 'FontSize',20);
% xlabel('x (distance)','FontSize',13);
% ylabel('y (distance)','FontSize',13);



% %-------------------------------------------------------------------------- 
% % To plot the current density:
% %-------------------------------------------------------------------------- 
% figure('Name','Current Density');  %Optional
% quiver(eFlowx', eFlowy');
% axis([0 L 0 W]);
% grid;
% title('Quiver plot of the current density (A/m^2)', 'FontSize',20);
% xlabel('x (distance)','FontSize',13);
% ylabel('y (distance)','FontSize',13);


end

