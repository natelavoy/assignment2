%% ELEC 4700 Assignment 2
%% Nathan Lavoy 
%% 100995612
%% Submitted: Feb. 22, 2019

%% Part One: Electrostatic Potential with Vo Fixed at Each End

% Simulation parameters
clear; clc;
L = 30;    % Length
W = 20;      % Width
Vo = 5;     % Fixed potential source

% using GV = F as the matrix solution
G = zeros(L*W);
F = zeros(1,L*W);
for i = 1:L
    for j = 1:W
        N = j + (i-1)*W;
        Nxm = j + (i-2)*W;
        Nxp = j + i*W;
        Nym = j-1 + (i-1)*W;
        Nyp = j+1 + (i-1)*W;
        if i == 1 
            G(N,N) = 1;
            F(N)=Vo; 
        elseif i == L
            G(N,N) = 1;
            F(N) = Vo;
        elseif j == 1
            G(N,N) = 1;
        elseif j == W
            G(N,N) = 1;
        else
            G(N,N)=-4;
            G(N,Nxm)=1;
            G(N,Nxp)=1;
            G(N,Nym)=1;
            G(N,Nyp)=1;
        end
    end
end
% Convert to sparse after allocating locations to improve efficiency 
G = sparse(G);
F = sparse(F);
Vsol = G\F';

% Reformat
Vmap = zeros(W,L);
for j = 1:W
    for i = 1:L
        n = j + (i - 1) * W;
        Vmap(j, i) = Vsol(n);
    end
end

figure(2)
surf(Vmap)
title(' Voltage Map for Numerical Solution')
view(10,45)

%%
% Numeric solution to laplace equation
V = zeros(L,W);
V(L,:) = Vo;
V(1,:) = Vo;
[X,Y] = meshgrid(1:W,1:L);  % Mesh grid for plotting
for x= (-L/2+1):L/2
    for y=1:W
        if (x == (-L/2+1) || x == (L/2))
            V(x+L/2,y) = Vo;
        elseif (y == 1 || y == W)
            V(x+L/2,y) = 0;
        else
            for N=1:2:251
                V(x+L/2,y) = V(x+L/2,y) + 4*Vo/pi*1/N*cosh(N*pi*x/W)/cosh(N*pi*0.5*L/W)*sin(N*pi*y/W);
            end
        end
    end
end

figure(2)
title('part1(b)');
surf(X,Y,V);
xlabel('X');
ylabel('Y');
%% Part One Conclusion
% The solution found in the meshing method is very similar to the solution
% developed in the analytical analysis. The problem with the analytic
% solution is that it requires a infinite sumation. I set n to 250 steps to
% get a clear value. This method slower than the meshing method (due to
% Matlab's sparse techniques) but can reproduce a more accurate result.
%% Part Two A: Current Flow
% Simulation parameters
clear; clc;
L = 30;                     % Length
W = 20;                     % Width
[X,Y] = meshgrid(1:L,1:W);  % Mesh grid for plotting
Vo = 5;                     % Fixed potential source
G = zeros(L*W,L*W);         % G matrix
F = zeros(L*W,1);           % F matrix
sigma1 = 1;                 % Sigma value 1    
sigma2 = 1e-2;              % Sigma value 2

% Solve Simulation
[S,V,J] = solveG(L,W,0.4,0.6,0.4,0.6,sigma1,sigma2,Vo);

% plot for sigma
figure(3)
surf(S);
title('Resistive Surface Plot');

% plot for voltage
figure(4);
surf(X,Y,V);
title('Surface Plot for Voltage');

% Electric fields plots
[Ex, Ey] = gradient(V);
Ex = -Ex;
Ey = -Ey;

figure(5);
surf(Ex);
title('Surface Plot of X-Component Electric Field');

figure(6);
surf(Ey);
title('Surface Plot of Y-Component Electric Field');

figure (7)
quiver(X,Y,Ex, Ey);
title('Electric Field Vector Plot');

% Current Densities
Jx = S.*Ex;
Jy = S.*Ey;

figure(8)
surf(J)
title('Surface Plot of Current Density');

figure(9)
quiver(X,Y,Jx, Jy);
title('Vector Plot of Current Density');
%% Part Two B: Mesh Density
% Simulation parameters
clear; clc;
Vo = 5;         % Fixed potential source
sigma1 = 1;     % Sigma value 1    
sigma2 = 1e-2;  % Sigma value 2
n = 1;          % Iteration
% Solve system using increasing mesh sizes
for i=20:10:100
    [S,V,J] = solveG((3/2)*i,i,0.4,0.6,0.4,0.6,sigma1,sigma2,Vo);
    C_columns = sum(J, 1);               
    C(n) = sum(C_columns);  % Current at size
    ms(n) = i;              % Mesh Size
    n = n+1;                % Iterate
end
% Plot total current vs mesh size    
figure(10)
plot(ms,C);
title('Current vs. Mesh Size');
%% Part Two C: Narrowing Bottle Neck
% Simulation parameters
clear; clc;
L = 30;                     % Length
W = 20;                     % Width
Vo = 5;                     % Fixed potential source
sigma1 = 1;                 % Sigma value 1    
sigma2 = 1e-2;              % Sigma value 2
n = 1;      
for i=0.1:0.01:0.9
    Lower = 0.5-(i/2);
    Upper = 0.5+(i/2);
    [S,V,J] = solveG(L,W,0.4,0.6,Lower,Upper,sigma1,sigma2,Vo);
    C_columns = sum(J, 1);  % sum of current density is current               
    C(n) = sum(C_columns);  % need to get a scalar
    Gap(n) = i;
    n = n+1;
end
    
figure(11)
plot(Gap,C);
title('Current vs. Bottle Neck')
%% Part Two D: Resistivity in Bottle Neck
% Simulation parameters
clear; clc;
L = 30;                     % Length
W = 20;                     % Width
Vo = 5;                     % Fixed potential source
sigma1 = 1;                 % Sigma value 1    
sigma2 = 1e-2;              % Sigma value 2
n = 1;      
for i=sigma2:0.01:1
    [S,V,J] = solveG(L,W,0.4,0.6,0.4,0.6,sigma1,i,Vo);
    C_columns = sum(J, 1);  % sum of current density is current               
    C(n) = sum(C_columns);  % need to get a scalar
    rs(n) = i;
    n = n+1;
end
figure(12)
plot(rs,C);
title('Current vs. sigma')
%% Additional Functions Used
function [S,V,J] = solveG(L,W,Lb1,Lb2,Wb1,Wb2,sigma1,sigma2,Vo)
    G = zeros(L*W); % G matrix
    F = zeros(1,L*W);   % F matrix
    S = zeros(W,L);    % Mapped sigmas
    % Fill a 2D array for the sigma surface plot
    for x = 1:L
        for y = 1:W
            if x >= Lb1*L && x <= Lb2*L && (y <= Wb1*W || y >= Wb2*W)
                S(y,x) = sigma2;
            else
                S(y,x) = sigma1;
            end
        end
    end
    for i = 1:L
        for j = 1:W
            N = j + (i-1)*W;
            Nxm = j + (i-2)*W;
            Nxp = j + i*W;
            Nym = j-1+(i-1)*W;
            Nyp = j+1+(i-1)*W;
            % Fixed current end (left) 
            if i == 1
                G(N, N) = 1;
                F(N) = Vo;
            % Right boundary
            elseif i == L
                G(N, N) = 1;
            % Lower boundary
            elseif j == 1  
                % Lower boundary in high resistive area
                if i > Lb1*L && i < Lb2*L 
                    G(N, N) = -3;
                    G(N, Nyp) = sigma2;
                    G(N, Nxp) = sigma2;
                    G(N, Nxm) = sigma2;
                % Lower boundary in low resistive area
                else
                    G(N, N) = -3;
                    G(N, Nyp) = sigma1;
                    G(N, Nxp) = sigma1;
                    G(N, Nxm) = sigma1;
                end
            % Upper Boundary    
            elseif j == W
                % Upper boundary in high resistive area
                if i > Lb1*L && i < Lb2*L 
                    G(N, N) = -3;
                    G(N, Nym) = sigma2;
                    G(N, Nxp) = sigma2;
                    G(N, Nxm) = sigma2;
                % Upper boundary in low resistive area
                else
                    G(N, N) = -3;
                    G(N, Nym) = sigma1;
                    G(N, Nxp) = sigma1;
                    G(N, Nxm) = sigma1;
                end
            % Inner area    
            else
                if i > Lb1*L && i < Lb2*L && (j < Wb1*W||j > Wb2*W)
                    % High resistive area
                    G(N, N) = -4;
                    G(N, Nyp) = sigma2;
                    G(N, Nym) = sigma2;
                    G(N, Nxp) = sigma2;
                    G(N, Nxm) = sigma2;
                else
                    % Low resistive area
                    G(N, N) = -4;
                    G(N, Nyp) = sigma1;
                    G(N, Nym) = sigma1;
                    G(N, Nxp) = sigma1;
                    G(N, Nxm) = sigma1;
                end
            end
        end
    end
    
    % Solve and reformat 
    G = sparse(G);
    F = sparse(F);
    Vsol = G\F';
    % Reformat
    V = zeros(W,L);
    for j = 1:W
        for i = 1:L
            n = j + (i - 1) * W;
            V(j, i) = Vsol(n);
        end
    end
    
    % Electric fields 
    [Ex, Ey] = gradient(V);
    Ex = -Ex;
    Ey = -Ey;
    
    % Current Densities
    Jx = S.*Ex;
    Jy = S.*Ey;
    J = sqrt(Jx.^2 + Jy.^2);
end
%% Part Two Conclusion
% The results from part 2a match what I was expecting from the simulation.
% Using the relative sensitivity plot as a baseline, the voltage is
% constant until the bottle neck where there is some leakage through the
% middle. This is also seen in the electric field vector plot as well as 
% the current density where it behaves according to the Maxwell 
% relationship. The current increases as the mesh size increases. This
% value should converge as the system should get more realisistic with the
% smaller squares. As the gap increased the current increased indicating 
% the overall resistance decreased. As the barrier conductivity increased 
% the current. Overall I believe my simulation has some bugs but overall
% conveys a realisitic simulation for the requirements. 
% increased. This is due to the overall decrease in resistance. 