clear all; close all; 

c = 1; %[m] chord
U_inf = 1; %[m/s] freestream velocity
N = 10; %Number of panels  defining the airfoil
delta_s = 0.25;
s_final = 1;
Ns = ceil(s_final/delta_s);

rho = 1.2;

studyCase = 'unsteady';

leadingEdge = -c/4;
panelLength = c/N;
panel = table();
panel.VPos = zeros(N,2); %Vortex Position on 2D plane (1/4c of the pannel)
panel.CPPos = zeros(N,2); %Control Point Position of each panel (3/4c of the pannel) on 2D plane 
panel.VPos(1,:) = [leadingEdge+panelLength/4,0]; %firstPanelVortex position
panel.CPPos(1,:) = [leadingEdge+panelLength*3/4,0]; %firstPanelVortex position
for ii = 2:N
    panel.VPos(ii,:) = [panel.VPos(ii-1,1)+panelLength,0];
    panel.CPPos(ii,:)= [panel.CPPos(ii-1,1)+panelLength,0];
end


switch studyCase
    case 'steady'
        theta =@(s) deg2rad(30); %angle of attack    
        s = s_final;
    case 'unsteady'
        K = 0.02; %0.02; 0.05, 0.1 %reduced frequency
        theta =@(s) 0.5*sin(K*2*pi*s); %[rad] pitching motion
        s = 0;
        panel.previousCirc = zeros(N,1); % We start at AOA=0° so there is initialy no circulation around the airfoil
        wake.Vpos = nan(Ns,2); 
        wake.circ = nan(Ns,1);
end

count = 0; %Iteration counter
while s<=s_final
    count=count+1;
    s = s+delta_s;
    % A and B matrices 
    A = buildA(N,panel);
    B = buildB(U_inf,theta,panel,wake);
    panel.circ = A\B(theta(s));
    
    % The change of the total circulation of the airfoil release an oposite vortex at c1/2 behind thetrailing edge
    wake.Vpos(count,:) = [c 0];
    wake.circ(count) = - (sum(panel.circ)-sum(panel.previousCirc));
    
    
    
    %% Plot velocity field
    %TODO, rewrite circprevious and convect wake vortices
    
    N_x = 20; %amount of discretization points in x direction
    N_y = 10; %amount of discretization points in x direction
    %If the interaction matrix already exist, I load it
    intMat = sprintf('./Saved_matrices/Plot_Interac_Nx_%.0f_Ny_%.0f_N_%.0f_%s.mat',N_x,N_y,N,studyCase);
    if exist(intMat,'file')
        load(intMat)
    else
        x_vec = linspace(-c,2*c,N_x); y_vec=linspace(-1.5*c,1.5*c,N_y); [X,Y]=meshgrid(x_vec,y_vec);
        Ax = nan(N_x*N_y,N+N_wake); % for x component of induced velocities
        Ay = nan(N_x*N_y,N+N_wake); % for y components of induced veloocities
        for ii = 1:N_y
            for jj = 1:N_x
                %warning('check for unsteady')
                for kk = 1:N %loop around all the panel vortices
                    gridPoint = [X(ii,jj),Y(ii,jj)]; 
                    d_vec = gridPoint - panel.VPos(kk,:); %distance between the panel vortex kk and the grid point
                    d_norm= norm(d_vec);
                    d_ang = atan2(d_vec(2),d_vec(1));
                    Ax((ii-1)*N_x+jj,kk) = 1/(2*pi*d_norm)* -sin(d_ang);
                    Ay((ii-1)*N_x+jj,kk) = 1/(2*pi*d_norm)* cos(d_ang);
                end
            end
        end
        % Save matrices
        save(intMat,'Ax','Ay','x_vec','y_vec','X','Y')
    end

    warning('change for unsteady case with all the vorticities');
    Vx = reshape((Ax*panel.circ), [N_x,N_y])' + U_inf*cos(theta(s));
    Vy = reshape((Ay*panel.circ), [N_x,N_y])' - U_inf*sin(theta(s));
    % Pressure calculation with Bernoulli
    P = 0.5*rho*U_inf^2 - 0.5*rho.*(Vx.^2+Vy.^2);

    % circulation distribution plot
    figure('Name','circulation'); bar(panel.VPos(:,1)./c,panel.circ); xlabel('x/c'); ylabel('\Gamma [m^2/s]');
    %Airfoil plot
    figure('Name','surfplot'); hold on; xlabel('x'); ylabel('y');
    contourf(x_vec,y_vec,P,20,'LineColor','none'); cb_handle = colorbar; cb_handle.Label.String = 'Pressure [Pa]';
    plot(panel.VPos(:,1),panel.VPos(:,2),'-ok')
    quiver(X,Y,Vx,Vy)
    set(gca, 'YDir', 'reverse');
end
function [A] = buildA(N,panel)
    filename = sprintf('./Saved_matrices/A/N=%.0f.mat',N);
    if exist(filename,'file')
        load(filename,'A');
    else
        A = nan(N);
        for ii = 1:N %panel vortex
            for jj = 1:N % control point
                d_vec = panel.CPPos(jj,:) - panel.VPos(ii,:);
                d_norm = norm(d_vec); %distance between vortex ii with control point jj
                d_ang = atan2(d_vec(2),d_vec(1));
                disp(d_ang)
                A(jj,ii) = 1/(2*pi*d_norm)* cos(d_ang);
            end
        end
        save(filename,'A');
    end
end
function [B] = buildB(U_inf,theta,studyCase,panel,wake)
    NCP = size(panel.CPPos,1); %Number of panels control points
    Nw = sum(~isnan(wake.circ)); % I count the number of vortices generated so far in the wake.
    B = ones(N,1).*U_inf*sin(theta);
    
    if strcmp(studyCase,'unsteady') && Nw>0
        for ii=1:NCP % Control points
            for jj =1:Nw %Wake vortices
                
                    d_vec = panel.CPpos(ii,:) - wake.VPos(jj,:); % 2D vector
                    V_abs = abs(wake.circ(jj)/2/pi/norm(d_vec));
                    V_dir = - cross(d_vec,wake.VPos(jj,:))/norm(cross(d_vec,wake.VPos(jj,:)));
                    V = V_abs.*V_dir;
                    
                    B(ii) = B(ii)-V(2); 
            end
        end
    end
end