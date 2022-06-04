clear all; close all; 
N = 10; %Number of panels  defining the airfoil
delta_s = 0.1;
s_final = 5; %amount of periods of the unsteadiness
studyCase = 'unsteady';

c = 1; %[m] chord
U_inf = 1; %[m/s] freestream velocity
Ns = ceil(s_final/delta_s);
rho = 1.2;
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
        K = 0.1; %0.02; 0.05, 0.1 %reduced frequency K = omega*c/2/U_inf
        omega = K*2*U_inf/c; %[rad/s]
        theta =@(s) 0.5*sin(s*2*pi*omega); %[rad] pitching motion
        s = 0;
        panel.previousCirc = zeros(N,1); % We start at AOA=0° so there is initialy no circulation around the airfoil
        wake.VPos = nan(Ns,2); 
        wake.circ = nan(Ns,1);
end

count = 0; %Iteration counter
while s<=s_final
    count=count+1;
    s = s+delta_s;
    % A and B matrices 
    A = buildA(N,panel);
    B = buildB(U_inf,theta(s),studyCase,panel,wake);
    panel.circ = A\B;
    
    % The change of the total circulation of the airfoil release an oposite vortex at c1/2 behind thetrailing edge
    wake.VPos(count,:) = [c 0];
    wake.circ(count) = - (sum(panel.circ)-sum(panel.previousCirc));
    
    [Vx,Vy,P] = VelocityPressureField(c,rho,U_inf,theta(s),panel,wake,'plot_on','save_on',count);
    
    % Convect wake vortices 
    wake.VPos = wake.VPos + [cos(theta(s)) , -sin(theta(s))].*U_inf*delta_s;
    
    wake.previousCirc = wake.circ;
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
    Nw = sum(~isnan(wake.circ)); % Number of vortices generated so far in the wake.
    
    B = ones(NCP,1).*U_inf*sin(theta);
    
    if strcmp(studyCase,'unsteady') && Nw>0
        for ii=1:NCP % Control points
            for jj =1:Nw %Wake vortices
                
                    d_vec = panel.CPPos(ii,:) - wake.VPos(jj,:); % 2D vector
                    V_abs = abs(wake.circ(jj)/2/pi/norm(d_vec));
                    V_dir = - cross([d_vec 0],[0 0 wake.circ(jj,:)])/norm(cross([d_vec 0],[0 0 wake.circ(jj,:)]));
                    V = V_abs.*V_dir;
                    
                    B(ii) = B(ii)-V(2); 
            end
        end
    end
end
function [Vx,Vy,P] = VelocityPressureField(c,rho,U_inf,theta,panel,wake,plot_flag,save_flag,counter)
    if strcmp(plot_flag,'plot_off'); close all; end
    N_x = 20; %amount of discretization points in x direction
    N_y = 10; %amount of discretization points in y direction

    x_vec = linspace(-c,2*c,N_x); y_vec=linspace(-1.5*c,1.5*c,N_y); [X,Y]=meshgrid(x_vec,y_vec);
    Vx = ones(N_y,N_x).* U_inf*cos(theta); 
    Vy = ones(N_y,N_x).*(-U_inf*sin(theta)); 
    % compilation of all the vortices
    vortices.circ = [panel.circ; wake.circ];
    vortices.VPos = [panel.VPos; wake.VPos];
    Nw = sum(~isnan(vortices.circ)); % Number of vortices generated so far in the wake.

    for ii = 1:N_y
        for jj = 1:N_x
            for kk = 1:Nw %loop around all the vortices
                gridPoint = [X(ii,jj),Y(ii,jj)]; 
                d_vec = gridPoint - vortices.VPos(kk,:); %distance between the panel vortex kk and the grid point

                V_abs = abs(vortices.circ(kk)/2/pi/norm(d_vec));
                V_dir = - cross([d_vec 0],[0 0 vortices.circ(kk,:)])/norm(cross([d_vec 0],[0 0 vortices.circ(kk,:)]));
                V = V_abs.*V_dir;

                Vx(ii,jj) = Vx(ii,jj) + V(1); 
                Vy(ii,jj) = Vy(ii,jj) + V(2); 
            end
        end
    end
        
    % Pressure calculation with Bernoulli
    P = 0.5*rho*U_inf^2 - 0.5*rho.*(Vx.^2+Vy.^2);
    
    % circulation distribution plot
    fig_1 = figure('Name','circulation'); bar(panel.VPos(:,1)./c,panel.circ); xlabel('x/c'); ylabel('\Gamma [m^2/s]');
    if save_flag;   saveas(gcf,sprintf('./figures/circulation/%.0f.jpg',counter));  end
    %Airfoil plot
    fig = figure('Name','surfplot'); hold on; xlabel('x'); ylabel('y');
    contourf(x_vec,y_vec,P,20,'LineColor','none'); cb_handle = colorbar; cb_handle.Label.String = 'Pressure [Pa]';
    plot(panel.VPos(:,1),panel.VPos(:,2),'-ok'); plot(wake.VPos(:,1),wake.VPos(:,2),'ro'); quiver(X,Y,Vx,Vy);
    set(gca, 'YDir', 'reverse'); title(sprintf('theta = %.2f',theta)); xlim([-1 8]);
    fig.Position = [100 100 1750 500];
    if save_flag;  mkdir('./figures/field'); saveas(gcf,sprintf('./figures/field/%.0f.jpg',counter));  end
    close(fig,fig_1)
end