clear; close all; addpath('aux_functions'); run('plot_settings.m');
N = 50; %Number of panels  defining the airfoil
delta_s = 0.1;
s_final = 5; %amount of periods of the unsteadiness
studyCase = 'steady';

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
        alpha_vec = -15:1:20;
        lift_vec = nan(size(alpha_vec)); %pre-allocation
        for ii = 1:length(alpha_vec)
            theta =@(s) deg2rad(alpha_vec(ii)); %angle of attack    
            s = s_final;
            wake.VPos = nan(Ns,2); 
            wake.circ = nan(Ns,1);
            [panel,wake] = getCirculation(s,delta_s,s_final,N,U_inf,c,theta,studyCase,panel,wake);
            total_circ = sum(panel.circ); %Total circulation around the airfoil
            lift_vec(ii) = rho*U_inf*total_circ;
        end
    case 'unsteady'
        K = 0.1; %0.02; 0.05, 0.1 %reduced frequency K = omega*c/2/U_inf
        omega = K*2*U_inf/c; %[rad/s]
        theta =@(s) 0.5*sin(s*2*pi*omega); %[rad] pitching motion
        s = 0;
        panel.previousCirc = zeros(N,1); % We start at AOA=0° so there is initialy no circulation around the airfoil
        wake.VPos = nan(Ns,2); 
        wake.circ = nan(Ns,1);
end

%% PLOTS
figure(); xlabel('$\alpha$ [deg]'); ylabel('$C_L$ [-]'); grid on; hold on
plot(alpha_vec,lift_vec,'displayName','in-house model');
lift_theory=@(alpha) 2*pi.*deg2rad(alpha);
plot(alpha_vec,lift_theory(alpha_vec),'displayName','$C_{L} =2 \pi \alpha$');
EXP = importExperimentalData('.\experimental_data\Pelletier-mueller-2000.csv');
plot(EXP.alpha,EXP.Cl,'ko','displayName','Experiment');
legend('location','best');

%% EXTRA FUNCTIONS
function [panel,wake] = getCirculation(s,delta_s,s_final,N,U_inf,c,theta,studyCase,panel,wake)
    count = 0; %Iteration counter
    while s<=s_final
        count=count+1;
        s = s+delta_s;
        % A and B matrices 
        A = buildA(N,panel);
        B = buildB(U_inf,theta(s),studyCase,panel,wake);
        panel.circ = A\B;
        
        if strcmp(studyCase,'unsteady')
            % The change of the total circulation of the airfoil release an oposite vortex at c1/2 behind thetrailing edge
            wake.VPos(count,:) = [c 0];
            wake.circ(count) = - (sum(panel.circ)-sum(panel.previousCirc));

            [Vx,Vy,P] = VelocityPressureField(c,rho,U_inf,theta(s),panel,wake,'plot_on','save_on',count);

            % Convect wake vortices 
            wake.VPos = wake.VPos + [cos(theta(s)) , -sin(theta(s))].*U_inf*delta_s;

            wake.previousCirc = wake.circ;
        end
    end
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
function PelletierMueller2000 = importExperimentalData(filename, dataLines)
    if nargin < 2
        dataLines = [1, Inf];
    end
    % Setup the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 2);
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ",";
    % Specify column names and types
    opts.VariableNames = ["alpha", "Cl"];
    opts.VariableTypes = ["double", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Import the data
    PelletierMueller2000 = readtable(filename, opts);
end