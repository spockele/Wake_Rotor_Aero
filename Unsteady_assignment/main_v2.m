clear; close all; addpath('aux_functions'); run('plot_settings.m');
N = 50; %Number of panels  defining the airfoil
studyCase = 'unsteady'; K=0.1; %0.02; 0.05, 0.1 %reduced frequency K = omega*c/2/U_inf

delta_s = 0.1;
s_final = 1/2/K; %final time in order to get one full period
c = 1; %[m] chord
U_inf = 1; %[m/s] freestream velocity
Ns = ceil(s_final/delta_s);
rho = 1.2; nu = 1.5e-5; %[m2/s]
Re = U_inf*c/nu;
leadingEdge = -c/4; %[m]
panelLength = c/N; %[m]
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
        alpha_vec = -20:1:20; %[deg]
        lift_vec = nan(size(alpha_vec)); textList={}; %pre-allocation
        for ii = 1:length(alpha_vec)
            theta =@(s) deg2rad(alpha_vec(ii)); %angle of attack    
            s = s_final;
            wake.VPos = nan(Ns,2); wake.circ = nan(Ns,1);
            [panel,wake] = getCirculation(s,delta_s,s_final,N,rho,U_inf,c,theta,studyCase,panel,wake,0);
            total_circ = sum(panel.circ); %Total circulation around the airfoil
            lift_vec(ii) = rho*U_inf*total_circ;
            textList{ii,1} = sprintf('./figures/field/%s/subplot_alpha_%.0f_1.jpg',studyCase,rad2deg(theta(2)));
        end
        Cl_vec = lift_vec./(0.5*rho*U_inf^2*c);
        save './figures/field/steady/subplot_figList.mat' textList
        save './figures/field/steady/Cl_vs_alpha.mat' Cl_vec alpha_vec
        % PLOT CL vs alpha
        plot_steady_Cl(alpha_vec,Cl_vec);
    case 'unsteady'
        omega = K*2*U_inf/c; %[rad/s]
        theta =@(s) deg2rad(20)*sin(s*2*pi*omega); %[rad] pitching motion
        s = 0;
        panel.previousCirc = zeros(N,1); % We start at AOA=0° so there is initialy no circulation around the airfoil
        wake.VPos = nan(Ns,2); 
        wake.circ = nan(Ns,1);
        
        [panel,wake,Cl_vec,theta_vec] = getCirculation(s,delta_s,s_final,N,rho,U_inf,c,theta,studyCase,panel,wake,K);
        
        % save Cl vs alpha
        save(sprintf('./figures/field/unsteady/K_%.2f/Cl_vs_alpha.mat',K),'Cl_vec','theta_vec');
        % plot Cl vs alpha
        plot_transient_Cl     
end

%% PLOTS
function []= plot_steady_Cl(alpha_vec,Cl_vec)
    figure(); xlabel('$\alpha$ [deg]'); ylabel('$C_L$ [-]'); grid on; hold on
    plot(alpha_vec,Cl_vec,'.-','displayName','in-house model');
    lift_theory=@(alpha) 2*pi.*deg2rad(alpha);
    plot(alpha_vec,lift_theory(alpha_vec),'displayName','$C_{L} =2 \pi \alpha$');
    EXP = importExperimentalData('.\experimental_data\Pelletier-mueller-2000.csv');
    plot(EXP.alpha,EXP.Cl,'ko','displayName','Experiment');
    legend('location','best');
    %--------------------SAVE IMAGE-----------------------------------%
    plotName = strcat('./figures/steady_Cl_vs_alpha.pdf');
    set(gcf, 'Position', 600.*[0.1 0.1 1.5 1]);
    set(gcf, 'PaperPosition', 10.*[0 0 1.5 1]); 
    set(gcf, 'PaperSize',  10.*[1.5 1]); 
    print(plotName,'-dpdf','-bestfit');
end
function []= plot_transient_Cl()
    figure(); xlabel('$\alpha$ [deg]'); ylabel('$C_L$ [-]'); grid on; hold on
    % unsteady results
    for K=[0.02, 0.05, 0.1]
        load(sprintf('./figures/field/unsteady/K_%.2f/Cl_vs_alpha.mat',K),'Cl_vec','theta_vec');
        plot(rad2deg(theta_vec),Cl_vec,'o-','displayName',sprintf('Unsteady model $K=%.2f$',K));
    end
    %steady results
    load './figures/field/steady/Cl_vs_alpha.mat' Cl_vec alpha_vec 
    plot(alpha_vec,Cl_vec,'-','displayName','Steady model');
    legend('location','best');
    %--------------------SAVE IMAGE-----------------------------------%
    plotName = strcat('./figures/unsteady_Cl_vs_alpha.pdf');
    set(gcf, 'Position', 600.*[0.1 0.1 1.5 1]);
    set(gcf, 'PaperPosition', 10.*[0 0 1.5 1]); 
    set(gcf, 'PaperSize',  10.*[1.5 1]); 
    print(plotName,'-dpdf','-bestfit');
end
%% EXTRA FUNCTIONS
function [panel,wake,Cl_vec,theta_vec] = getCirculation(s,delta_s,s_final,N,rho,U_inf,c,theta,studyCase,panel,wake,K)
    count = 0; %Iteration counter
    %pre-allocation
    Cl_vec  = nan(N,1); 
    theta_vec = nan(N,1);
    
    while s<=s_final
        count=count+1;
        s = s+delta_s;
        % A and B matrices 
        A = buildA(N,panel);
        B = buildB(U_inf,theta(s),studyCase,panel,wake);
        panel.circ = A\B;
        
        switch studyCase
            case 'steady'
                VelocityPressureField(c,rho,U_inf,theta(s),panel,wake,'plot_on','save_on',studyCase,count,K);
            case 'unsteady'
                % The change of the total circulation of the airfoil release an oposite vortex at c1/2 behind thetrailing edge
                wake.VPos(count,:) = [c 0];
                wake.circ(count) = - (sum(panel.circ)-sum(panel.previousCirc));

                VelocityPressureField(c,rho,U_inf,theta(s),panel,wake,'plot_on','save_on',studyCase,count,K);

                % Convect wake vortices 
                wake.VPos = wake.VPos + [cos(theta(s)) , -sin(theta(s))].*U_inf*delta_s;
                wake.previousCirc = wake.circ;
                % compute the Cl and save it in a matrix
                total_circ = sum(panel.circ); %Total circulation around the airfoil
                Cl_vec(count) = rho*U_inf*total_circ / (0.5*rho*U_inf^2*c);
                theta_vec(count) = theta(s);
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
                B(ii) = B(ii) - V(2);
            end
        end
    end
end
function [Vx,Vy,P] = VelocityPressureField(c,rho,U_inf,theta,panel,wake,plot_flag,save_flag,studyCase,counter,K)
    if strcmp(plot_flag,'plot_off'); close all; end
    N_x = 60; %amount of discretization points in x direction
    N_y = 10; %amount of discretization points in y direction
    
    x_vec = linspace(-c,6*c,N_x); y_vec=linspace(-1.5*c,1.5*c,N_y); [X,Y]=meshgrid(x_vec,y_vec);
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
                % De-singularization of the velocity 
                d_min = 1e-2;
                if abs(d_vec)<d_min
                    d_vec = d_min * sign(d_vec);
                end
                    
                V_abs = abs(vortices.circ(kk)/2/pi/norm(d_vec));
                if norm(cross([d_vec 0],[0 0 vortices.circ(kk,:)])) ~= 0
                    V_dir = - cross([d_vec 0],[0 0 vortices.circ(kk,:)])/norm(cross([d_vec 0],[0 0 vortices.circ(kk,:)]));
                else
                    V_dir = - cross([d_vec 0],[0 0 vortices.circ(kk,:)]);
                end
                V = V_abs.*V_dir;

                Vx(ii,jj) = Vx(ii,jj) + V(1); 
                Vy(ii,jj) = Vy(ii,jj) + V(2); 
            end
        end
    end
        
    % Pressure calculation with Bernoulli
    P = 0.5*rho*U_inf^2 - 0.5*rho.*(Vx.^2+Vy.^2);
    % SUB-PLOT
    fig=figure(); subplot(1,4,1);
    bar(panel.VPos(:,1)./c,panel.circ); ylim([-0.18 0.18]); xlabel('$x/c$ [-]'); ylabel('$\Gamma$ [m$^2$/s]');
    %Airfoil plot
    subplot(1,4,2:4); hold on; xlabel('x'); ylabel('y');
    contourf(x_vec,y_vec,P,20,'LineColor','none'); cb_handle = colorbar; 
    switch studyCase
        case 'steady'
            caxis([-1 0.3]); 
        case 'unsteady'
            caxis([-1 0.3]); 
    end
    cb_handle.Label.String = 'Pressure [Pa]';
    plot(panel.VPos(:,1),panel.VPos(:,2),'-ok'); plot(wake.VPos(:,1),wake.VPos(:,2),'ro'); quiver(X,Y,Vx,Vy);
    set(gca, 'YDir', 'reverse'); title(sprintf('$\\theta = %.1f^\\circ$',rad2deg(theta))); xlim([-1 6]);
    fig.Position = [100 100 1700 500];
    if strcmp(save_flag,'save_on') && strcmp(studyCase,'steady')
        saveas(gcf,sprintf('./figures/field/%s/subplot_alpha_%.0f_%.0f.jpg',studyCase,rad2deg(theta),counter));  
    elseif strcmp(save_flag,'save_on') && strcmp(studyCase,'unsteady')
        mkdir(sprintf('./figures/field/%s/K_%.2f',studyCase,K));
        saveas(gcf,sprintf('./figures/field/%s/K_%.2f/%.0f.jpg',studyCase,K,counter));
    end
    close(fig)
%     % circulation distribution plot
%     fig_1 = figure('Name','circulation'); bar(panel.VPos(:,1)./c,panel.circ); xlabel('$x/c$ [-]'); ylabel('$\Gamma$ [m$^2$/s]'); title(sprintf('theta = %.0f°',rad2deg(theta)));
%     if save_flag;   saveas(gcf,sprintf('./figures/circulation/%s/aplha_%.0f_%.0f.jpg',studyCase,rad2deg(theta),counter));  end
%     %Airfoil plot
%     fig = figure('Name','surfplot'); hold on; xlabel('x'); ylabel('y');
%     contourf(x_vec,y_vec,P,20,'LineColor','none'); cb_handle = colorbar; cb_handle.Label.String = 'Pressure [Pa]';
%     plot(panel.VPos(:,1),panel.VPos(:,2),'-ok'); plot(wake.VPos(:,1),wake.VPos(:,2),'ro'); quiver(X,Y,Vx,Vy);
%     set(gca, 'YDir', 'reverse'); title(sprintf('theta = %.0f°',rad2deg(theta))); xlim([-1 6]);
%     fig.Position = [100 100 1310 500];
%     if save_flag; saveas(gcf,sprintf('./figures/field/%s/alpha_%.0f_%.0f.jpg',studyCase,rad2deg(theta),counter));  end
%     close(fig,fig_1)
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