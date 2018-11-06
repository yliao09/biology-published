% This is the reaction-diffusion model for the MinDE oscillator published
% in Liao and Rust, Cell Systems (2018).


% INPUTS: 

% L: cell length in um

% tEnd: duration of simulation in seconds

% x_gridSz_um: spatial grid size in um

% viszKymo: set to 1 to visualize simulation results.

% simParams: biochemical parameters for the model. 


% OUTPUTS: 

% c_D_total_timeAvg: time-average profile of total MinD (cytosolic, membrane-bound and MinE-bound) concentration

% x: spatial axis of the simulated kymograph

% t: temporal axis of the simulated kymograph

% c_Dd: kymograph of total MinD concentration

% The rest of the outputs (c_ADP, c_ATP, c_E, c_d, c_de and c_e) are
% kymographs of individual species.


% Sample calling syntax: 

% L = 5; tEnd = 60*60; x_gridSz_um = 0.1; viszKymo = 1;

% simParams = [12.5, 9.5, 0.6, 0.4, 0.05, 7.6, 0.52, 2.59, 0.0843, 0.649, 0.00935, 0.000417, 0.0457, 850,890];

% [c_D_total_timeAvg, x, t, c_Dd, c_ADP, c_ATP, c_E, c_d, c_de, c_e] = sixSpecies_YL2_MCMC(L, tEnd, x_gridSz_um, viszKymo, simParams);



function [c_D_total_timeAvg, x, t, c_Dd, c_ADP, c_ATP, c_E, c_d, c_de, c_e] = MinDE_YL_MCMC(L, tEnd, x_gridSz_um, viszKymo, simParams)

% Additional simulation parameters
xGridN = round(L/x_gridSz_um);
tGridN = 360;

% Plotting and formatting parameters
tDisplay = 300:tGridN; % Unit: grid (i.e., 1 means t = 0);
L_tick = 2; % Unit: um
labelFontSz = 15;
lWidth = 2;

% Kinetic parameters
c_ADP_init = simParams(14);
c_ATP_init = 0;
c_d_init = 0;
c_E_init = simParams(15);
c_e_init = 0;
c_de_init = 0;

D_D = simParams(1);
D_E = simParams(2);
D_d = simParams(3);
D_e = simParams(4);
D_de = simParams(5);

omega_1 = simParams(6);
omega_2 = simParams(7);
omega_3 = simParams(8);
omega_4 = simParams(9);
omega_5 = simParams(10);
omega_6 = simParams(11);
omega_7 = simParams(12);
omega_8 = simParams(13);

% Prepare for PDE solver
m = 0;
x = linspace(0, L, xGridN);
t = linspace(0, tEnd, tGridN);
sol = pdepe(m, @pde_MinC_1D, @pde_MinC_1D_ic, @pde_MinC_1D_bc, x, t);
c_ADP = sol(:, :, 1)';
c_ATP = sol(:, :, 2)';
c_d = sol(:, :, 3)';
c_E = sol(:, :, 4)';
c_e = sol(:, :, 5)';
c_de = sol(:, :, 6)';

c_Dd = c_ADP + c_ATP + c_d + c_de;
c_D_total_timeAvg = mean(c_Dd(:, tDisplay), 2);

% -----------------------------------------------------------------------
% Visualize MinDE kymographs
%
if viszKymo
    t_minutes = t/60;
    figure;
    set(gcf, 'Color', 'w'); colormap('jet');
    
    subplot(4, 2, 1);
    pcolor(t_minutes(tDisplay), x, c_ADP(:, tDisplay));
    shading interp;
    c_D_colorbar = colorbar; c_D_colorbar.Label.String = '[MinD-ADP]_{c} ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 2);
    pcolor(t_minutes(tDisplay), x, c_ATP(:, tDisplay));
    shading interp;
    c_D_colorbar = colorbar; c_D_colorbar.Label.String = '[MinD-ATP]_{c} ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 3);
    pcolor(t_minutes(tDisplay), x, c_d(:, tDisplay));
    shading interp;
    c_d_colorbar = colorbar; c_d_colorbar.Label.String = '[MinD]_{m} ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 4);
    pcolor(t_minutes(tDisplay), x, c_E(:, tDisplay));
    shading interp;
    c_E_colorbar = colorbar; c_E_colorbar.Label.String = '[MinE]_{c} ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 5);
    pcolor(t_minutes(tDisplay), x, c_e(:, tDisplay));
    shading interp;
    c_e_colorbar = colorbar; c_e_colorbar.Label.String = '[MinE]_{m} ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 6);
    pcolor(t_minutes(tDisplay), x, c_de(:, tDisplay));
    shading interp;
    c_de_colorbar = colorbar; c_de_colorbar.Label.String = '[MinDE] ({\it\mu}m^{-1})';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 7);
    pcolor(t_minutes(tDisplay), x, c_Dd(:, tDisplay));
    shading interp;
    c_Dd_colorbar = colorbar; c_Dd_colorbar.Label.String = '{\bf[MinD]_{total} ({\it\mu}m^{-1})}';
    xlabel('Time (min)', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick', 0:L_tick:L,'LineWidth', lWidth, 'FontSize', labelFontSz);
    
    subplot(4, 2, 8);
    plot(mean(c_Dd(:, tDisplay), 2), x, 'k-', 'LineWidth', lWidth);
    xlabel('<[MinD]_{total}>_t', 'FontSize', labelFontSz); ylabel('{\itL} ({\it\mu}m)', 'FontSize', labelFontSz); set(gca, 'YTick',0:L_tick:L, 'YLim', [0 L],'LineWidth', lWidth, 'FontSize', labelFontSz); grid on;
    pbaspect([1.5 1 1])
end % viszKymo
% -----------------------------------------------------------------------
% Nested functions -- parameters are provided by the outer function.
%
    function [c, f, s] = pde_MinC_1D(x,t,u,DuDx)
        c = [1; 1; 1; 1; 1; 1];
        f = [D_D; D_D; D_d; D_E; D_e; D_de].*DuDx;
        
        
        A1 = omega_1*u(1); % Nucleotide exchange
        A2 = omega_2*u(2); % Direct MinD recruitment to membrane
        A3 = omega_3*u(2)*u(3); % Cooperative MinD recruitment to membrane by membrane-bound MinD
        A4 = omega_4*u(3)*u(4); % Formation of MinDE by membrane-bound MinD + cytosolic MinE
        A5 = omega_5*u(3)*u(5); % Formation of MinDE by membrane-bound MinD + membrane-bound MinE
        A6 = omega_6*u(6); % Dissociation of MinDE to MinD-ADP + lingering MinE
        A7 = omega_7*u(6); % Dissociation of MinDE to MinD-ADP + cytosolic MinE
        A8 = omega_8*u(5); % Release of membrane-bound MinE back to cytosol
        
        
        s = [A6 + A7 - A1; ... % MinD-ADP
            A1 - A2 - A3; ... % MinD-ATP
            A2 + A3 - A4 - A5; ... % membrane-bound MinD
            A7 + A8 - A4; ... % cytosolic MinE
            A6 - A5 - A8; ... % membrane-bound MinE
            A4 + A5 - A6 - A7]; % MinDE
        
    end
    function u0 = pde_MinC_1D_ic(x)
        u0 = [c_ADP_init; c_ATP_init; c_d_init; c_E_init; c_e_init; c_de_init].*(x/L)*2.0;
    end
    function [pl, ql, pr, qr] = pde_MinC_1D_bc(xl, ul, xr, ur, t) % Neumann Boundary Conditions
        pl = [0; 0; 0; 0; 0; 0];
        ql = [1; 1; 1; 1; 1; 1];
        pr = [0; 0; 0; 0; 0; 0];
        qr = [1; 1; 1; 1; 1; 1];
    end
end % End of function 'MinDE_YL_MCMC'