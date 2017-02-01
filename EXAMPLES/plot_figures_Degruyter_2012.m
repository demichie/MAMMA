clear all; close all;
[filename, pathname] = uigetfile( '*.bak','Pick a file');

read_bak(pathname,filename)

filename = strrep(filename,'.bak','_p.std');

data = importdata(strcat(pathname,filename));

L = length(data);
data_reshaped = reshape(data(1:L)',8+2*N_CRY+4*N_GAS+4,[]);

n_grid = size(data_reshaped,2);

comp_cells = size(data_reshaped,2);

subgrid_idx = round(linspace(1,n_grid,25));

zeta_grid = data_reshaped(1,:);

z0 = Z0;
zN = ZN;

alfa_2(1:N_GAS,:) = data_reshaped(1+1:1+N_GAS,:);
alfa_1(1,:) = 1.D0 - sum(alfa_2,1);

p_1 = data_reshaped(1+N_GAS+1,:);
p_2 = data_reshaped(1+N_GAS+2,:);
u_1 = data_reshaped(1+N_GAS+3,:);
u_2 = data_reshaped(1+N_GAS+4,:);
T =   data_reshaped(1+N_GAS+5,:);

beta = zeros(N_CRY,comp_cells);

for i=1:N_CRY,

    beta(i,:) = data_reshaped(1+N_GAS+5+i,:);

end

x_d = zeros(N_GAS,comp_cells);

for i=1:N_GAS,

    x_d(i,:) = data_reshaped(1+N_GAS+5+N_CRY+i,:);

end

rho_1 = data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1,:);

rho_2 = zeros(N_GAS,comp_cells);

for i=1:N_GAS,

    rho_2(i,:) = data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+i,:);

end

beta_eq = zeros(N_CRY,comp_cells);

for i=1:N_CRY,

    beta_eq(i,:) = data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+i,:);

end

x_d_eq = zeros(N_GAS,comp_cells);

for i=1:N_GAS,

    x_d_eq(i,:) = data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+N_CRY+i,:);

end

visc =  data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+N_CRY+N_GAS+1,:);

visc_melt =  data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+N_CRY+N_GAS+2,:);

visc_rel_crystals =  data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+N_CRY+N_GAS+3,:);

visc_rel_bubbles =  data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1+N_GAS+N_CRY+N_GAS+4,:);

radius = data_reshaped(8+2*N_CRY+4*N_GAS+4,:);

rho_mix = alfa_1 .* rho_1 + sum(alfa_2 .* rho_2 , 1);

c_1 = alfa_1 .* rho_1 ./ rho_mix;
c_2 = 1.0 - c_1;

p_mix = alfa_1 .* p_1 + sum(alfa_2,1) .* p_2;

u_mix = c_1 .* u_1 + c_2 .* u_2;
u_rel = u_2 - u_1;

mass_flow_rate = pi * radius.^2 .* ( rho_mix .* u_mix );

beta_tot = sum(beta);

x_d_tot = sum(x_d,1);

beta_eq_tot = sum(beta_eq,1);
 
x_d_eq_tot = sum(x_d_eq,1);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 895]);

set(gcf,'color','w');

color_list=['b','g','r','c','m','y','k'];
color_list_eq=['bo','go','ro','co','mo','yo','ko'];
color_list_eq2=['b--','g--','r--','c--','m--','y--','k--'];

zeta_grid_reverse = zN - zeta_grid;

subplot(2,2,1)
semilogx(p_1,zeta_grid_reverse,'k');
title('(a)');
xlabel('Pressure (Pa)');
ylim([-50,ZN]);
ylabel('Depth (m)');
set(gca,'YDir','reverse')
xlim([1e5,1e9]);

subplot(2,2,2)
semilogx(p_1,alfa_2,'k');
title('(b)');
xlabel('Pressure (Pa)');
ylim([0,1]);
xlim([1e5,1e9]);

subplot(2,2,3)
loglog(p_1,u_1,'k');
hold on;
loglog(p_1,u_2,'k--');
title('(c)');
xlabel('Pressure (Pa)');
ylabel('Velocity (m/s)');
xlim([1e5,1e9]);
legend('u_{melt}','u_{gas}');
ylim([1e-3,1e3]);

radius_bubble = ( (1.0 - alfa_1) ./ ( 4.0 / 3.0 * pi *  10.^LOG10_BUBBLE_NUMBER_DENSITY .* ...
    ( alfa_1 ) ) ).^( 1.0 / 3.0 ) ;

throat_radius = radius_bubble * THROAT_BUBBLE_RATIO;

k1 = 0.1250 * throat_radius.^2.0 .* (1.0 - alfa_1).^TORTUOSITY_FACTOR;

k2 = throat_radius ./ FRICTION_COEFFICIENT .* (1.0 - alfa_1).^ ...
    ( ( 1.0 + 3.0 * TORTUOSITY_FACTOR ) ./ 2.0 ) ;

subplot(2,2,4)
[hAx,hLine1,hLine2] = plotyy(p_1, k1, p_1, k2, @loglog);
set(hLine1,'LineStyle','-');
set(hLine2,'LineStyle','--');

set(hLine1,'color','k');
set(hLine2,'color','k');

set(hAx,{'ycolor'},{'k';'k'}) 

set(hAx(1),'YLim',[1e-17 1e-9]) 
set(hAx(2),'YLim',[1e-14 1e-6]) 

ylabel(hAx(1),'Permeability (m^2)') 
ylabel(hAx(2),'Inertial Permeability (m)') 
legend('Permeability','Inertial Permeability')
box on;
hold all;
