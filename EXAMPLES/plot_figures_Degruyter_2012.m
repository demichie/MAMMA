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

%frag_eff = data_reshaped(1+N_GAS+5+N_CRY+N_GAS+1,:);

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

subplot(1,3,1)
semilogx(p_1,zeta_grid_reverse,'k');
hold on;
semilogx(p_2,zeta_grid_reverse,'k:');
semilogx(p_mix,zeta_grid_reverse,'r');
title('(a)');
xlabel('Pressure (Pa)');
ylim([-50,ZN]);
ylabel('Depth (m)');
set(gca,'YDir','reverse')
xlim([1e5,1e9]);
legend('p_1','p_2','p_{mix}');
hold all;

subplot(1,3,2)
semilogx(p_1,alfa_2,'k');
hold on;
semilogx(p_2,alfa_2,'k:');
semilogx(p_mix,alfa_2,'r');
title('(b)');
xlabel('Pressure (Pa)');
ylim([0,1]);
legend('p_1','p_2','p_{mix}');
hold all;
xlim([1e5,1e9]);

subplot(1,3,3)
loglog(p_1,u_mix,'k');
hold on;
loglog(p_2,u_mix,'k:');
loglog(p_mix,u_mix,'r');
title('(c)');
xlabel('Pressure (Pa)');
ylabel('Mixture velocity (m/s)');
xlim([1e5,1e9]);
legend('p_1','p_2','p_{mix}');
ylim([1e-3,1e3]);
