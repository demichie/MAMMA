clear all
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

subplot(2,4,1)
plot(T,zeta_grid_reverse,'k');
title('(a)');
xlabel('Temperature (K)');
ylabel('Depth (m)');
ylim([-50,ZN]);
set(gca,'YDir','reverse')
box on;
hold all;

subplot(2,4,5)
hold all

for i=1:N_GAS,
    plot(alfa_2(i,:),zeta_grid_reverse,color_list(i));
end

if ( N_GAS > 1)
    i = N_GAS + 1;
    plot(sum(alfa_2,1),zeta_grid_reverse,'r');
    legend('Type 1','Type 2' ,'Total','Location','SouthEast')
    
end

title('(e)');
ylim([Z0,ZN]);
set(gca,'YDir','reverse')
ylim([-50,ZN]);
xlabel('Volatile volume fractions');
ylabel('Depth (m)');
box on;
hold all;

subplot(2,4,3)
hold all

for i=1:N_CRY,
    plot(beta(i,:),zeta_grid_reverse,color_list(i));    
end

if ( N_CRY > 1)
    i = N_CRY + 1; 
    plot(beta_tot,zeta_grid_reverse,'k');
end

if(N_CRY==3)
    legend('Type 1','Type 2','Type 3','Total');
end
 
title('(c)');
xlabel('Crystal content (vol. %)');
ylim([-50,ZN]);
ylabel('Depth (m)');
set(gca,'YDir','reverse')
box on;
hold all;

subplot(2,4,6)
hold all

for i=1:N_GAS,

    plot(x_d(i,:),zeta_grid_reverse,color_list(i));
    
end

if ( N_GAS >1)
    
    legend('Type 1','Type 2','Location','SouthEast')

end

title('(f)');
ylim([-50,ZN]);
set(gca,'YDir','reverse')
xlabel('Dissolved volatiles content');
ylabel('Depth (m)');
box on;
hold all;

subplot(2,4,7)
semilogx(u_mix,zeta_grid_reverse,'k');
title('(g)');
xlabel('Mixture velocity (m/s)');
ylim([-50,ZN]);
set(gca,'YDir','reverse')
ylabel('Depth (m)');
box on;

subplot(2,4,8)

if ( min(u_rel) > 0 )

    semilogx(u_rel,zeta_grid_reverse,'k');
else
    plot(u_rel,zeta_grid_reverse,'k');
    
end
    
title('(h)');
xlabel('Relative velocity (m/s)');
ylim([-50,ZN]);
set(gca,'YDir','reverse')
ylabel('Depth (m)');
box on;
hold all;

subplot(2,4,2)
hold all
plot(p_1,zeta_grid_reverse,'k');
plot(p_2,zeta_grid_reverse,'k:');
plot(p_mix(subgrid_idx),zeta_grid_reverse(subgrid_idx),'r');
title('(b)');
xlabel('Pressure (Pa)');
ylim([-50,ZN]);
ylabel('Depth (m)');
set(gca,'YDir','reverse')
box on;
legend('p_1','p_2','p_{mix}');
hold all;

subplot(2,4,4)
plot(visc,zeta_grid_reverse,'r');
hold all
plot(visc_melt,zeta_grid_reverse,'b');
plot(visc_rel_crystals,zeta_grid_reverse,'g');
legend('Liquid viscosity','Melt viscosity','Relative viscosity','Location','SouthEast')
set(gca,'XScale','log')
title('(d)');
xlabel('Viscosity (Pa\cdot s)');
set(gca,'YDir','reverse')
ylim([-50,ZN]);
box on;
hold all;
