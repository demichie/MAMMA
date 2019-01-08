clear all
[filename, pathname] = uigetfile( '*.bak','Pick a file');

read_bak(pathname,filename)

filename = strrep(filename,'.bak','_p.std');

data = importdata(strcat(pathname,filename));

L = length(data);

if(METHOD_OF_MOMENTS_FLAG == 'T')

	data_reshaped = reshape(data(1:L)',12 + 2 * N_CRY + 4 * N_GAS + N_COMPONENTS + 2 * N_MOM * N_CRY ,[]);

	n_grid = size(data_reshaped,2);

	comp_cells = size(data_reshaped,2);

	subgrid_idx = round(linspace(1,n_grid,25));

	zeta_grid = data_reshaped(1,:);

	z0 = Z0;

	zN = ZN;

	alfa_2(1:N_GAS,:) = data_reshaped( 6 + N_GAS + 1 : 6 + N_GAS + N_GAS,:);

	alfa_1(1,:) = 1.D0 - sum(alfa_2,1);

	p_1 = data_reshaped(2,:);

	p_2 = data_reshaped(3,:);

	u_1 = data_reshaped(4,:);

	u_2 = data_reshaped(5,:);

	T =   data_reshaped(6,:);

	rho_1 = data_reshaped(1 + 6 + 2 * N_GAS + 2 * N_CRY * N_MOM + N_COMPONENTS , :);

	rho_2 = zeros( N_GAS , comp_cells );

	for i=1:N_GAS,

	    rho_2(i,:) = data_reshaped(1 + 6 + 2 * N_GAS + 2 * N_CRY * N_MOM + N_COMPONENTS + i,:);

	end

	x_d = zeros(N_GAS,comp_cells);

	for i=1:N_GAS,

	    x_d(i,:) = data_reshaped( 6 + i,:);

	end

	x_d_eq = zeros( N_GAS , comp_cells );

	for i=1:N_GAS,

	    x_d_eq(i,:) = data_reshaped( 7 + 3 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY + i,:);

	end

	components = zeros(N_COMPONENTS,comp_cells);

	for i=1:N_COMPONENTS,

	    components(i,:) = data_reshaped( 6 + 2 * N_GAS + 2 * N_CRY * N_MOM + i,:);

	end

	moms = zeros(2 * N_CRY * N_MOM , comp_cells);

	for i=1:2*N_CRY*N_MOM,

	    moms(i,:) = data_reshaped(6 + 2 * N_GAS + i,:);

	end

	beta = zeros(N_CRY,comp_cells);

	for i=1:N_CRY,

	    beta(i,:) =  data_reshaped(11 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY + i,:);

	end

	visc =  data_reshaped(8 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY,:);

	visc_melt =  data_reshaped( 9 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY,:);

	visc_rel_crystals =  data_reshaped(10 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY,:);

	visc_rel_bubbles =  data_reshaped(11 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY,:);

	radius = data_reshaped(12 + 4 * N_GAS + 2 * N_MOM * N_CRY + N_COMPONENTS + N_CRY + N_CRY,:);

	rho_mix = alfa_1 .* rho_1 + sum(alfa_2 .* rho_2 , 1);

	c_1 = alfa_1 .* rho_1 ./ rho_mix;

	c_2 = 1.0 - c_1;

	p_mix = alfa_1 .* p_1 + sum(alfa_2,1) .* p_2;

	u_mix = c_1 .* u_1 + c_2 .* u_2;

	u_rel = u_2 - u_1;

	mass_flow_rate = pi * radius.^2 .* ( rho_mix .* u_mix );

	beta_tot = sum(beta);

	x_d_tot = sum(x_d,1);

	x_d_eq_tot = sum(x_d_eq,1);

	color_list = ['b','r','k','g','c','m','y'];

	color_list_eq = ['b--','r--','k--','g--','c--','m--','y--'];

	zeta_grid_reverse = zN - zeta_grid;

else

	data_reshaped = reshape(data(1:L)',12 + 3 * N_CRY + 4 * N_GAS, [] );

	n_grid = size(data_reshaped, 2);

	comp_cells = size(data_reshaped, 2);

	subgrid_idx = round( linspace(1, n_grid , 25) );

	zeta_grid = data_reshaped(1,:);

	z0 = Z0;
	zN = ZN;

	alfa_2(1:N_GAS,:) = data_reshaped(6 + N_GAS + 1: 6 + N_GAS + N_GAS,:);

	alfa_1(1,:) = 1.D0 - sum(alfa_2,1);

	p_1 = data_reshaped(2,:);

	p_2 = data_reshaped(3,:);

	u_1 = data_reshaped(4,:);

	u_2 = data_reshaped(5,:);

	T =   data_reshaped(6,:);

	beta = zeros( N_CRY , comp_cells);

	for i=1:N_CRY,

	    beta(i,:) = data_reshaped( 6 + 2 * N_GAS + i,:);

	end

	x_d = zeros( N_GAS , comp_cells );

	for i=1:N_GAS,

	    x_d(i,:) = data_reshaped( 6 + i,:);

	end

	rho_1 = data_reshaped(7 + 2*N_GAS + N_CRY,:);

	rho_2 = zeros(N_GAS,comp_cells);

	for i=1:N_GAS,

	    rho_2(i,:) = data_reshaped( 7 + 2 * N_GAS + N_CRY + i , : );

	end

	beta_eq = zeros(N_CRY,comp_cells);

	for i=1:N_CRY,

	    beta_eq(i,:) = data_reshaped( 7 + 3 * N_GAS + N_CRY + i,:);

	end

	x_d_eq = zeros(N_GAS,comp_cells);

	for i=1:N_GAS,

	    x_d_eq(i,:) = data_reshaped( 7 + 3 * N_GAS + 2 * N_CRY + i,:);

	end

	visc =  data_reshaped( 7 + 4 * N_GAS + 2 * N_CRY + 1,:);

	visc_melt = data_reshaped( 7 + 4 * N_GAS + 2 * N_CRY + 2,:);

	visc_rel_crystals = data_reshaped( 7 + 4 * N_GAS + 2 * N_CRY + 3,:);

	visc_rel_bubbles =  data_reshaped( 7 + 4 * N_GAS + 2 * N_CRY + 4,:);

	radius =data_reshaped( 12 + 4 * N_GAS + 3 * N_CRY ,:);
	
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
 
	color_list = ['b','r','k','g','c','m','y'];

	color_list_eq = ['b--','r--','k--','g--','c--','m--','y--'];

	zeta_grid_reverse = zN - zeta_grid;

end

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 895]);

set(gcf,'color','w');

color_list=['b','g','r','c','m','y','k'];
color_list_eq=['bo','go','ro','co','mo','yo','ko'];
color_list_eq2=['b--','g--','r--','c--','m--','y--','k--'];

zeta_grid_reverse = zN - zeta_grid;

subplot(3,2,1)
semilogx(p_1,zeta_grid_reverse,'k');
title('(a)');
xlabel('Pressure (Pa)');
ylim([-50,ZN]);
ylabel('Depth (m)');
set(gca,'YDir','reverse')
xlim([1e5,1e9]);

subplot(3,2,2)
semilogx(p_1,alfa_2,'k');
title('(b)');
xlabel('Pressure (Pa)');
ylim([0,1]);
xlim([1e5,1e9]);
ylabel('Gas volume fraction');

subplot(3,2,3)
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

subplot(3,2,4)
loglog(p_1, k1, 'k');
ylabel('Permeability (m^2)')
title('(d)');
ylim([1e-17 1e-9]);
xlabel('Pressure (Pa)');
box on;
hold all;

subplot(3,2,5)
loglog(p_1, k2, 'k');
ylabel('Inertial Permeability (m)')
title('(e)');
ylim([1e-14 1e-6]);
xlabel('Pressure (Pa)');
box on;
hold all;
