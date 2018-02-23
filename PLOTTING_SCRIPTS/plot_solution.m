clear all
[filename, pathname] = uigetfile( '*.bak','Pick a file');

read_bak(pathname,filename)

filename = strrep(filename,'.bak','_p.std');

data = importdata(strcat(pathname,filename));

L = length(data);

if(METHOD_OF_MOMENTS_FLAG == 'T')

	data_reshaped = reshape(data(1:L)',12+N_CRY+4*N_GAS+N_COMPONENTS+2*N_MOM*N_CRY,[]);

	n_grid = size(data_reshaped,2);

	comp_cells = size(data_reshaped,2);

	subgrid_idx = round(linspace(1,n_grid,25));

	zeta_grid = data_reshaped(1,:);

	z0 = Z0;
	zN = ZN;

	alfa_2(1:N_GAS,:) = data_reshaped(6+N_GAS+1:6+N_GAS+N_GAS,:);
	alfa_1(1,:) = 1.D0 - sum(alfa_2,1);

	p_1 = data_reshaped(1+1,:);
	p_2 = data_reshaped(1+2,:);
	u_1 = data_reshaped(1+3,:);
	u_2 = data_reshaped(1+4,:);
	T =   data_reshaped(1+5,:);

	rho_1 = data_reshaped(7+2*N_GAS+2*N_CRY*N_MOM + N_COMPONENTS,:);
	rho_2 = zeros(N_GAS,comp_cells);

	for i=1:N_GAS,
	    rho_2(i,:) = data_reshaped(7+2*N_GAS+2*N_CRY*N_MOM + N_COMPONENTS+i,:);
	end

	beta = zeros(N_CRY,comp_cells);
	for i=1:N_CRY,
	    beta(i,:) = data_reshaped(6+2*N_GAS+i,:);
	end

	x_d = zeros(N_GAS,comp_cells);
	for i=1:N_GAS,
	    x_d(i,:) = data_reshaped(6 + i,:);
	end
	x_d_eq = zeros(N_GAS,comp_cells);
	for i=1:N_GAS,
	    x_d_eq(i,:) = data_reshaped(7+3*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY+i,:);
	end

	components = zeros(N_COMPONENTS,comp_cells);

	for i=1:N_COMPONENTS,
	    components(i,:) = data_reshaped(6+2*N_GAS+2*N_CRY*N_MOM+i,:);
	end

	moms = zeros(2*N_CRY*N_MOM,comp_cells);

	for i=1:2*N_CRY*N_MOM,
	    moms(i,:) = data_reshaped(6+2*N_GAS+i,:);
	end

	visc =  data_reshaped(8+4*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY,:);
	visc_melt =  data_reshaped(9+4*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY,:);
	visc_rel_crystals =  data_reshaped(10+4*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY,:);
	visc_rel_bubbles =  data_reshaped(11+4*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY,:);
	radius =data_reshaped(12+4*N_GAS+2*N_MOM*N_CRY+N_COMPONENTS+N_CRY,:);
	rho_mix = alfa_1 .* rho_1 + sum(alfa_2 .* rho_2 , 1);
	c_1 = alfa_1 .* rho_1 ./ rho_mix;
	c_2 = 1.0 - c_1;
	p_mix = alfa_1 .* p_1 + sum(alfa_2,1) .* p_2;
	u_mix = c_1 .* u_1 + c_2 .* u_2;
	u_rel = u_2 - u_1;
	mass_flow_rate = pi * radius.^2 .* ( rho_mix .* u_mix );

	beta_tot = sum(beta);
	x_d_tot = sum(x_d,1);

	FigHandle = figure;
	set(FigHandle, 'Position', [100, 100, 1049, 895]);
	set(gcf,'color','w');
	color_list=['b','g','r','c','m','y','k'];
	color_list_eq=['b:','g:','r:','c:','m:','y:','k:'];
	color_list_eq2=['b--','g--','r--','c--','m--','y--','k--'];
	zeta_grid_reverse = zN - zeta_grid;

	subplot(2,5,1)
	plot(rho_1,zeta_grid_reverse);
	title('rho_1');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(2,5,2)
	hold all

	for i=1:N_GAS,
	    plot(rho_2(i,:),zeta_grid_reverse,color_list(i));
	end
	title('rho_2');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,3)
	plot(T,zeta_grid_reverse);
	title('T');
	xlabel('K');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,4)
	hold all
	for i=1:N_GAS,
	    plot(alfa_2(i,:),zeta_grid_reverse,color_list(i));
	end
	if ( N_GAS > 1)
	    i = N_GAS + 1;
	    plot(sum(alfa_2,1),zeta_grid_reverse,color_list(i));
	    legend('H20','CO2' ,'Tot','Location','SouthEast')
	end
	title('alfa_g');
	xlim([0,1]);
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,5)
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 7 ,:)./rho_mix,zeta_grid_reverse,color_list_eq(i)); hold on;
	end
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 8 ,:)./rho_mix,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m(3)/rho');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(2,5,6)
	hold all
	for i=1:N_GAS,
	    plot(x_d(i,:),zeta_grid_reverse,color_list(i));
	    plot(x_d_eq(i,subgrid_idx),zeta_grid_reverse(subgrid_idx),...
	    color_list_eq(1+2*(i-1):2*i) );
	end
	if ( N_GAS >1)
	    i = N_GAS + 1;
	    plot(x_d_tot,zeta_grid_reverse,color_list(i));
	    plot(x_d_eq_tot(subgrid_idx),zeta_grid_reverse(subgrid_idx),...
	    color_list_eq(1+2*(i-1):2*i) );
	end
	title('x_{d}');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,7)
	semilogx(u_1,zeta_grid_reverse);
	hold all
	semilogx(u_2,zeta_grid_reverse);
	title('phase vel.');
	xlabel('m/s');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;

	subplot(2,5,8)
	hold all
	plot(p_1,zeta_grid_reverse);
	plot(p_2,zeta_grid_reverse);
	title('pressures');
	xlabel('Pa');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,9)
	plot(visc,zeta_grid_reverse);
	hold all
	plot(visc_melt,zeta_grid_reverse);
	plot(visc_rel_crystals,zeta_grid_reverse);
	plot(visc_rel_bubbles,zeta_grid_reverse,'k');
	legend('Mix','Melt','Crys','Bubbles','Location','SouthEast')
	set(gca,'XScale','log')
	title('Viscosity');
	xlabel('Pa\cdot s');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,5,10)
	plot(mass_flow_rate ,zeta_grid_reverse);
	title('Mass flow rate');
	xlabel('kg/s');
	xlimits = xlim;
	xlim([0.99*xlimits(1),1.01*xlimits(2)]);
	set(gca,'YDir','reverse')
	ylim([Z0,ZN]);
	box on;
	hold all;

	FigHandle = figure;
	set(FigHandle, 'Position', [100, 100, 1049, 895]);
	set(gcf,'color','w');
	color_list=['b','g','r','c','m','y','k'];
	color_list_eq=['b:','g:','r:','c:','m:','y:','k:'];
	color_list_eq2=['b--','g--','r--','c--','m--','y--','k--'];
	zeta_grid_reverse = zN - zeta_grid;

	subplot(1,5,1);
	for i=1:N_COMPONENTS,
	    plot(components(i,:)./rho_mix,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('components/rho');
	%xlabel('-');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(1,5,2);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 7 ,:)./rho_mix,zeta_grid_reverse,color_list_eq(i)); hold on;
	end
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 8 ,:)./rho_mix,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m(3)/rho');
	%xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(1,5,3);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 7 ,:)./moms(N_CRY*N_MOM*(i-1) + 5,:),zeta_grid_reverse,color_list_eq(i)); hold on;
	end
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 8 ,:)./moms(N_CRY*N_MOM*(i-1) + 6 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m(3)/m(2)');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(1,5,4);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 5 ,:)./moms(N_CRY*N_MOM*(i-1) + 3,:),zeta_grid_reverse,color_list_eq(i)); hold on;
	end
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 6 ,:)./moms(N_CRY*N_MOM*(i-1) + 4 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m(2)/m(1)');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(1,5,5);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 3 ,:)./moms(N_CRY*N_MOM*(i-1) + 1,:),zeta_grid_reverse,color_list_eq(i)); hold on;
	end
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 4 ,:)./moms(N_CRY*N_MOM*(i-1) + 2 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m(1)/m(0)');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;




else

	data_reshaped = reshape(data(1:L)',8+2*N_CRY+4*N_GAS+4,[]);

	n_grid = size(data_reshaped,2);

	comp_cells = size(data_reshaped,2);

	subgrid_idx = round(linspace(1,n_grid,25));

	zeta_grid = data_reshaped(1,:);

	z0 = Z0;
	zN = ZN;

	alfa_2(1:N_GAS,:) = data_reshaped(6+N_GAS+1:6+N_GAS+N_GAS,:);
	alfa_1(1,:) = 1.D0 - sum(alfa_2,1);

	p_1 = data_reshaped(1+1,:);
	p_2 = data_reshaped(1+2,:);
	u_1 = data_reshaped(1+3,:);
	u_2 = data_reshaped(1+4,:);
	T =   data_reshaped(1+5,:);

	beta = zeros(N_CRY,comp_cells);

	for i=1:N_CRY,

	    beta(i,:) = data_reshaped(6+2*N_GAS+i,:);

	end

	x_d = zeros(N_GAS,comp_cells);
	for i=1:N_GAS,
	    x_d(i,:) = data_reshaped(1+5+i,:);
	end

	rho_1 = data_reshaped(7+2*N_GAS+N_CRY,:);
	rho_2 = zeros(N_GAS,comp_cells);
	for i=1:N_GAS,
	    rho_2(i,:) = data_reshaped(7+2*N_GAS+N_CRY+i,:);
	end

	beta_eq = zeros(N_CRY,comp_cells);
	for i=1:N_CRY,
	    beta_eq(i,:) = data_reshaped(7+3*N_GAS+N_CRY+i,:);
	end

	x_d_eq = zeros(N_GAS,comp_cells);

	for i=1:N_GAS,
	    x_d_eq(i,:) = data_reshaped(7+3*N_GAS+2*N_CRY+i,:);
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

	subplot(2,6,1)
	plot(rho_1,zeta_grid_reverse);
	title('rho_1');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('depth (m)')
	box on;
	hold all;

	subplot(2,6,2)
	hold all

	for i=1:N_GAS,

	    plot(rho_2(i,:),zeta_grid_reverse,color_list(i));

	end

	title('rho_2');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,3)
	plot(T,zeta_grid_reverse);
	title('T');
	xlabel('K');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,4)
	hold all
	for i=1:N_GAS,
	    plot(alfa_2(i,:),zeta_grid_reverse,color_list(i));
	end
	if ( N_GAS > 1)
	    i = N_GAS + 1;
	    plot(sum(alfa_2,1),zeta_grid_reverse,color_list(i));
	    legend('H20','CO2' ,'Tot','Location','SouthEast')
	end
	title('alfa_g');
	xlim([0,1]);
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;
	
	subplot(2,6,5)
	hold all

	for i=1:N_CRY,
	     plot(beta(i,:),zeta_grid_reverse,color_list(i));
	     plot(beta_eq(i,subgrid_idx),zeta_grid_reverse(subgrid_idx),...
	     color_list_eq(1+2*(i-1):2*i));
	end

	if ( N_CRY > 1)
	    i = N_CRY + 1;
	    plot(beta_tot,zeta_grid_reverse,color_list(i));
	    plot(beta_eq_tot(subgrid_idx),zeta_grid_reverse(subgrid_idx),...
    	    color_list_eq(1+2*(i-1):2*i));
	end
	title('beta');
	xlim([0,1]);
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;
	
	subplot(2,6,6)
	hold all
	
	for i=1:N_GAS,
	    plot(x_d(i,:),zeta_grid_reverse,color_list(i));
	    plot(x_d_eq(i,subgrid_idx),zeta_grid_reverse(subgrid_idx),...
	    color_list_eq(1+2*(i-1):2*i) );
	end

	if ( N_GAS >1)
	    i = N_GAS + 1;
	    plot(x_d_tot,zeta_grid_reverse,color_list(i));
	    plot(x_d_eq_tot(subgrid_idx),zeta_grid_reverse(subgrid_idx),...
	    color_list_eq(1+2*(i-1):2*i) );
	end
	title('x_{d}');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;
	
	subplot(2,6,7)
	semilogx(u_1,zeta_grid_reverse);
	hold all
	semilogx(u_2,zeta_grid_reverse);
	semilogx(u_mix(subgrid_idx),zeta_grid_reverse(subgrid_idx),'ro');
	title('phase vel.');
	xlabel('m/s');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	
	subplot(2,6,8)
	if ( min(u_rel) > 0 )
	    semilogx(u_rel,zeta_grid_reverse);
	else
	    plot(u_rel,zeta_grid_reverse);  
	end
	title('u_{gas} - u_{liq}');
	xlabel('m/s');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,9)
	hold all
	plot(p_1,zeta_grid_reverse);
	plot(p_2,zeta_grid_reverse);
	plot(p_mix(subgrid_idx),zeta_grid_reverse(subgrid_idx),'ro');
	title('pressures');
	xlabel('Pa');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,10)
	plot(visc,zeta_grid_reverse);
	hold all
	plot(visc_melt,zeta_grid_reverse);
	plot(visc_rel_crystals,zeta_grid_reverse);
	plot(visc_rel_bubbles,zeta_grid_reverse,'k');
	legend('Mix','Melt','Crys','Bubbles','Location','SouthEast')
	set(gca,'XScale','log')
	title('Viscosity');
	xlabel('Pa\cdot s');
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,11)
	hp(11) = plot(sum(rho_2.*alfa_2,1),zeta_grid_reverse);
	title('ex.gas bulk density');
	xlabel('kg/m^3');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(2,6,12)
	plot(mass_flow_rate ,zeta_grid_reverse);
	title('Mass flow rate');
	xlabel('kg/s');
	xlimits = xlim;
	xlim([0.99*xlimits(1),1.01*xlimits(2)]);
	set(gca,'YDir','reverse')
	ylim([Z0,ZN]);
	box on;
	hold all;

end
