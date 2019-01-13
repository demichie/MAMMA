clear all
[filename, pathname] = uigetfile( '*.bak','Pick a file');
close all
read_bak(pathname,filename)

csf = [0.8 0.7];

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

	FigHandle = figure;
	set(FigHandle, 'Position', [100, 100, 1049, 895]);
	set(gcf,'color','w');
	color_list=['k','r','g','c','m','y','b'];
	color_list_eq=['k:','r:','g:','c:','m:','y:','b:'];
	color_list_eq2=['k--','r--','g--','c--','m--','y--','b--'];
	zeta_grid_reverse = zN - zeta_grid;
        ZN = ZN./1000;
	zeta_grid_reverse=zeta_grid_reverse./1000;

	subplot(3,6,1)
	plot(rho_1,zeta_grid_reverse,'k');
	title('Phase 1 density');
	xlabel('[kg/m^3]');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]')
	box on;
	hold all;

	subplot(3,6,2)
	hold all
	for i=1:N_GAS,
	    plot(rho_2(i,:),zeta_grid_reverse,color_list(i));
	end
	title('Phase 2 density');
	xlabel('[kg/m^3]');
	ylabel('Depth [km]')
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(3,6,3)
	plot(T,zeta_grid_reverse,'k');
	title('Temperature');
	xlabel('[K]');
	ylabel('Depth [km]')
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(3,6,4)
	hold all
	for i=1:N_GAS,
	    plot(100.*alfa_2(i,:),zeta_grid_reverse,color_list(i));
	end
	if (N_GAS > 1)
	    i = N_GAS + 1;
	    plot(100.*sum(alfa_2,1),zeta_grid_reverse,color_list(i));
	    legend('H20','CO2' ,'Tot','Location','SouthEast')
	end
	title('Exsolved gas');
	xlim([0,100]);
	ylim([Z0,ZN]);
	ylabel('Depth [km]')
	xlabel('[vol.%]');
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(3,6,5)
	plot(100.*beta(1,:),zeta_grid_reverse, 'b'); hold on;
	plot(100.*beta(2,:),zeta_grid_reverse, 'r'); hold on;
	plot(100.*beta(1,:) + 100.*beta(2,:),zeta_grid_reverse, 'k'); hold on;
	title('Crystallinity');
	xlabel('[vol.%]');
        legend('Pl','Cpx','Total');
	ylim([Z0,ZN]);
        xlim([0,100])
	set(gca,'YDir','reverse')
	ylabel('Depth [km]')
	box on;
	hold all;

	subplot(3,6,6)
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
	title('Dissolved gas in melt');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Actual value','Equilibrium');
	xlabel('[wt.%]')
	box on;
	hold all;

	subplot(3,6,7)
	semilogx(u_1,zeta_grid_reverse,'k');
	hold all
	semilogx(u_2,zeta_grid_reverse,'k:');
	title('Velocity');
	xlabel('[m/s]');
	ylabel('Depth [km]');
        legend('Phase 1','Phase 2');
	ylim([Z0,ZN]);
	set(gca,'YDir','reverse')
	box on;

	subplot(3,6,8)
	hold all
	plot(p_1,zeta_grid_reverse,'k');
	plot(p_2,zeta_grid_reverse,'k:');
	title('Pressure');
	xlabel('[Pa]');
	ylabel('Depth [km]');
	ylim([Z0,ZN]);
        legend('Phase 1','Phase 2');
	set(gca,'YDir','reverse')
	box on;
	hold all;

	subplot(3,6,9)
	plot(mass_flow_rate ,zeta_grid_reverse,'k');
	title('Mass flow rate');
	xlabel('[kg/s]');
	xlimits = xlim;
	ylabel('Depth [km]');
	xlim([0.8*xlimits(1),1.2*xlimits(2)]);
	set(gca,'YDir','reverse')
	ylim([Z0,ZN]);
	box on;
	hold all;

	subplot(3,6,10);
	for i=1:N_COMPONENTS,
	    plot(100.*components(i,:)./rho_mix,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('Components/\rho');
	%xlabel('-');
	ylim([Z0,ZN]);
        xlim([0,100]);
	xlabel('[wt.%]');
	set(gca,'YDir','reverse')
	ylabel('Depth [m]');
	box on;
	hold all;
        legend('Ab','An','Di');

	subplot(3,6,11);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 7 ,:)./moms(N_CRY*N_MOM*(i-1) + 5,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(3)}/m^{(2)} (Mcx)');
	xlabel('[m]');
	ylim([Z0,ZN]);
        legend('Pl','Cpx');
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
	box on;
	hold all;

	subplot(3,6,12);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 8 ,:)./moms(N_CRY*N_MOM*(i-1) + 6 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(3)}/m^{(2)} (Pcx)');
	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Pl','Cpx');
	box on;
	hold all;

	subplot(3,6,13);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) +5 ,:)./moms(N_CRY*N_MOM*(i-1) + 3,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(2)}/m^{(1)} (Mcx)');

	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Pl','Cpx');
	box on;
	hold all;

	subplot(3,6,14);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) + 6,:)./moms(N_CRY*N_MOM*(i-1) + 4 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(2)}/m^{(1)} (Pcx)');
	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Pl','Cpx');
	box on;
	hold all;


	subplot(3,6,15);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) +3 ,:)./moms(N_CRY*N_MOM*(i-1) + 1,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(1)}/m^{(0)} (Mcx)');

	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Pl','Cpx');
	box on;
	hold all;

	subplot(3,6,16);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) +4,:)./moms(N_CRY*N_MOM*(i-1) +2 ,:),zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(1)}/m^{(0)} (Pcx)');
	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse');
        legend('Pl','Cpx');
	ylabel('Depth [km]');
	box on;
	hold all;

	subplot(3,6,17);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) +1 ,:).*u_1.*radius.^2,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(0)}uR^2 (Mcx)');

	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse')
	ylabel('Depth [km]');
        legend('Pl','Cpx');
	box on;
	hold all;

	subplot(3,6,18);
	for i=1:N_CRY,
	    plot(moms(N_CRY*N_MOM*(i-1) +2,:).*u_1.*radius.^2,zeta_grid_reverse,color_list(i)); hold on;
	end
	title('m^{(0)}uR^2 (Pcx)');
	xlabel('[m]');
	ylim([Z0,ZN]);
        xlimits = xlim;
        xlim([0 xlimits(2).*1.1]);
	set(gca,'YDir','reverse');
        legend('Pl','Cpx');
	ylabel('Depth [km]');
	box on;
	hold all;

end
