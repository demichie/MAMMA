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

zeta_grid_reverse = zN - zeta_grid;

% Mechanical stability
c = 5.0*1e6;
phi = 38*pi/180;
rho_cr = 2600;
h_g = 1800;

omz = 0*rho_mix; omt = 0*rho_mix; omr = 0*rho_mix; 
Pmin = NaN*ones(length(rho_mix),1); PminB = Pmin;
p_lith_v=0*rho_mix;  p_lith_h=0*rho_mix; Pporo = 0*rho_mix;

ap=2 * c * cos(phi);
bp = sin(phi);
q = tan(pi/4+phi/2) * tan(pi/4+phi/2);
C0 = 2 * c * cos(phi) / (1-sin(phi));

for i= 1:length(p_lith_h)
    p_lith_v(i) = GRAV * rho_cr * (ZN - zeta_grid(i));
    p_lith_h(i) = GRAV * h_g * (ZN - zeta_grid(i));
    omz(i) = p_lith_v(i);
    omr(i) = p_1(i);
    omt(i) = 2 * p_lith_h(i) - p_1(i);
    fA = 2 * p_lith_h(i); 
    fB = p_lith_v(i);
    fH = fA * fA * (4 * bp * bp - 3) + (fB * fB - fA * fB) * (4 * bp * bp - 12);
    fK = ap + bp * (fB - 2 * Pporo(i));
    fG = fK + bp * fA;    
    fC = C0 - Pporo(i)*(q-1);
    if(omz(i)>=omt(i) && omt(i)>=omr(i))
        Pmin(i,1)=max(0,(fB-fC(1))/q);
    elseif(omt(i)>=omz(i) && omz(i)>=omr(i))
        Pmin(i,1)=max(0,(fA-fC(1))/(1+q));
    elseif(omt(i)>=omr(i) && omr(i)>=omz(i))
        Pmin(i,1)=max(0,fA-fC(1)-q*fB);
    end
    if(omz(i) >= omt(i) && omt(i) >= omr(i))
        PminB(i,1) = (1/(6-2*bp*bp))*((3*fA+2*bp*fK(:,1))-sqrt(fH+12*(fK(:,1)*fK(:,1)+bp*fA*fK(:,1))));
        PminB(i,1) = max(0,PminB(i,1));
    elseif(omt(i) >= omz(i) && omz(i) >= omr(i))
        PminB(i,1) = fA/2-(1/6)*sqrt( 12*( ap + bp*(fA - 2*Pporo(i) ) ).^2 -   3*(fA-2*fB).^2 );
        PminB(i,1) = max(0,PminB(i,1));
    elseif(omt(i) >= omr(i) && omr(i) >= omz(i))
        PminB(i,1) = (1/(6-2*bp*bp))*((3*fA-2*bp*fG)-sqrt(fH+12*(fG*fG-bp*fA*fG)));
        PminB(i,1) = max(0,PminB(i,1));
    end            
end

figure();

hold all
plot(p_1,zeta_grid_reverse,'k');
plot(p_2,zeta_grid_reverse,'k:');
plot(Pmin,zeta_grid_reverse,'r--');
plot(PminB,zeta_grid_reverse,'r:');
legend('Phase 1','Phase 2','Mohr - Coulomb','Mogi - Coulomb');
title('Pressures');
xlabel('Pressure (Pa)');
ylabel('Depth (m)');
ylim([Z0,ZN]);
set(gca,'YDir','reverse')
box on;
hold all;
