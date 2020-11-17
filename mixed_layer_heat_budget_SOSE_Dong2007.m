%Compute mixed layer heat budget using method from Dong et al (2007)

% Tendency = dT/dt
% Surface heat flux = (Qnet-q)/rho0*Cp*h
% Horizontal advection = u*grad(T)
% Vertical entrainment  = wvel * deltaT/h


%Load dimension data
load('/data/oceans_open/racnde/SOSE_run_Abernathey/grid.mat','XC','YC','RC','RF');
%XC = longitude, YC = latitude, RC = depth at cell edge, RF = depth at cell centre
load /data/oceans_open/racnde/sose/data/5day_ave/all_dates.mat  %dd = decimal time

%Select area over which to calculate budget
ILON = find(XC(:,1)>=120 & XC(:,1)<=300); ILAT = find(YC(1,:)>=-80 & YC(1,:)<=-40);

XC = XC(ILON,1); DXG = DXG(ILON,ILAT);
YC = YC(1,ILAT); DYG = DYG(ILON,ILAT);



%% TENDENCY TERM = dT/dt
%Load temperature averaged over mixed layer
MLT    = ncread('ML_ave_properties.nc','MLT',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);   %Mixed layer average temperature

dt = diff(dd.*(24*60*60));
heat_tendency = zeros(size(MLT));
for t = 2:length(dd)-1
    
    heat_tendency(:,:,t) = (MLT(:,:,t+1) - MLT(:,:,t-1))./(dt(t-1) + dt(t));
end
clear t dt MLT


%% SURFACE FLUX TERM
%Load data for area ILON, ILAT
MLD = ncread('MLD.nc','MLD RF',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);  %Mixed layer depth
tflux  = ncread('surface_fluxes.nc','Tflux',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);    %Net heat flux
sw_ra  = ncread('surface_fluxes.nc','oceqsw',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);   %Shortwave radiation

heat_surface_flux = zeros(length(ILON),length(ILAT),length(dd));

%Constants
rho0 = 1027;    %Seawater reference density
cp = 4000;      %Specific heat capacity of sewater

%Constants defined by turbidity
R = 0.67;
gamma1 = 1;
gamma2 = 17;


for ii = 1:length(XC)
    for jj = 1:length(YC)
        for kk = 1:length(dd)
            
            q0 = sw_rad(ii,jj,kk);   %Downward surface rad flux
            
            %Downward radiative heat flux at ML base (generally <20)
            q = q0*(R*exp(-MLD(ii,jj,kk)/gamma1) + (1-R)*exp(-MLD(ii,jj,kk)/gamma2));
            
            heat_surface_flux(ii,jj,kk) = (tflux(ii,jj,kk) - q)/(rho0*cp*MLD(ii,jj,kk));
            
        end
    end
end
clear ii jj kk R cp gamma1 gamma2 q q0 tflux sw_rad rho0



%% ADVECTIVE TERM = u*dt/dx + v*dt/dy
dTdx = zeros(size(MLT));
dTdy = zeros(size(MLT));

for jj = 1:length(dd)
    
    for yy = 1:length(YC)
        dTdx(:,yy,jj) = gradient(squeeze(MLT(:,yy,jj)),x(:,yy));
    end
    
    for xx = 1:length(XC)
        dTdy(xx,:,jj) = gradient(squeeze(MLT(xx,:,jj)),DYG(1,1));
    end
end


%Load uvel, vvel
MLuvel = ncread('ML_ave_properties.nc','MLuvel',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);  %Mixed layer averaged zonal velocity
MLvvel = ncread('ML_ave_properties.nc','MLvvel',[ILON(1) ILAT(1) 1],[length(ILON) length(ILAT) length(dd)]);  %Mixed layer averaged meridional velocity

%Put uvel, vvel onto XC,YC,RC grid - to be same as MLT
uvel_at_theta = zeros(size(MLT));
vvel_at_theta = zeros(size(MLT));

for jj = 1:length(dd)
    
    uvel_at_theta(1:end-1,:,jj) = (squeeze(MLuvel(1:end-1,:,jj)) + squeeze(MLuvel(2:end,:,jj)))./2;
    vvel_at_theta(:,1:end-1,jj) = (squeeze(MLvvel(:,1:end-1,jj)) + squeeze(MLvvel(:,2:end,jj)))./2;
end

heat_advection = -(uvel_at_theta.*dTdx + vvel_at_theta.*dTdy);
clear MLuvel MLvvel jj uvel_at_theta vvel_at_theta xx yy jj



%% DIFFUSIVE TERM = kappa*u*d2t/dx2 + v*d2t/dy2
kappa = 500;    %Eddy diffusivity
d2T = zeros(length(XC),length(YC),length(dd));

for jj = 1:length(dd)
    
    d2T(:,:,jj) = divergence(y,x,dTdy(:,:,jj),dTdx(:,:,jj));
end
heat_diffusion = kappa.*d2T;
clear kappa jj dTdy dTdx d2T x y



%% ENTRAINMENT TERM = We*deltaT

theta = ncread('theta.nc','THETA',[ILON(1) ILAT(1) 1 1],[length(ILON) length(ILAT) length(RC) length(dd)]);         %Potential temperature
wvel = ncread('ocean_velocity.nc','WVEL',[ILON(1) ILAT(1) 1 1],[length(ILON) length(ILAT) length(RF) length(dd)]);  %Vertical velocity


heat_entrainment = zeros(length(XC),length(YC),length(dd));

for ii = 1:length(XC)
    for jj = 1:length(YC)
        for kk = 1:length(dd)
            
            max_MLD = find(RF==MLD(ii,jj,kk));
            We = wvel(ii,jj,max_MLD,kk);    %Entrainment velocity = velocity at base of ML
            deltaT = MLT(ii,jj,kk) - theta(ii,jj,max_MLD+1);   %DeltaT = MLT - theta in cell immediately below MLD
            
            heat_entrainment(ii,jj,kk) = -((We.*deltaT)./MLD(ii,jj,kk));
            clear We deltaT 
        end
    end
end 
clear ii jj kk max_MLD MLD theta MLT


heat_budget_sum = heat_surface_flux + heat_advection + heat_diffusion + heat_entrainment;
heat_budget_residual = heat_tendency - heat_budget_sum;