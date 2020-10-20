% Tendency = dT/dt
% Surface heat flux = (Qnet-q)/rho0*Cp*h
% Horizontal advection = u*grad(T)
% Vertical entrainment  = dh/dt * deltaT/h
% Lateral induction = u*grad(h) * deltaT/h
% Vertical diffusion = K/h*dT/dz


%Load dimension data
load /data/oceans_open/racnde/sose/data/5day_ave/all_dates.mat
load('/data/oceans_open/racnde/sose/data/5day_ave/grid.mat','XC','YC','RC','RF','DRC','DRF','DXC','DYC')

%Check depths are positive
if sum(RC)<0
    RC = -RC;
end
if sum(RF)<0
    RF = -RF;
end

XC = XC(:,1); YC = YC(1,:);

%XC  = longitude at cell centre
%YC  = latitude at cell centre
%RC  = depth at cell centre
%RF  = depth at cell edge
%DXC = longitudinal distance between cell centres (m)
%DYC = longitudinal distance between cell edges(m)
%DXG = latitudinal distance between cell centres (m)
%DYG = latitudinal distance between cell edges (m)
%dd  = decimal date


years = 1992:2015;

%Check RC, RF are positive
if sum(RC)<0
    RC = -RC;
end
if sum(RF)<0
    RF = -RF;
end



%% TENDENCY TERM = dT/dt
%Load temperature averaged over mixed layer
MLT = ncread('/data/oceans_open/racnde/SOSE_run_Abernathey/ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLT');

dt = diff(dd.*(24*60*60));
heat_tendency = zeros(size(MLT));
for t = 2:length(dd_all)-1
    
    heat_tendency(:,:,t) = (MLT(:,:,t+1) - MLT(:,:,t-1))./(dt(t-1) + dt(t));
end
clear count dt MLT




%% SURFACE FLUX TERM = (Qnet-q)/rho0*Cp*h
heat_surface_flux = nan(length(XC),length(YC),438);

for t = 1:length(dd)
    
    %Load MLD
    MLD = ncread('MLD.nc','MLD RF',[1 1 t],[length(XC) length(YC) 1]);
    
    %Load net surface heat flux
    Qnet = ncread('surface_fluxes.nc','Tflux',[1 1 t],[length(XC) length(YC) 1]);
    
    %Load shortwave radiation term
    q_sw = ncread('surface_fluxes.nc','oceqsw',[1 1 t],[length(XC) length(YC) 1]);
    
    
    %Constants
    rho0 = 1035;    %Seawater reference density
    Cp = 3994;      %Specific heat capacity of sewater
    
    fac = 1 ./ (Cp*rho0.*MLD);
    surface_heating_term = fac.*Qnet;
    
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            
            if ~isnan(MLD(ii,jj))
                
                RF_index = find(RF==MLD(ii,jj));   %Depth equal to MLD
                
                %Proportion of SW radiation entering each cell
                SWFRAC =      [   0.7695        0.0975      0.0600   0.0349       0.0192 ...
                    0.0104        0.0050      0.0022   8.6887e-04   2.9287e-04 ...
                    8.3983e-05    2.5737e-05 ];
                SWFRAC(1,13:42) = 0;
                
                %Surface heating term
                heat_surface_flux(ii,jj,t) = surface_heating_term(ii,jj) - (q_sw(ii,jj)*nansum(SWFRAC(RF_index+1:end)))./(Cp*rho0*MLD(ii,jj));
            end
        end
    end
    
    clear Qs Qnet q_sw fac SWFRAC MLD
end




%% ADVECTION TERM = u*grad(T)
heat_advection = nan(length(XC),length(YC),438);

for t = 1:length(dd)
    
    %Mixed layer average temperature
    MLT = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLT',[1 1 t],[length(XC) length(YC) 1]);
    
    %Mixed layer average ocean velocity
    MLuvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLuvel',[1 1 t],[length(XC) length(YC) 1]);
    MLvvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLvvel',[1 1 t],[length(XC) length(YC) 1]);
    
    
    %Calculate temperature gradients
    dTdx = nan(length(XC),length(YC));
    dTdy = nan(length(XC),length(YC));
    
    for jj = 1:length(YC)
        dTdx(2:end,jj) = diff(squeeze(MLT(:,jj)))./DXC(2:end,jj);
    end
    
    for ii = 1:length(XC)
        dTdy(ii,2:end) = diff(squeeze(MLT(ii,:)))./DYC(ii,2:end);
    end
    clear MLT
    
    %Put uvel, vvel onto XC,YC,RC grid - to be same as theta
    uvel_at_theta(1:end-1,:) = (squeeze(MLuvel(1:end-1,:)) + squeeze(MLuvel(2:end,:)))./2;
    dTdx_at_theta(1:end-1,:) = (squeeze(dTdx(1:end-1,:)) + squeeze(dTdx(2:end,:)))./2;
    
    vvel_at_theta(:,1:end-1) = (squeeze(MLvvel(:,1:end-1)) + squeeze(MLvvel(:,2:end)))./2;
    dTdy_at_theta(:,1:end-1) = (squeeze(dTdy(:,1:end-1)) + squeeze(dTdy(:,2:end)))./2;
    clear MLuvel MLvvel dTdy dTdy
    
    heat_advection(:,:,t) = -(uvel_at_theta.*dTdx_at_theta + vvel_at_theta.*dTdy_at_theta);
    
    clear uvel_at_theta vvel_at_theta dTdy_at_theta dTdx_at_theta
end





%% ENTRAINMENT TERM = dh/dt * deltaT/h

%Heat budget entrainment terms
heat_entrainment = nan(length(XC),length(YC),438);
heat_lat_ind = nan(length(XC),length(YC),438);

%Calculate rate of change of MLD
MLD = ncread('MLD.nc','MLD RF');
dhdt = nan(size(MLD));

dt = diff(dd_all.*(24*60*60));

for ii = 2:length(dd_all)-1
    
    dhdt(:,:,ii) = (MLD(:,:,ii+1) - MLD(:,:,ii-1))./(dt(ii-1) + dt(ii));
end
clear ii dt

dhdt(dhdt<0) = 0; %dh/dt only entrains, no detrainment



for t = 1:length(dd)
    
    MLD = ncread('MLD.nc','MLD RF',[1 1 t],[length(XC) length(YC) 1]);
    
    theta = ncread('theta.nc','THETA',[1 1 1 t],[length(XC) length(YC) 42 1]);
    MLT = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLT',[1 1 t],[length(XC) length(YC) 1]);
    
    deltaT = nan(length(XC),length(YC));
    
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            
            if ~isnan(MLD(ii,jj))
                
                RF_index = find(RF==MLD(ii,jj));
                
                if ~isempty(RF_index)
                    
                    T_MLD = theta(ii,jj,RF_index);    %Temp below MLD
                    deltaT(ii,jj) = MLT(ii,jj) - T_MLD; %Difference between ave ML theta & theta below
                end
            end
        end
    end
end
clear uvel vvel theta T_MLD RF_index MLT

heat_entrainment(:,:,t) = -dhdt(:,:,t).*(deltaT./MLD);





%% LATERAL INDUCTION TERM = u*grad(h) * deltaT/h

%Calculate grad of MLD for lateral induction term
dhdx = nan(length(XC),length(YC));
dhdy = nan(length(XC),length(YC));

for jj = 1:length(YC)
    dhdx(2:end,jj) = diff(squeeze(MLD(:,jj)))./DXC(2:end,jj);
end

for ii = 1:length(XC)
    dhdy(ii,2:end) = diff(squeeze(MLD(ii,:)))./DYC(ii,2:end);
end
clear ii jj


%Put uvel, vvel onto XC,YC,RC grid - to be same as theta
MLuvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLuvel',[1 1 t],[length(XC) length(YC) 1]);
MLvvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLvvel',[1 1 t],[length(XC) length(YC) 1]);

uvel_at_theta(1:end-1,:) = (squeeze(MLuvel(1:end-1,:)) + squeeze(MLuvel(2:end,:)))./2;
dhdx_at_theta(1:end-1,:) = (squeeze(dhdx(1:end-1,:)) + squeeze(dhdx(2:end,:)))./2;

vvel_at_theta(:,1:end-1) = (squeeze(MLvvel(:,1:end-1)) + squeeze(MLvvel(:,2:end)))./2;
dhdy_at_theta(:,1:end-1) = (squeeze(dhdy(:,1:end-1)) + squeeze(dhdy(:,2:end)))./2;
clear dhdx dhdy ML_uvel ML_vvel

u_dhdx = -(uvel_at_theta.*dhdx_at_theta);
v_dhdy = -(vvel_at_theta.*dhdy_at_theta);
clear uvel_at_theta vvel_at_theta dhdx_at_theta dhdy_at_theta

u_gradh = u_dhdx + v_dhdy;
clear u_dhdx v_dhdy

heat_lat_ind(:,:,t) = u_gradh.*(deltaT./MLD);
clear MLD u_gradh deltaT





%% VERTICAL DIFFUSION TERM = K/h*dT/dz

%Vertical diffusivity coefficient
load('Kz_mask'); %Kz = 0.0450e-3 in SIZ, 0.135e-3 south of SAF, 0.05e-3 north of SAF

%Heat diffusion term
heat_diffusion   = nan(length(XC),length(YC),length(dd));


for t = 1:length(dd)
    
    %MLD
    MLD = ncread('MLD.nc','MLD RF',[1 1 t],[length(XC) length(YC) 1]);
    
    %Temperature
    theta = ncread('theta.nc','THETA',[1 1 1 t],[length(XC) length(YC) 42 1]);
    %Mixed layer average temperature
    MLT = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLT',[1 1 t],[length(XC) length(YC) 1]);
    
    
    dT = nan(length(XC),length(YC));
    dz = nan(length(XC),length(YC));
    
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            
            if ~isnan(MLD(ii,jj))
                
                if ~isempty(RF_index)
                    dT(ii,jj) = MLT(ii,jj) - theta(ii,jj,RF_index);     %Difference between ave ML theta & theta below                    
                    
                    %Distance between midpoint of ML and 50m below ML
                    dz(ii,jj) = RC(RF_index) - (MLD(ii,jj)./2);
                end
                clear RF_index
            end
        end
    end
    clear theta MLT ii jj
    
    
    %Diffusion = 1/H * Kz * dT/dz
    heat_diffusion(:,:,t) = -Kz./MLD.*(dT./dz);
end





%% Calculate residual in budget - Sum of budget terms - equal to tendency if no residual
heat_sum = heat_surface_flux + heat_advection + heat_entrainment + heat_lat_ind + heat_diffusion;

%Residual/error in budget
heat_residual = heat_tendency_term - heat_sum;

