% Tendency = dS/dt
% Surface salt flux = Qs/rho0*h
% Horizontal advection = u*grad(S)
% Vertical entrainment  = dh/dt * deltaS/h 
% Lateral induction = u*grad(h) * deltaS/h 
% Vertical diffusion = K/h*dS/dz 


%Load dimension data
load /data/oceans_open/racnde/sose/data/5day_ave/all_dates.mat
load('/data/oceans_open/racnde/sose/data/5day_ave/grid.mat','XC','YC')

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



%% TENDENCY TERM = dS/dt
%Load salinity averaged over mixed layer
MLS = ncread('/data/oceans_open/racnde/SOSE_run_Abernathey/ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLS');

dt = diff(dd.*(24*60*60));
salt_tendency = zeros(size(MLS));
for t = 2:length(dd_all)-1
    
    salt_tendency(:,:,t)  = (MLS(:,:,t+1) - MLS(:,:,t-1))./(dt(t-1) + dt(t));
end
clear count dt MLS




%% SURFACE FLUX TERM = Qfw/rho0*h
salt_surface_flux = nan(length(XC),length(YC),438);

for t = 1:length(dd)
     
    %Load MLD
    MLD = ncread('MLD.nc','MLD RF',[1 1 t],[length(XC) length(YC) 1]);
    
    %Load surface salt flux term
    Qs = ncread('surface_fluxes.nc','Sflux',[1 1 t],[length(XC) length(YC) 1]);
       
    %Constants
    rho0 = 1035;    %Seawater reference density
    
    %Surface fw term calculated from salt flux
    salt_surface_flux(:,:,t) = Qs./(rho0.*MLD);
    
    clear Qs MLD
end




%% ADVECTION TERM = u*grad(S)
salt_advection = nan(length(XC),length(YC),438);

for t = 1:length(dd)
    
    %Mixed layer average salinity    
    MLS = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLS',[1 1 t],[length(XC) length(YC) 1]);
    
    %Mixed layer average ocean velocity
    MLuvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLuvel',[1 1 t],[length(XC) length(YC) 1]);
    MLvvel = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLvvel',[1 1 t],[length(XC) length(YC) 1]);
        
    dSdx = nan(length(XC),length(YC));
    dSdy = nan(length(XC),length(YC));
    
    for jj = 1:length(YC)
        dSdx(2:end,jj) = diff(squeeze(MLS(:,jj)))./DXC(2:end,jj);
    end
    
    for ii = 1:length(XC)
        dSdy(ii,2:end) = diff(squeeze(MLS(ii,:)))./DYC(ii,2:end);
    end
    clear MLS
    
    %Put uvel, vvel onto XC,YC,RC grid - to be same as salinity
    uvel_at_salinity(1:end-1,:) = (squeeze(MLuvel(1:end-1,:)) + squeeze(MLuvel(2:end,:)))./2;
    dSdx_at_salinity(1:end-1,:) = (squeeze(dSdx(1:end-1,:)) + squeeze(dSdx(2:end,:)))./2;
    
    vvel_at_salinity(:,1:end-1) = (squeeze(MLvvel(:,1:end-1)) + squeeze(MLvvel(:,2:end)))./2;
    dSdy_at_salinity(:,1:end-1) = (squeeze(dSdy(:,1:end-1)) + squeeze(dSdy(:,2:end)))./2;
    clear MLuvel MLvvel dSdx dSdy
    
    salt_advection(:,:,t) = -(uvel_at_salinity.*dSdx_at_salinity + vvel_at_salinity.*dSdy_at_salinity);
    
    clear uvel_at_salinity vvel_at_salinity dSdy_at_salinity dSdx_at_salinity
end





%% ENTRAINMENT TERM = dh/dt * deltaS/h

%Salt budget entrainment terms
salt_entrainment = nan(length(XC),length(YC),438);
salt_lat_ind = nan(length(XC),length(YC),438);

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
    
    salinity = ncread('salinity.nc','SALT',[1 1 1 t],[length(XC) length(YC) 42 1]);
    MLS = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLS',[1 1 t],[length(XC) length(YC) 1]);
        
    deltaS = nan(length(XC),length(YC));
    
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            
            if ~isnan(MLD(ii,jj))
                
                RF_index = find(RF==MLD(ii,jj));
                
                if ~isempty(RF_index)
                    
                    S_MLD = salinity(ii,jj,RF_index); %Salinity below MLD
                    deltaS(ii,jj) = MLS(ii,jj) - S_MLD; %Difference between ave ML salinity & salinity below
                end
            end
        end
    end
end
clear uvel vvel salinity S_MLD RF_index MLS

salt_entrainment(:,:,t) = -dhdt(:,:,t).*(deltaS./MLD);





%% LATERAL INDUCTION TERM = u*grad(h) * deltaS/h

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

salt_lat_ind(:,:,t) = u_gradh.*(deltaS./MLD);
clear MLD u_gradh deltaS






%% VERTICAL DIFFUSION TERM = K/h*dS/dz

%Vertical diffusivity coefficient
load('Kz_mask'); %Kz = 0.0450e-3 in SIZ, 0.135e-3 south of SAF, 0.05e-3 north of SAF

%Salt diffusion term
salt_diffusion   = nan(length(XC),length(YC),length(dd));

for t = 1:length(dd)
    
    %MLD
    MLD = ncread('MLD.nc','MLD RF',[1 1 t],[length(XC) length(YC) 1]);
    
    %Salinity
    salinity = ncread('salinity.nc','SALT',[1 1 t],[length(XC) length(YC) 42 1]);
    %Mixed layer average salinity
    MLS = ncread('./ML_ave_properties/ML_ave_properties_pacific_whole_SO_full_cells.nc','MLS',[1 1 t],[length(XC) length(YC) 1]);    
    
    dS = nan(length(XC),length(YC));
    dz = nan(length(XC),length(YC));
    
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            
            if ~isnan(MLD(ii,jj))
                
                if ~isempty(RF_index)
                    dS(ii,jj) = MLS(ii,jj) - salinity(ii,jj,RF_index);  %Difference between ave ML salinity & salinity below
                    
                    %Distance between midpoint of ML and 50m below ML
                    dz(ii,jj) = RC(RF_index) - (MLD(ii,jj)./2);
                end
                clear RF_index
            end
        end
    end
    clear salinity MLS ii jj
        
    %Diffusion = 1/H * Kz * dS/dz
    salt_diffusion(:,:,t) = -Kz./MLD.*(dS./dz);
end



%% Calculate residual in budget - Sum of budget terms equal to tendency if no residual
salt_sum = salt_surface_flux + salt_advection + salt_eddy_advection + salt_entrainment + salt_horizontal_diffusion + salt_vertical_diffusion;

%Residual/error in budget
salt_residual = salt_tendency_term - salt_sum;
