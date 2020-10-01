% Tendency = dT/dt
% Surface heat flux = (Qnet-q)/rho0*Cp*h
% Horizontal advection (plus eddy advection) = u*grad(T)
% Vertical entrainment  = dh/dt * deltaT/h 
% Lateral induction = u*grad(h) * deltaT/h 
% Horizontal diffusiond = k*d2T/dx2 
% Vertical diffusion = K/h*dT/dz 


%Load dimension data
load('ECCOv4_grid.mat','XC','YC','RC','RF','DXC','DYC','DXG','DYG','dd')

years = 1992:2015;

%Check RC, RF are positive
if sum(RC)<0
    RC = -RC;
end
if sum(RF)<0
    RF = -RF;
end



%% TENDENCY TERM = dT/dt

MLT = ncread('ML_properties.nc','MLT');	%Mixed layer averaged temperature
dt = diff(dd.*(24*60*60));

heat_tendency_term = zeros(size(MLT));

for count = 2:length(dd)-1
    heat_tendency_term(:,:,count) = (MLT(:,:,count+1) - MLT(:,:,count-1))./(dt(count-1) + dt(count));
end
clear count dt MLT 



%% SURFACE FLUX TERM = (Qnet-q)/rho0*Cp*h

%Method adapted from ftp://ecco.jpl.nasa.gov/Version4/Release3/doc/evaluating_budgets_in_eccov4r3.pdf

heat_surface_flux = nan(length(XC),length(YC),length(dd));

for year = 1:length(years)
    
    year_index = find(dd>=datenum(years(year),01,01,00,00,00)&dd<datenum(years(year)+1,01,01,00,00,00));
    
    %MLD
    MLD = ncread('MLD.nc','MLD RF',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);

    %Net heat flux
    Qnet = ncread(strcat('oceQnet.',num2str(years(year)),'.nc'),'oceQnet');
    
    %Shortwave radiation
    Qsw = ncread(strcat('oceQsw.',num2str(years(year)),'.nc'),'oceQsw');
    
    %Set constants
    rho0 = 1029; %Reference density
    Cp = 3994;   %Specific heat capacity
    
    %Constants describing turbidity
    R = 0.62;
    zeta1 = 0.6; 
    zeta2 = 20;
    
    nLevels = 20;   %Number of depth levels
    
    q1=R*exp(1/zeta1*-RF(1:nLevels)) + (1-R)*exp(1/zeta2*-RF(1:nLevels));
    q2=R*exp(1/zeta1*-RF(2:(nLevels+1))) + (1-R)*exp(1/zeta2*-RF(2:(nLevels+1)));
    
    zCut = find(RC>200,1);  %No SW penetration below 200m in ECCOv4
    q1(zCut:nLevels) = 0;
    q2((zCut-1):nLevels)=0;
    
    %Fraction of shortwave passing through each cell of mixed layer
    SWFRAC    = zeros(length(RF));
    SWFRAC(1) = (1-(q1(1)-q2(1)));
    
    for kk = 2:nLevels
        SWFRAC(kk) = (q1(kk) - q2(kk));
    end
 
    fac = 1 ./ (Cp*rho0.*MLD);
    surface_heating_term = fac.*Qnet;
        
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            for month = 1:12
                
                MLD_index = find(RF==MLD(ii,jj,month));
                                
                if ~isnan(MLD_index)
                    %Surface heating term
                    heat_surface_flux(ii,jj,month+12*(year-1)) = surface_heating_term(ii,jj,month) - ...
                    	(Qsw(ii,jj,month)*sum(SWFRAC(MLD_index:end))./(Cp*rho0*MLD(ii,jj,month)));
                end
            end
        end
    end
    clear Qs Qsw Qnet SWFRAC MLD_RF MLD_index Cp rho0 ii jj month surface_heating_term
end



%% ADVECTION TERM = u*grad(T)

heat_advection = nan(length(XC),length(YC),length(dd));
heat_eddy_advection = nan(length(XC),length(YC),length(dd));

for year = 1:length(years)
    
    year_index = find(dd>=datenum(years(year),01,01,00,00,00)&dd<datenum(years(year)+1,01,01,00,00,00));
    
    %Mixed layer average temperature
    MLT = ncread('ML_ave_properties.nc','MLT',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    
    %Mixed layer averaged velocities
    MLuvel = ncread('ML_ave_properties.nc','MLuvel',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    MLvvel = ncread('ML_ave_properties.nc','MLvvel',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
   
    %Mixed layer averaged eddy velocities
    MLuvel_star = ncread('ML_ave_properties.nc','MLuvel_star',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    MLvvel_star = ncread('ML_ave_properties.nc','MLvvel_star',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);


    %Calculate horizontal temperature gradients
    dTdx = nan(length(XC),length(YC),12);
    dTdy = nan(length(XC),length(YC),12);
    
    for month = 1:12
        for jj = 1:length(YC)
            dTdx(2:end,jj,month) = diff(squeeze(MLT(:,jj,month)))./DXC(2:end,jj);
        end
        for ii = 1:length(XC)
            dTdy(ii,2:end,month) = diff(squeeze(MLT(ii,:,month)))./DYC(ii,2:end);
        end
    end
    clear ii jj kk MLT
    

    %Put dTdy, dTdx back onto XC, YC, RC grid
    dTdx_at_theta = nan(length(XC),length(YC),length(year_index));
    dTdy_at_theta = nan(length(XC),length(YC),length(year_index));
        
    for ii = 1:12
        dTdx_at_theta(1:end-1,:,ii) = (squeeze(dTdx(1:end-1,:,ii)) + squeeze(dTdx(2:end,:,ii)))./2;
        dTdy_at_theta(:,1:end-1,ii) = (squeeze(dTdy(:,1:end-1,ii)) + squeeze(dTdy(:,2:end,ii)))./2;
    end
    clear dTdx dTdy
        
    heat_advection(:,:,year_index) = -(MLuvel.*dTdx_at_theta + MLvvel.*dTdy_at_theta);
    heat_eddy_advection(:,:,year_index) = -(MLuvel_star.*dTdx_at_theta + MLvvel_star.*dTdy_at_theta);

    clear dTdx dTdy dTdx_at_theta dTdy_at_theta MLvvel MLuvel
end

heat_total_advection = heat_advection + heat_eddy_advection;



%% ENTRAINMENT TERM = dh/dt * deltaT/h

heat_lateral_induction = nan(length(XC),length(YC),length(dd));
heat_entrainment = nan(length(XC),length(YC),length(dd));


%Calculate rate of change of MLD (entrainment velocity)
MLD = ncread('MLD.nc','MLD RF');

dt = diff(dd.*(24*60*60));
dhdt = nan(size(MLD));

for ii = 2:length(dd)-1
    dhdt(:,:,ii) = (MLD(:,:,ii+1) - MLD(:,:,ii-1))./(dt(ii-1) + dt(ii));
end
clear ii dt MLD

dhdt(dhdt<0) = 0; %dh/dt only entrains, no detrainment - negative entrainment velocities set to 0


for year = 1:length(years)
    
    year_index = find(dd>=datenum(years(year),01,01,00,00,00)&dd<datenum(years(year)+1,01,01,00,00,00));

	%Load MLD, MLT
    MLD = ncread('MLD.nc','MLD RF',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    MLT = ncread('ML_ave_properties.nc','MLT',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    
  	%Load potential temperature 
    theta = ncread(strcat('THETA.',num2str(years(year)),'.nc'),'THETA');    
    
    %deltaT = MLT - temperature in first grid cell below MLD
    deltaT = NaN(length(XC),length(YC),length(year_index));
    for ii = 1:length(XC)
        for jj = 1:length(YC)
            for kk = 1:12
                
                if ~isnan(MLD(ii,jj,kk))
                    RF_index = find(RF==MLD(ii,jj,kk));
                    
                    if ~isempty(RF_index)
                        deltaT(ii,jj,kk) = MLT(ii,jj,kk) - theta(ii,jj,RF_index+1,kk);
                    end
                end
            end
        end
    end
    clear ii jj kk RF_index MLT theta
    
    heat_entrainment(:,:,year_index) = -dhdt(:,:,year_index).*(deltaT(:,:,year_index)./MLD);
   
    
        
    %% LATERAL INDUCTION TERM = u*grad(h) * deltaT/h
        
    %Load mixed layer averaged velocities
    MLuvel = ncread('ML_ave_properties.nc','MLuvel',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    MLvvel = ncread('ML_ave_properties.nc','MLvvel',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
   
    
    %Put uvel, vvel onto XC,YC,RC grid - to be same as theta
    %u*grad(H)
    uvel_at_H = nan(length(XC),length(YC),12);
    vvel_at_H = nan(length(XC),length(YC),12);
    
    MLD_at_uvel = nan(length(XC),length(YC),12);
    MLD_at_vvel = nan(length(XC),length(YC),12);
    
    for ii = 1:12
        uvel_at_H(1:end-1,:,ii) = (squeeze(MLuvel(1:end-1,:,ii)) + squeeze(MLuvel(2:end,:,ii)))./2;
        MLD_at_uvel(1:end-1,:,ii) = (squeeze(MLD(1:end-1,:,ii)) + squeeze(MLD(2:end,:,ii)))./2;
        
        vvel_at_H(:,1:end-1,ii)   = (squeeze(MLvvel(:,1:end-1,ii)) + squeeze(MLvvel(:,2:end,ii)))./2;
        MLD_at_vvel(:,1:end-1,ii) = (squeeze(MLD(:,1:end-1,ii)) + squeeze(MLD(:,2:end,ii)))./2;
    end
    
    
    %Calculate horizontal gradient in mixed layer depth
    dHdx = nan(length(XC),length(YC),12);
    dHdy = nan(length(XC),length(YC),12);
    
    for kk = 1:12
        for jj = 1:length(YC)
            dHdx(2:end,jj,kk) = diff(squeeze(MLD_at_uvel(:,jj,kk)))./DXG(2:end,jj);
        end
        for ii = 1:length(XC)
            dHdy(ii,2:end,kk) = diff(squeeze(MLD_at_vvel(ii,:,kk)))./DYG(ii,2:end);
        end
    end
    clear ii jj kk MLD_at_vvel MLD_at_uvel
    
    u_gradh = (uvel_at_H.*dHdx) + (vvel_at_H.*dHdy);
    clear uvel_at_H vvel_at_H dHdx dHdy
    
    heat_lateral_induction(:,:,year_index) = -u_gradh.*(deltaT(:,:,year_index)./MLD);
    clear u_gradh deltaT MLD
end

clear dhdt year_index year years



%% HORIZONTAL DIFFUSION TERM = k*d2T/dx2

kappa = 500;    %Eddy diffusivity - choose value that gives lowest overall resiudal in budget

heat_horizontal_diffusion = nan(length(XC),length(YC),length(dd));

for year = 1:length(years)
    
    year_index = find(dd>=datenum(years(year),01,01,00,00,00)&dd<datenum(years(year)+1,01,01,00,00,00));
        
    %Mixed layer average temperature
    MLT = ncread('ML_ave_properties.nc','MLT',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);    
   
    dTdx = nan(length(XC),length(YC),12);
    dTdy = nan(length(XC),length(YC),12);
    
    for month = 1:12
        for jj = 1:length(YC)
            dTdx(2:end,jj,month) = diff(squeeze(MLT(:,jj,month)))./DXC(2:end,jj);
        end
    end
    clear ii jj kk MLT
    
    
    %Horizontal diffusion term = kappa*u*d2t/dx2 + v*d2t/dy2
    d2Tdx2 = zeros(length(XC),length(YC),12);
    d2Tdy2 = zeros(length(XC),length(YC),12);
        
    for month = 1:12
        
        for jj = 1:length(YC)
            d2Tdx2(2:end,jj,month) = diff(squeeze(dTdx(:,jj,month)))./DXG(1:end-1,jj);
        end
        for xx = 1:length(XC)
            d2Tdy2(ii,2:end,month) = diff(squeeze(dTdy(ii,:,month)))./DYG(ii,1:end-1);
        end
    end
    clear ii jj dTdx dTdy
    
    heat_horizontal_diffusion(:,:,year_index) = kappa.*(d2Tdx2 + d2Tdy2);
    clear d2Tdx2 d2Tdy2
end



%% VERTICAL DIFFUSION TERM = K/h*dT/dz
    
%Vertical diffusivity coefficient
Kz = 1e-4; %Choose coefficient that gives lowest residual in budget

%Heat diffusion term
heat_vertical_diffusion = nan(length(XC),length(YC),length(dd));

for year = 1:length(years)
    
    year_index = find(dd>=datenum(years(year),01,01,00,00,00)&dd<datenum(years(year)+1,01,01,00,00,00));

    %Load MLD, MLT
    MLD = ncread('MLD.nc','MLD RF',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
    MLT = ncread('ML_ave_properties.nc','MLT',[1 1 year_index(1)],[length(XC) length(YC) length(year_index)]);
        
    %Load potential temperature data
    theta = ncread(strcat('THETA.',num2str(years(year)),'.nc'),'THETA');
    
    dT = nan(length(XC),length(YC),length(year_index));
    dz = nan(length(XC),length(YC),length(year_index));

    for ii = 1:length(XC)
        for jj = 1:length(YC)
            for kk = 1:length(year_index)
                
                if ~isnan(MLD(ii,jj,kk))
                    
                    RF_index = find(RF==MLD(ii,jj,kk));
                    RC_index = RF_index;
                    
                    if ~isempty(RC_index)
                        dT(ii,jj,kk) = MLT(ii,jj,kk) - theta(ii,jj,RC_index+1,kk);     %Difference between ave ML theta & theta below

                        %Distance between midpoint of ML and mid-point of cell below
                        dz(ii,jj,kk) = RC(RC_index+1) - (MLD(ii,jj,kk)./2);
                    end
                    clear RF_index RC_index
                end
            end
        end
    end
    clear theta MLT ii jj kk
    
    
    %Diffusion = Kz/H * dT/dz
    heat_vertical_diffusion(:,:,year_index) = -Kz./MLD.*(dT./dz);
    
    clear ii dT dz MLD year_index
end



%% Calculate residual in budget - Sum of budget terms - equal to tendency if no residual
heat_sum = heat_surface_flux + heat_advection + heat_eddy_advection + heat_entrainment + heat_horizontal_diffusion + heat_vertical_diffusion;

%Residual/error in budget
heat_residual = heat_tendency_term - heat_sum;

