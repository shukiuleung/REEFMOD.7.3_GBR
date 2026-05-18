%__________________________________________________________________________
%
% REEFMOD-GBR combined Hindcast (2008-2024.5) and Forecast (2025-2100)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 12/2019. Last update: 02/2025
% Revised: 02/2026
%__________________________________________________________________________

%% Populate initial coral cover for each reef (Feb 2026)
load('CORAL_INIT_2007.mat') % Feb 2026: contains the average percent cover of each taxonomic group for intialisation.
% Derived from monitoring data (transects + manta tow between 2006-2008.5), see GENERATE_INITIAL_CORAL_COVER.
init_coral_cover = init_coral_cover(META.reef_ID,:); % only keep values for the simulated reefs

varSDcoral = 0.2*ones(1,META.nb_coral_types); % assumed variance for each taxonomic group to generate random noise

% Create a randomised matrix of proportional cover (between 0-1 for each group, sum per reef < 1)
for s=1:6
    init_coral_cover(:,s) = normrnd(init_coral_cover(:,s), varSDcoral(1,s)*init_coral_cover(:,s))/100;
end

init_coral_cover(init_coral_cover<0.005)=0.005; % minimum 0.5% cover for each group
TCmax = 0.8; %maximum initial total cover is 80%
select_TCmax = find(sum(init_coral_cover,2)>TCmax);
adjust_TCmax = sum(init_coral_cover(select_TCmax,:),2);
init_coral_cover(select_TCmax,:) = TCmax*init_coral_cover(select_TCmax,:)./adjust_TCmax(:,ones(size(init_coral_cover,2),1));

%% Populate random rubble and sand covers
varSDother = 0.2;

init_rubble = 0.1; % 10% rubble on average for every reef

init_rubble_cover = normrnd(init_rubble*100, varSDother*init_rubble*100, length(META.reef_ID), 1)/100 ;
init_rubble_cover(init_rubble_cover<0.01) = 0.01;

% Update Nov 2022: now sand determined by Roelfsema et al (2021)'s geomorphic and benthic maps
init_sand_cover = normrnd(GBR_REEFS.UNGRAZABLE(META.reef_ID)*100, varSDother*GBR_REEFS.UNGRAZABLE(META.reef_ID)*100)/100 ;
init_sand_cover(init_sand_cover<0.05) = 0.05; % impose minimum of 5%

% Adjust if too high (coral + sand > 95%). We leave 5% free space for a safe intialisation
CHECK = sum(init_coral_cover,2) + init_sand_cover;
init_sand_cover(CHECK>0.95) = 0.95 - sum(init_coral_cover(CHECK>0.95,:),2);

%% Populate initial algal cover
init_algal_cover = 0.001*ones(META.nb_reefs,META.nb_algal_types);

clear init_rubble relative_cover varSDcoral adjust_TCmax select_TCmax TCmax varSDother CHECK

%% COTS densities - Feb 2026: now input as CoTS per tow (NOT density)
% Derived from monitoring data (manta tow between 2006-2008.5), see GENERATE_INITIAL_COTS_PER_TOW.
% For each monitored reef we keep the original observation (maximum over the selected period).
% For non-monitored reefs (0% Far North, 20% Cairns/Cooktown, 5% Townsville/Whitsunday, 2% Mackay/Capricorn -> 155 reefs
% in total, 4% of the GBR) we generated 1/0 based on the frequency distribution of non-zero/zero observations in each region,
% then we predicted density based on gamma distribution fitted to the non-zero observations.
load('COTS_INIT_2007.mat')
% Contains 'CoTS_per_tow' which gives predicted+observed CoTS per tow value for every of the 3806 reefs in 2007
% and 'observed' which identifies whether each CoTS per tow value comes from a prediction (0) or observation (1).

% For every new run, we shuffle the predicted values within each region (and always keep the original observations)
ListRegions = unique(GBR_REEFS.AREA_DESCR);
for this_region = 1:length(ListRegions)
    ID_region = find(strcmp(GBR_REEFS.AREA_DESCR, ListRegions(this_region)) & COTS_INIT.observed ==0); % flag which values were predicted
    COTS_INIT.CoTS_per_tow(ID_region) = COTS_INIT.CoTS_per_tow(ID_region(randperm(length(ID_region)))); % shuffle these values
end

init_COTS_per_tow = COTS_INIT.CoTS_per_tow(META.reef_ID); % define initial coTS conditions for the simulated reefs

% Feb 2026: no need to force the CoTS hindcast with past observations, thanks to Suki's optimised calibration of CoTS mortality
% + larval-stock recruitment. Only need to initialise in 2007.

clear COTS_INIT ID_region ListRegions this_region

%% SCENARIO OF THERMAL STRESS
if META.doing_genetics == 0 % NO GENETIC ADAPTATION (WITH THERMAL OPTIMUM)
    
    DHW = zeros(length(META.reef_ID),META.nb_time_steps);
    
    %% PAST THERMAL STRESS REGIME
    % Using NOAA Coral Reef Watch 5km product: max DHW every year from 1985 to 2023 from closest 5x5 km pixel (updated Jan 2024)   
    load('GBR_past_DHW_CRW_5km_1985_2026.mat')
    GBR_PAST_DHW = GBR_PAST_DHW(:,24:end); % col 24 is for 2008
  
    end_hindcast = size(GBR_PAST_DHW,2);
    DHW(:,1:2:end_hindcast*2) = GBR_PAST_DHW(META.reef_ID,:);

    %% FUTURE THERMAL STRESS (CMIP6)
    start_future = 35; % Forecast starts in 2025 (step 35)
    % Need to be odd number (always start future in summer)
    
    if META.nb_time_steps >= start_future
        
        % The DHW matrices start in 2000 (column 6) so year 2025 is column 31
        
        if OPTIONS.SSP=='119' && ismember(OPTIONS.GCM, {'GFDL-ESM4', 'MPI-ESM1-2-HR', 'NorESM2-LM'})==1           
            error('##### REEFMOD ERROR #### SSP1-1.9 is not available for the climate models FDL-ESM4, MPI-ESM1-2-HR and NorESM2-LM')     
        end
        
        % Load the selected forecast scenario of DHW
        load([char(OPTIONS.GCM) '_' char(OPTIONS.SSP) '_annual_DHW_max.mat'])
       
        if META.nb_time_steps-start_future-1 > 2*size(max_annual_DHW(:,31:end),2)   
            error('##### REEFMOD ERROR #### Simulated timeframe inconsistent with the DHW forecast. NB_TIME_STEPS has to be <= 31+155')
        else
            
            % Select the specified timeframe in the available forecast
            DHW_FORECAST = max_annual_DHW(META.reef_ID,31:(31+(META.nb_time_steps-start_future-1)/2));
            % DHW_FORECAST_shuffled = DHW_FORECAST; % deactivate the shuffling
            % Let's shuffle available years within each decade
            DHW_FORECAST_shuffled = nan(size(DHW_FORECAST));
            start = 1; % first year of the selected forecast
            remain = size(DHW_FORECAST,2); % number of years still available

            while remain > 10

                sample = randperm(10); % sample at random within the decade
                DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample); % assign the shuffled years
                start = start + length(sample);
                remain = remain - length(sample);
            end

            % Last years available to be shuffled as well (Less than a decade s available)
            sample = randperm(remain);
            DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample);
            
            % Finally assign to the DHW matrix (only in summers)
            DHW(:,start_future:2:end) = DHW_FORECAST_shuffled;
        end
    end
    
else
    error('Cannot run genetic adaptation with this version of ReefMod')
end

clear GBR_PAST_DHW end_hindcast max_annual_DHW DHW_FORECAST_shuffled DHW_FORECAST sample remain start start_future 

%% PAST STORM REGIME 2008 to 2024
% Category cyclones (Saffir-Simpson scale) assigned to the spatial prediction of >4m wave height from Puotinen et al (2016).
% Original set of predictions is matrice of 0/1 compiled by Rob for the entire GBR (3806 reefs). 
% Then, occurence of damaging wave is blended with category cyclone estimated from maximum sustained wind and distance to cyclone track
% (with 10% increase to meet US standards) measured along each real cyclone track (data extracted from BoM Database of past cyclone tracks).
% Feb 2025: now modelling the wind field to estimate past cyclone categories
% load('GBR_cyclones_2008-2024.mat') % Using Holland's model
% load('GBR_cyclones_2008-2024_BOOSE.mat')
load('GBR_cyclones_2008-2026.mat')% Using Holland's model with asymmetry and cyclone tracks with hourly interpolation

CYCLONE_CAT = zeros(length(META.reef_ID),2*size(GBR_PAST_CYCLONES,2)) ;
CYCLONE_CAT(:,1:2:end) = GBR_PAST_CYCLONES(META.reef_ID,:);

% Set the mitigating of cyclone on bleaching to none for the hindcast
% Cyclones and bleaching occurred during the same season in 2017 and 2024. 
% In 2017, ~997 reefs were exposed to both, with Debbie occurring in late March
% Hughes et al. 2019: severe tropical cyclone Debbie crossed the southern Great Barrier Reef at approximately 20° S on 27–28 March 2017.
% However, the resulting wind, cloud and rain was 4–6 weeks too late and too far south to moderate the second bout of severe bleaching.
% In 2024, cyclone Jasper occurredn in dec and Kirrily in Jan, before the mass bleaching.
META.allow_cyclone_cooling = ones(1, META.nb_time_steps);
META.allow_cyclone_cooling(1:17) = 0;
    
%% FUTURE STORM REGIME (from 2025 onwards)
if META.nb_time_steps > 34
    
    load('Reef_Cyclone_TimeSeries_Count_Cat.mat')
    
    FUTURE_CYCLONE_CAT =  zeros(length(META.reef_ID),2*size(Cyc_cat,2)) ;
    FUTURE_CYCLONE_CAT(:,1:2:end) = Cyc_cat(GBR_REEFS.Nick_ID(META.reef_ID),:,simul);
    
    CYCLONE_CAT = [CYCLONE_CAT FUTURE_CYCLONE_CAT];
    
end

clear GBR_PAST_CYCLONES FUTURE_CYCLONE_CAT Cyc_cat

%% Update March 2024: add probability of incidence of cyclone striking a reef (based on orbital velocity thresholds)
% This only works if META.randomize_hurricane_strike set to 1 in f_multiple_reef
load('GBR_coral_breakage_probability.mat') % BREAKAGE_PROBABILITY = probability of exceeding threshold of breakage
BREAKAGE_PROBABILITY = BREAKAGE_PROBABILITY(META.reef_ID,:);
STRIKE_INCIDENCE = BREAKAGE_PROBABILITY.AvgPr_arbo_lge_041; % select the risk of cyclone strike based on available thresholds
% STRIKE_INCIDENCE = BREAKAGE_PROBABILITY.AvgPr_acro_lge05_18; % select the risk of cyclone strike based on available thresholds
STRIKE_INCIDENCE(isnan(STRIKE_INCIDENCE)==1)=1; % if not defined (NaN), assume proba=1

clear yr rel nb_time_steps