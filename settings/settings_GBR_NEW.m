%__________________________________________________________________________
%
% REEFMOD-GBR combined Hindcast (2008-2024.5) and Forecast (2025-2100)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 12/2019. Last update: 02/2025
%__________________________________________________________________________

% load('GBR_REEF_POLYGONS_2026.mat') % should be already loaded in f_multiple_reef

%% INITIALISE REEF STATES WITH AIMS LTMP TRANSECT AND MANTA TOW DATA

%% 1) Determine average community composition
load('AIMS_LTMP_Transect_means.mat') % only used for community composition, not total coral cover
% Complete (1992-2017) LTMP transect dataset at 5-6m averaged per reef (2-3 sites per reef).
% Total cover (TOT) + cover of each coral group (G1....G6).
% Reefs assigned to ReefMod-GBR polygons based on shortest distance between survey and centroid coordinates
% (see 'RE-ARRANGE_LTMP_TRANSECTS.m'). There are discrepancies in reef names and shelf position 
% between AIMS database and Karlo's list of 3806. using Karlo's designations here for consistency
% 'YEAR_Reefmod' is modelled year (surveys assigned to each season from their actual dates) 

survey_year = LTMP_Transects_means.YEAR_Reefmod;
select_LTMP_compo = find(survey_year>=2004 & survey_year<2009); % select transect data from 2004 to 2008 to initialise community compo

TOT_COVER = table2array(LTMP_Transects_means(select_LTMP_compo,9)); % 341 reefs surveyed between 2004-08 (123 unique reefs)
Z = table2array(LTMP_Transects_means(select_LTMP_compo,10:15))./TOT_COVER(:,ones(6,1)); % relative proportion of each group (each row sums to 1)

% Extract AIMS sector (2), shelf (3), reef name (4), ReefMod year (5) and ReefMod reef_ID (18) for the selected years
LTMP_Transects_compos = [LTMP_Transects_means(select_LTMP_compo,[16 2 3 4 5]) array2table(Z)]; 
% Rename for consistency across tables
LTMP_Transects_compos.Properties.VariableNames{1} = 'Reef_ID';
LTMP_Transects_compos.Properties.VariableNames{5} = 'YearReefMod';
% Order by increasing Reef_ID
[~,new_order] = sort(LTMP_Transects_compos.Reef_ID);
LTMP_Transects_compos = LTMP_Transects_compos(new_order,:);
% Assign shelf position as designated in ReefMod
NEW_SHELF_DESIGNATION = table2array(innerjoin(LTMP_Transects_compos,GBR_REEFS,'LeftKeys', 'Reef_ID','RightKeys', 'Reef_ID', 'LeftVariables',[1],'RightVariables',[ 10 11]));
LTMP_Transects_compos.Shelf_position = NEW_SHELF_DESIGNATION(:,2);
LTMP_Transects_compos.AIMS_sector = NEW_SHELF_DESIGNATION(:,3);

% Calculate average community composition (relative props) for each sector(1-11)*shelf positon(1-3)
omean = @(x) mean(x,'omitnan');
AIMS_AVERAGECOMPO = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables',{'AIMS_sector','Shelf_position'},...
    'InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'})); % 20 combinations of sector*shelf (out of 33)

% Add the missing sector/shelf combinations. For missing compo take the one of the preceding sector (north) for the same shelf position
SECTOR = AIMS_AVERAGECOMPO(:,1);
SHELF = AIMS_AVERAGECOMPO(:,2);

% Note there is no inshore reef in the Swain sector (sector 10)!
Missings = [1 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;...
            1 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            1 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            2 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;...
            2 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            2 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            3 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;... 
            5 3 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==3,4:9) ;...          
            7 1 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==1,4:9) ;...
            7 2 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==2,4:9) ;...
            7 3 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==3,4:9) ;...
            9 3 NaN AIMS_AVERAGECOMPO(SECTOR==8&SHELF==3,4:9) ];
            
AIMS_AVERAGECOMPO = array2table([AIMS_AVERAGECOMPO ; Missings]) ; % Now 32 -> all possible sector*shelf combinations populated

% Assign to each of the 3806 reefs the average community compo representative of shelf and sector
relative_cover = table2array(outerjoin(GBR_REEFS,AIMS_AVERAGECOMPO, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4:9]));

[~,I]=sort(relative_cover(:,1));
relative_cover = relative_cover(I,:); % sort by increasing reef ID

% Then swap for the surveyed reefs (n=123) which compo was actually recorded (averaged over the selected period (2006-2008)
AIMS_AVERAGECOMPO_REEF = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables','Reef_ID','InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'}));

relative_cover(AIMS_AVERAGECOMPO_REEF(:,1),4:end) = AIMS_AVERAGECOMPO_REEF(:,3:end); % including 123 reefs with actual community composition recorded

%% 2) Assign to each of the 3806 reefs the average total cover representative of shelf and sector between 2006-2008
% Load the  linear model allowing to convert manta tow into transect equivalent estimates
load('LTMP_Transect2Tow_2024.mat')

% Load last updated AIMS observations including manta and fixed transect (LTMP+MMP) [generated after data extraction from RimRep-DMS]
load('GBR_AIMS_OBS_CORAL_COVER_2024.mat')
all_project_codes = [1:3]; % 1 for 

DATA_TR_LTMP = GBR_AIMS_OBS_CORAL_COVER(strcmp(string(GBR_AIMS_OBS_CORAL_COVER.project_code),'LTMP')==1 & ...
    strcmp(string(GBR_AIMS_OBS_CORAL_COVER.data_type),'photo-transect')==1,:);
DATA_TR_LTMP.CCOVER = 100*DATA_TR_LTMP.mean; % calculate percent coral cover

DATA_TR_MMP = GBR_AIMS_OBS_CORAL_COVER(strcmp(string(GBR_AIMS_OBS_CORAL_COVER.project_code),'MMP')==1 & ...
    strcmp(string(GBR_AIMS_OBS_CORAL_COVER.data_type),'photo-transect')==1 & ...
    GBR_AIMS_OBS_CORAL_COVER.depth>=5,:); % exclude shallow observations (2m depth)
DATA_TR_MMP.CCOVER = 100*DATA_TR_MMP.mean; % calculate percent coral cover

DATA_MT_LTMP = GBR_AIMS_OBS_CORAL_COVER(strcmp(string(GBR_AIMS_OBS_CORAL_COVER.project_code),'LTMP')==1 & ...
    strcmp(string(GBR_AIMS_OBS_CORAL_COVER.data_type),'manta')==1,:);
DATA_MT_LTMP.CCOVER = predict(LTMP_Transect2Tow_Model2, 100*DATA_MT_LTMP.mean); % calculate percent coral cover (transect equivalent

% Group all data
DATA_ALL = [DATA_TR_LTMP ; DATA_TR_MMP ; DATA_MT_LTMP];

StartDate = 2006 ; % earliest date used to initialize coral cover
EndDate = 2009 ; % exclusive

select_DATA = find(DATA_ALL.YearReefMod >= StartDate & DATA_ALL.YearReefMod < EndDate);

% Calculate the mean total coral cover per reef across all surveys
DATA_MEAN = varfun(omean, DATA_ALL(select_DATA,:),'GroupingVariables',{'Reef_ID';'Shelf_position';'AIMS_sector'}, 'InputVariables','CCOVER'); %average for each reef

% Average per shelf position and sector 
CORAL_AVERAGESTATE = table2array(varfun(omean, DATA_MEAN,'GroupingVariables', {'AIMS_sector','Shelf_position'}, 'InputVariables',{'Fun_CCOVER'})); 
% Exclude sectors represented by less than 3 reefs
CORAL_AVERAGESTATE = CORAL_AVERAGESTATE(CORAL_AVERAGESTATE(:,3)>=3,:); % 19 combinations of sector*shelf (out of 32)

% Add the missing sector/shelf combinations
SECTOR = CORAL_AVERAGESTATE(:,1);
SHELF = CORAL_AVERAGESTATE(:,2);

% Note there is no inshore reef in the Swain sector!
Missings = [1 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...
    1 2 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    1 3 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    2 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...
    2 2 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    2 3 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    3 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...  
    5 3 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==3,4) ;...
    7 1 NaN CORAL_AVERAGESTATE(SECTOR==6&SHELF==1,4) ;...
    7 3 NaN CORAL_AVERAGESTATE(SECTOR==6&SHELF==3,4) ;...
    9 1 NaN CORAL_AVERAGESTATE(SECTOR==8&SHELF==1,4) ;...
    9 3 NaN CORAL_AVERAGESTATE(SECTOR==8&SHELF==3,4) ;...
    11 2 NaN CORAL_AVERAGESTATE(SECTOR==9&SHELF==2,4) ] ;
        
CORAL_AVERAGESTATE = array2table([CORAL_AVERAGESTATE ; Missings],'VariableNames', {'AIMS_sector';'Shelf_position';'N';'TOT'} ) ;

init_coral = table2array(outerjoin(GBR_REEFS,CORAL_AVERAGESTATE, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4]));

[~,I]=sort(init_coral(:,1));
init_coral = init_coral(I,:); % sort by increasing reef ID (just in case)

% Then for each surveyed reef swap the assigned average by the actual value of total cover
% This gives x reefs initialised with the total cover observed, and y reefs initialised with a regional average
init_coral(DATA_MEAN.Reef_ID,4) = DATA_MEAN.Fun_CCOVER;

clear AIMS_AVERAGECOMPO AIMS_AVERAGECOMPO_REEF DATA_ALL DATA_MEAN DATA_MT_LTMP DATA_TR_LTMP DATA_TR_MMP
clear GBR_AIMS_OBS_CORAL_COVER LTMP_Transects_means LTMP_Transects_compos LTMP_Transect2Tow_Model2
clear CORAL_AVERAGESTATE select_DATA select_LTMP_compo NEW_SHELF_DESIGNATION

%% 3) Now populate each reef with a random cover (normally distributed around specified mean)
X = init_coral(META.reef_ID,4); % Percent cover of all corals
Y = relative_cover(META.reef_ID,4:end);
varSDcoral = 0.2*ones(1,META.nb_coral_types);

init_coral_cover = X(:,ones(size(Y,2),1)).* Y ;

for s=1:6
    init_coral_cover(:,s) = normrnd(init_coral_cover(:,s), varSDcoral(1,s)*init_coral_cover(:,s))/100;
end

init_coral_cover(init_coral_cover<0.005)=0.005; % minimum 0.5% cover for each group
TCmax = 0.8; %maximum initial total cover is 80%
select_TCmax = find(sum(init_coral_cover,2)>TCmax);
adjust_TCmax = sum(init_coral_cover(select_TCmax,:),2);
init_coral_cover(select_TCmax,:)=TCmax*init_coral_cover(select_TCmax,:)./adjust_TCmax(:,ones(size(init_coral_cover,2),1));

%% Now populate random rubble and sand covers
init_rubble = 0.1; % 10% on average 
init_sand = 0.3; % 30% sand on average

varSDother = 0.2;

init_rubble_cover = normrnd(init_rubble*100, varSDother*init_rubble*100, length(META.reef_ID), 1)/100 ;
init_rubble_cover(init_rubble_cover<0.01) = 0.01;

% Update Nov 2022: now sand determined by Roelfsema et al (2021)'s geomorphic and benthic maps
init_sand_cover = normrnd(GBR_REEFS.UNGRAZABLE(META.reef_ID)*100, varSDother*GBR_REEFS.UNGRAZABLE(META.reef_ID)*100)/100 ;
init_sand_cover(init_sand_cover<0.05) = 0.05; % impose minimum of 5%

% Adjust if too high (coral + sand > 95%). We leave 5% free space for a safe intialisation
CHECK = sum(init_coral_cover,2) + init_sand_cover;
init_sand_cover(CHECK>0.95) = 0.95 - sum(init_coral_cover(CHECK>0.95,:),2);

%% Populate initial algal cover
init_algal_cover = 0.001*ones(META.nb_reefs,META.nb_coral_types);

clear init_coral init_sand init_rubble relative_cover X Y varSD adjust_TCmax select_TCmax TCmax varSDother CHECK

%% COTS densities - INPUT DENSITIES MUST BE PER GRID (400m2), NOT PER TOW
REEF_COTS.densities_M = nan(META.nb_reefs,META.nb_time_steps+1);
REEF_COTS.densities_SD = nan(META.nb_reefs,META.nb_time_steps+1);

% First, initialize CoTS populations with CoCoNet hindcast densities of adult CoTS (mean and SD per grid) provided by Scott
load('GBR_past_COTS_CoCoNet.mat')
% Matrices of mean density and SD already converted into nb CoTS per grid for the period 1985-2017.
% Values are mean based on 50 replicate runs, interpolated from 3638 reefs - see 'REEF_COTS_INTERPOLATION_1985_2017.m'
% So 33 years * 2 seasons = 66 columns
CoCoNet_years = 1985:0.5:2017.5;

start_year = find(CoCoNet_years==2007.5) ; %first column is 1985, so 2007 is 1+22
REEF_COTS.densities_M(:,1) = PAST_COTS_M_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid
REEF_COTS.densities_SD(:,1) = PAST_COTS_SD_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid

% OCT 2025: forcing of CoTS hindcast density with past observations now disabled thanks to new calibration of CoTS mortality
% and larval-stock recruitment relationships (Suki's optimisation against CoTS observations with GBRLUP connectivity
% We only use available observations at initial step (replacing CoCONet predictions where observations available between 2006-2007.5
% Mean coTS per tow from AIMS LTMP (ID=1), FMP (ID=2) and Control Program (ID=3). Missing values identified as NaN.
load('GBR_PAST_COTS_1992_2025.mat') % 67 columns from 1992 to 2022.5 inclusive
% 2025, 2024.5, 2024, 2023.5: only CP ; 2023: CP+FMP ; 2022.5: CP+FMP+LTMP

GBR_PAST_COTS_PER_GRID = (GBR_PAST_COTS_NEW(META.reef_ID,:)/0.22)*(1500/2500); % Conversion: 1500 COTS per km2 ~0.22 per manta tow (Moran and De'ath 1986)
% which is further divided by 2500 to get number of COTS per 400m2 (reef grid size)

COTS_years = 1992:0.5:2025;
start_step = find(COTS_years==2007.5);
end_step = size(GBR_PAST_COTS_PER_GRID,2);
length_history = length(start_step:end_step);

RECENT_COTS = max(GBR_PAST_COTS_PER_GRID(:,(start_step-3):start_step),[],2); % maximum density observed from 2006 to 2007.5

% Replace CoCoNet predictions with recent observations (2006-2007.5)
REEF_COTS.densities_M(find(isnan(RECENT_COTS)==0),1)= RECENT_COTS(find(isnan(RECENT_COTS)==0));
REEF_COTS.densities_SD(find(isnan(RECENT_COTS)==0),1)= 0.28 * RECENT_COTS(find(isnan(RECENT_COTS)==0)); 
% (with 0.28 = slope of linear model between SD and mean CoCoNet predictions)

% Populate from 2008 onwards with observations
% (OCT 2025: needs to be commented if using GBRLUP CoTS connectivity with new params BH.beta and BH.alpha in setting_COTS)
REEF_COTS.densities_M(:,2:length_history) = GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;
REEF_COTS.densities_SD(:,2:length_history) = 0.28 * GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;

clear GBR_PAST_COTS about_GBR_COTS PAST_COTS_M_COCONET PAST_COTS_SD_COCONET CoCoNet_years end_step
clear GBR_PAST_COTS_NEW GBR_PAST_COTS_NEW_ID GBR_PAST_COTS_PER_GRID about_GBR_PAST_COTS RECENT_COTS

%% SCENARIO OF THERMAL STRESS
if META.doing_genetics == 0 % NO GENETIC ADAPTATION (WITH THERMAL OPTIMUM)
    
    DHW = zeros(length(META.reef_ID),META.nb_time_steps);
    
    %% PAST THERMAL STRESS REGIME
    % Using NOAA Coral Reef Watch 5km product: max DHW every year from 1985 to 2023 from closest 5x5 km pixel (updated Jan 2024)   
    load('GBR_past_DHW_CRW_5km_1985_2024.mat')
    GBR_PAST_DHW = GBR_PAST_DHW(:,24:end); % col 24 is for 2008
  
    end_hindcast = size(GBR_PAST_DHW,2);
    DHW(:,1:2:end_hindcast*2) = GBR_PAST_DHW(META.reef_ID,:);
    clear GBR_PAST_DHW

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


%% PAST STORM REGIME 2008 to 2024
% Category cyclones (Saffir-Simpson scale) assigned to the spatial prediction of >4m wave height from Puotinen et al (2016).
% Original set of predictions is matrice of 0/1 compiled by Rob for the entire GBR (3806 reefs). 
% Then, occurence of damaging wave is blended with category cyclone estimated from maximum sustained wind and distance to cyclone track
% (with 10% increase to meet US standards) measured along each real cyclone track (data extracted from BoM Database of past cyclone tracks).
% Feb 2025: now modelling the wind field to estimate past cyclone categories
% load('GBR_cyclones_2008-2024.mat') % Using Holland's model
% load('GBR_cyclones_2008-2024_BOOSE.mat')
load('GBR_cyclones_2008-2024_HOURLY.mat')% Using Holland's model with asymmetry and cyclone tracks with hourly interpolation

CYCLONE_CAT = zeros(length(META.reef_ID),2*size(GBR_PAST_CYCLONES,2)) ;
CYCLONE_CAT(:,1:2:end) = GBR_PAST_CYCLONES(META.reef_ID,:);
clear GBR_PAST_CYCLONES

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

%% Update March 2024: add probability of incidence of cyclone striking a reef (based on orbital velocity thresholds)
% This only works if META.randomize_hurricane_strike set to 1 in f_multiple_reef
load('GBR_coral_breakage_probability.mat') % BREAKAGE_PROBABILITY = probability of exceeding threshold of breakage
BREAKAGE_PROBABILITY = BREAKAGE_PROBABILITY(META.reef_ID,:);
STRIKE_INCIDENCE = BREAKAGE_PROBABILITY.AvgPr_arbo_lge_041; % select the risk of cyclone strike based on available thresholds
% STRIKE_INCIDENCE = BREAKAGE_PROBABILITY.AvgPr_acro_lge05_18; % select the risk of cyclone strike based on available thresholds
STRIKE_INCIDENCE(isnan(STRIKE_INCIDENCE)==1)=1; % if not defined (NaN), assume proba=1

clear GBR_DHM_Max DHM_Max FUTURE_COTS_M FUTURE_COTS_SD CURRENT_CORAL_COVER CURRENT_COTS_DENSITIES
clear CURRENT_RUBBLE_COVER CURRENT_NONGRAZABLE ALL_CYCLONE_CAT Cyc_cat DHWtmp Model_reef_indices ModelReef_ID ModelReef_ID_label
clear SST_GBRtmp Years_unique yr nb_hindcast_runs rel T_opt
clear nb_time_steps Reef_SST_Mean_Yr total_coral_cover Yr Topt_scenario_YM varSD resize_Topt_list
clear FUTURE_CYCLONE_CAT AIMS_ALL AIMS_ALLtmp AIMS_MT AIMS_TR max_annual_DHW LTMP_Transect2Tow_Model
