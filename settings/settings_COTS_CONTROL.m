%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tina Skinner, MSEL, June 2023.
%
% Last edited: Suki Leung, MSEL, March 2026
% Parametrisation of COTS control settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

META.COTS_reef_trigger = 0; % whether (1) or not (0) the trigger of control is reef-level CoTS above threshold of 0.22 CoTS per tow
% If set to 0, the trigger is CoTS density per site, not per reef. The alternative (1) hasn't been tested yet.

META.COTS_control_start = 23; % timestep in 6-month intervals when control should start = 23 (summer 2019), 
%%% CHOOSING REEFS TO CONTROL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YM 07/25: now defining the region where CoTS control is deployed (region names to match 'New_regions_TS.mat').
% With 4 regions, a total of 15 geographic domains can be defined, including the previous stategies #1 (GBR-wide) to #8 (only South). 
% Note this also enables other stategies (#9 to #18) to perform in a specific geographic context (i.e., GBR-wide vs other specified domains)
% If only a subset of reefs if selected for simulation, make sure to specify the corresponding region, otherwise no reef might be visited. 
META.COTS_cull_region = {'FN';'N';'C';'S'};  % All regions -> GBR-wide (previously strategy 1)
% META.COTS_cull_region = {'FN'};  % Far North (previously strategy #2)
% META.COTS_cull_region = {'N'};  % North (previously strategy #4)
% META.COTS_cull_region = {'C'};  % Centre (previously strategy #6)
% META.COTS_cull_region = {'S'};  % South (previously strategy #8)
% META.COTS_cull_region = {'N';'C';'S'}; 
% Suki May 2026: note updated New_regions_TS_SL includes 'out'.

META.doing_COTS_heat_stress_detectibility = 1; % Uses evidance from Cook et al. in review

if META.doing_COTS_heat_stress_detectibility == 1 
    load("log_culling_dhw_effects.mat"); % Results from Cook et al. 2025 - How does culling change under DHW per size class (roughly mapped from the 4 size classes to the 16 here)
    META.log_culling_dhw_effects = log_culling_dhw_effects; % Note: will only work with 16 size classes
end

% ========================================================================
% Suki Feb 2026
% Create pair-wise reef distance matrix 
reef_coords = [GBR_REEFS.LON, GBR_REEFS.LAT];
% Preallocate distance matrix
META.distance_matrix = zeros(size(reef_coords, 1));

% Compute distance between all reef pairs using lldistkm (Haversine formula)
for i = 1:size(reef_coords, 1)
    for j = i+1:size(reef_coords, 1)
        d1km = lldistkm([reef_coords(i,2), reef_coords(i,1)], ...
                        [reef_coords(j,2), reef_coords(j,1)]);
        META.distance_matrix(i,j) = d1km;
        META.distance_matrix(j,i) = META.distance_matrix(i,j); % Make symmetric
    end
end

% ========================================================================

% Choose culling strategy (detail in f_makeReefList_NEW)
META.COTS_reefs2cull_strat = 1;
% 1 - GBRMPA strategy that goes to Target reefs first, then Priority reefs, then Non Priority reefs
% 9 - Outbreak front (latitude): GBRMPA strategy that goes to Target reefs first, then also goes to 0.5° lat (~50 km) from target reefs with outbreaks - whole GBR.
% 10 - Outbreak front (sector): look for the AIMS sector (1-11) that has the highest density of COTS on ALL reefs, start control there, then remaining.
% 11 - Outbreak front (sector/priority): look for the AIMS sector (1-11) that has the highest proportion of PRIORITY reefs with outbreaks, start control there, then remaining.
% 12 - Connectivity-based: effort sink - at each timestep, still make priority reef list same as before but don't include reefs that have COTS > 3 per tow at that timestep.
% 13 - Connectivity/coral-based: effort sink with coral cover minimum - at each timestep, still make priority reef list same as before but don't include reefs that have COTS > 3 per tow at that timestep unless they have > threshold coral cover.
% 14 - Protection status: at each timestep, still make priority reef list etc but only include green zone reefs.
% 15 - Protection status: at each timestep, still make priority reef list etc but don't include any green zone reefs.
% 16 - GreenZone weighting with CoTS connectivity and coral cover.
% 17 - BlueZone weighting with CoTS connectivity and coral cover.
% 18 - CoTS connec and coral cover, same as case 16 and 17 above, but no preferential weighting to blue or green zones.
% 19 - Control maximise benefits (WIP)
% 20 - Adaptive Control: default (cat 1 cull, cat 2 +monitor, cat 3 +culling, cat 4 relocate 250km, cat 5 relocate 1000km)
% 21 - Adaptive Control shifted +1 cat (cat 1-2 cull, cat 3 +monitor, cat 4 +culling, cat 5 relocate 1000km)
% 22 - Adaptive Control shifted -1 cat (cat 1 +monitor, cat 2 +culling, cat 3-4 relocate 250km, cat 5 relocate 1000km)
% 23 - Adaptive Control with alternative regional bleaching threshold: skip region when >=60% of reefs are severely bleached
% 24 - Adaptive Control with alternative regional bleaching threshold: skip region when >=30% of reefs are severely bleached
% 25 - Adaptive Control with alternative regional bleaching decisions: expand effort to NP reefs when at least 3 regions are skip (>=3 unworkable -> no redistribution).
% 26 - Adaptive Control with alternative regional bleaching decisions: always expand to NP reefs (never redistribute, regardless of how many regions are skip).
% 27 - Adaptive Control with alternative knowledge level: perfect knowledge a - score based on CoTS risk (density) + predicted manta tow density
% 28 - Adaptive Control with alternative knowledge level: perfect knowledge b - score based on CoTS benefits (culled density) + predicted manta tow density
% 29 - Adaptive Control and maximise benefits  (WIP)


% Suki April 2026 WIP - dynamic strategy schedule. This allows switching between strategies at different time points in the simulation, which could be used to simulate adaptive management (e.g. start with strategy 1, then switch to strategy 19 when bleaching starts).
% Strategy schedule: Nx2 matrix [t_start, strategy], rows in ascending order of t_start.
% Each row activates 'strategy' from timestep t_start onwards, overriding META.COTS_reefs2cull_strat.
% Leave empty ([]) to use META.COTS_reefs2cull_strat for all timesteps.
% Example: run strat 1 from control start, then switch to strat 19 from t=35 onwards:
% META.COTS_strat_schedule = [META.COTS_control_start, 1; 35, 19];
META.COTS_strat_schedule = [META.COTS_control_start, 1; 37, 20];

% Updates from Tina (July 2023) - now have fixed target reef list. Control at T, then P, then N.
load('New_regions_TS_SL.mat') %this has been updated for new GBRMPA 2023 PR list and includes target reefs now
% Suki May 05: updated again to match Isobel's definition of in/out of
% GBRMP
% Tina Jul 25: no out (outside) anymore, all reefs in the Top North now included in the MP
META.COTS_cull_reeflist = nregions(META.reef_ID,:);

% Number of cull sites per reef: Tina April 2023 new based on Geom_CH_2D_km2. For each reef, 2D coral habitat area in first column, no. cull sites in second. 
load('COTS_sites_new.mat'); %Note, all cull sites integers or model crashes. All 0's rounded up to 1 as number of sites per reef should not be 0, given it's calculated on the area of 2D coral habitat.
META.COTS_cull_reeflist.nb_sites = COTS_sites(META.reef_ID,2); % only retain the number of sites for each reef included in simulation
META.COTS_cull_reeflist.AIMS_sector = GBR_REEFS.AIMS_sector(META.reef_ID); % add AIMS sector
META.COTS_cull_reeflist.GreenZone = GBR_REEFS.GreenZone(META.reef_ID); % add AIMS sector

%% Build the starting priority list for this run, adjusted to the specified geographic domain
is_focus = ismember(META.COTS_cull_reeflist.Region, META.COTS_cull_region); % Find all reefs in the focused region(s)
focus = META.COTS_cull_reeflist(is_focus, :); % Only keep reefs from selected region

is_target = ismember(focus.reef_type, 'T'); % Find all Target Reefs (T)
target = focus(is_target, :);   %Select only TR from selected region

is_priority = ismember(focus.reef_type, 'P'); %Find all subsequent Priority Reefs (P)
priority = focus(is_priority, :);   %Select only PR from selected region

is_nonpriority = ismember(focus.reef_type, 'N'); %Find all Non Prority Reefs (N)
nonpriority = focus(is_nonpriority, :);     %Select all N from selected region

% Now shuffle each list at least once, such that each run starts with a different prioritisation
% META.COTS_fixed_list will tell whether this list is shuffled again at every time step
META.COTS_cull_reeflist_targetRUN = target(randperm(size(target,1)),:); %Randomly shuffle the target list
META.COTS_cull_reeflist_priorityRUN = priority(randperm(size(priority,1)),:); %Randomly shuffle the PR list
META.COTS_cull_reeflist_nonpriorityRUN = nonpriority(randperm(size(nonpriority,1)),:); %Randomly shuffle the NPR list

%% Thermal refugia/hotspot status for each reef and GCM which may determine control priority - Tina Jan 2026
load('refugia_lookup.mat');

% Filter lookup table for this GCM and SSP, and refugia = 1
is_refugia = refugia_lookup.GCM == OPTIONS.GCM & refugia_lookup.SSP == OPTIONS.SSP & refugia_lookup.refugia == 1;

% Set the reefs that are thermal refugia for this scenario
META.refugia = refugia_lookup.ReefID(is_refugia);

% Filter lookup table for this GCM and SSP, and hotspot = 1
is_hotspot = refugia_lookup.GCM == OPTIONS.GCM & refugia_lookup.SSP == OPTIONS.SSP & refugia_lookup.hotspot == 1;

% Set the reefs that are hotspots for this scenario
META.hotspot = refugia_lookup.ReefID(is_hotspot);

%% CoTS risk rank (coral cover lost to CoTS from counterfactual across six climate futures) for each reef to determine control priority - Tina Jan 2026
load('CoTSriskrank_lookup.mat'); % aggregated across time 2025-2075
CoTSreefloss = readtable('CoTSReefLoss_ByYear_ByClimateFuture'); % yearly for regional effort allocation
CoTSreefloss = convertvars(CoTSreefloss, ["GCM", "SSP"], "string");

% Filter for the current GCM and SSP
is_current_rank = CoTSriskrank.GCM == OPTIONS.GCM & CoTSriskrank.SSP == OPTIONS.SSP;
is_current_loss = CoTSreefloss.GCM == OPTIONS.GCM & CoTSreefloss.SSP == OPTIONS.SSP;

% Extract the reef IDs and their ranks for this scenario
META.cots_risk_rank = CoTSriskrank.CoTSRiskRank(is_current_rank);
META.cots_reef_loss = CoTSreefloss(is_current_loss,:);

%% Key source reef (KSR) reef sizes (if running scenarios 25-30 and using KSR ranking then need reef sizes) - Tina Feb 2026
load('KSR_reefsizes.mat');

% Set the CoTS risk ranking of all reefs 
META.KSR_reefsizes = KSR_reefsizes;

%% CoTS benefits 2025 - 2075. Suki March 2026. reef_coral_cover_deltas provided by Tina (see also Skinner et al. 2025)
META.cots_benefits = readtable('reef_coral_cover_deltas.csv');
META.cots_benefits.t = 23 + 2 * (META.cots_benefits.Year - 2019); % convert year to timestep
META.cots_benefits.SSP = extractAfter(META.cots_benefits.SSP, 3);

% Filter for the current GCM and SSP
is_current = META.cots_benefits.GCM == OPTIONS.GCM & META.cots_benefits.SSP == OPTIONS.SSP;

% Extract the reef IDs and their reef coral cover deltas for this scenario
META.cots_benefits = META.cots_benefits(is_current, :);

%% Prioritisation options
META.COTS_cull_fixed_reeflist = 0; % Specifies whether boat order for reef visitation is fixed (1) or randomised (0) at every time step
META.min_control_cover = 20; % minimum coral cover needed for control to still happen when high COTS in f_make_ReefList_TS
META.max_COTS = 3; % Specifies max COTS per tow above which control wouldn't happen as too many.


%% Workable reefs within GBRMP
META.workable = readtable('GBR_REEFS_GBRMP_Workable_SL.csv');

%% Nearest port distances between adjacent regions (for inter-port travel cost when redistributing effort)
% Columns: IN_REGION, NEAR_REGION, IN_FID, NEAR_FID, NEAR_DIST_KM (regions: 1=FN, 2=N, 3=C, 4=S)
META.port_distances = readtable('nearest_ports.csv');

%% TEMP FOR DEBUGGING
% META.COTS_fixed_list = 1; % required for each boat (consider deleting, might be redundant with META.COTS_cull_fixed_reeflist)
% META.cntrl_sites = [META.reef_ID META.COTS_cull_reeflist.nb_sites];
% META.COTS_control_sites = META.cntrl_sites;
% META.cntrl_reefID = META.reef_ID;

%% CONTROL PROGRAM EFFORT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.COTS_cull_boats = 5; % Specify the available effort in number of boats and days per boat; develop surveys later
META.COTS_cull_days = 90;  %Number of boat days at sea - 100 days/6months per AMPTO; some also lost to travel
META.COTS_cull_voyages = 9;  %number of discrete voyages per boat per 6 months - 13 per AMPTO - reduced to 9 to align with effort from actual Control Program
META.divers = 8; % number of divers per vessel
META.diven = 4; % number of dives per day

%META.COTS_pref_coral_groups=1:4;%which coral groups are taken into account for ET
% META.COTS_postcontrol_proportions=[ 1 1 (1-META.COTS_detectability(3:end))./(sum(1-META.COTS_detectability(3:end)))];%this is the population structure at ET after control; based on detectability, i.e. survivng population=1-detectability

%ecological threshold above which a reef must be to be treated - not currently used
%META.COTS_ecological_threshold=0.22;

%%% GIVING BOATS PROPERTIES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign properties to each boat
boats = META.COTS_cull_boats; %Specify boats using number of cull boats

for i = 1:boats %here you could remove the loop, if numb is a scalar integer
    %identifyTheBoat
    META.boatProperties.boatID = 1:boats ;
    %alotted boat days per 6 months
    META.boatProperties.boatDays = repmat(META.COTS_cull_days,1,boats);%
    %alotted voyages per 6 months
    META.boatProperties.voyages = repmat(META.COTS_cull_voyages,1,boats);%
    %number of divers onboard 
    META.boatProperties.divers = repmat(META.divers,1,boats); %
    %boat visit reefs in specific fixed order (1) of their ranking, or choose randomly (0) from the top X/reefs2cull reefs where to go first
    % META.boatProperties.fixedOrder=repmat(META.COTS_fixed_list,1,boats);
end

%%% CULLING EFFORT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating culling effort for each vessel

META.control_effort_allocation = 0.9; %Proportion of control effort to go to culling after some allocated to mantatowing/RHIS. Here 10% to manta towing/RHIS

for i=1:boats %for each boat
        
    META.boatProperties.totalTeamDives(i)=META.boatProperties.boatDays(i)*META.diven*META.control_effort_allocation;%total dives a team can make; each should be on a new site
    META.boatProperties.totalInidvDives(i)=META.boatProperties.boatDays(i)*META.diven*META.boatProperties.divers(i)*META.control_effort_allocation;%tdives of individual divers; note that this assumes divers from the same boat can all be on different sites on the same reef during the same dive which is probably unrealistic

end

META.max_dives_per_site = 300; %For high effort reefs, specify stopping rule threshold number of dives at cull site level - for hours convert * (40/60) = 200 hours

META.max_dives_per_reef = 3000; % For high effort reefs, specify stopping rule threshold number of dives at reef level - for hours convert * (40/60) = 2000 hours

% Suki April 2026 WIP - static effort allocation across region.
% Compute regional effort allocation proportional to mean CoTS loss area (counterfactual 2025-2075) per region.
% Regions: 1=Far Northern (FN), 2=Cairns/Cooktown (N), 3=Townsville/Whitsunday (C), 4=Mackay/Capricorn (S)
% Reefs with AREA_DESCR == 'out' (outside Marine Park) are excluded from this calculation.
cots_risk_current = CoTSriskrank(CoTSriskrank.GCM == OPTIONS.GCM & CoTSriskrank.SSP == OPTIONS.SSP, :);

% Match cots_risk_current.ReefID to GBR_REEFS.REEF_ID to get the positional index in GBR_REEFS
[~, gbr_idx] = ismember(cots_risk_current.ReefID, GBR_REEFS.Reef_ID);
% Keep only reefs that are part of this simulation (META.reef_ID)
in_sim = ismember(gbr_idx, META.reef_ID);
% Look up AREA_DESCR and loss values for simulated reefs only
sim_area_descr = GBR_REEFS.AREA_DESCR(gbr_idx(in_sim));
sim_loss = cots_risk_current.mean_COTS_loss_area(in_sim);

% group by regions avaialble in the simulation and sum the loss
[G, region_names] = findgroups(sim_area_descr);
region_loss = splitapply(@sum, sim_loss, G);

region_order = categorical({ ...
    'Far Northern', ...
    'Cairns/Cooktown', ...
    'Townsville/Whitsunday', ...
    'Mackay/Capricorn'});

region_loss_ordered = zeros(4,1);

for i = 1:4
    idx = region_names == region_order(i);
    if any(idx)
        region_loss_ordered(i) = region_loss(idx);
    end
end

% Adjust weighting based on current capacity (Reef Authority)
region_loss_ordered = region_loss_ordered .* [0.1; 0.5; 0.2; 0.2];

META.regional_effort_allocation = (region_loss_ordered ./ sum(region_loss_ordered))'; % proportion of total CoTS risk per region [FN, N, C, S]
% when aggregated across time = [0.338 0.242 0.163 0.257]
% According to Dave, N region currently holds ~50% of total capacity
% (highest among regions) - which isn't reflective of CoTS risk.

% Toggle for effort allocation method:
%   0 = static  - proportional to mean CoTS loss area (counterfactual, computed above)
%   1 = dynamic - recomputed every timestep in f_runmodel from current CoTS densities
META.effort_alloc_method = 1;

META.use_regional_dives = 1;   % 1 = enforce per-region dive budgets (new); 0 = global budget only, original behaviour
META.use_site_skip_find = 1;   % 1 = skip expensive sites and try cheaper ones on same reef (new); 0 = break site loop on first budget hit, original behaviour

%%% NOT IMPLEMENTED YET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify whether full state of the system is known; implement partial knowledge later - not currently developed
%META.COTS_state_known=1; 

% cull surrounding reefs - not currently developed
%META.COTS_cull_surrounding = 0; 

%Number of survey boats - not developed yet
%META.COTS_survey_boats = 0;  

%Number of survey voyages - not developed yet
%META.COTS_survey_voyages = 0;  

%Is this vessel a cull boat (1) or a survey vessel (0) currently not developed
%META.boatProperties.mission = repmat(1,1,boats);%randi([0 1],1,boats);

%the location on the coast that the vessel is stationed: 1 = Port Douglas, 2 = Cairns, 3 = Townsville, 4 = Mackay
%homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extend this vector if your wish to have more than 8 vessels
%homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extended for modelling up to 40 vessels
%META.boatProperties.homePort = homeports(1:boats);%randi([1 4],1,boats); Doesn't seem to be used?

%META.boatProperties.maxDays_atSea = repmat(13,1,boats);%   %maximum days the vessel can be at sea at one time

%%% OBSOLETE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No record of these parameters throughout the model, have commented them and left them here for now in case 

%META.listFileName = 1;   

%choose control strategy - obsolete? Not used anymore it seems.
%META.COTS_control_strat = 3;   

%picking the highest priority reef
%META.top_reef_picks=1;  %What is this?
%META.COTS_top_reef_picks = 1;  

% META.calculate_effort = 2;   %Not sure what this is ??

%effort quota, or area cleaned by boat per 6 months; 11.52km2 per 6
%months * proportion effort available
%META.boatProperties.effortQuota = zeros(1, length(META.boatProperties.boatID));
%META.boatProperties.effortQuota(i) = area_day*boat_days*META.control_effort_allocation;

%average distance of 20 min swim in m; calculated from timed swims, assume AMPTO moves at same speed, although culls probably slower
%swimd=480;
%average width covered during swim in m; manual reach; no changes due to habitat complexity etc, no slowdown due to high densities
%swimw=3;
%average area covered by a diver per hour in m2; 4320m2, or 66x66m
%swima=swimd*swimw*3;
%number of hours dived per dive; AMPTO
%divet=2/3;%was 2/3 originally, maybe still is?
%area covered by boat per day; 115200m2, or 340x340m, or 0.1152km2
%area_day=swima*META.divers*divet*diven;
%META.boatProperties.areaPerDay(i) = area_day;
