%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tina Skinner, MSEL, June 2023.
%
% Parametrisation of COTS control settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

META.COTS_reef_trigger = 0; % whether (1) or not (0) the trigger of control is reef-level CoTS above threshold of 0.22 CoTS per tow
% If set to 0, the trigger is CoTS density per site, not per reef. The alternative (1) hasn't been tested yet.

META.COTS_control_start = 23; % timestep in 6-month intervals when control should start = 23 (summer 2019)

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

META.doing_COTS_heat_stress_detectibility = 1; %Uses evidance from Cook et al. in review

if META.doing_COTS_heat_stress_detectibility == 1 
    load("log_culling_dhw_effects.mat"); %Results from Cook et al. 2025 - How does culling change under DHW per size class (roughly mapped from the 4 size classes to the 16 here)
    META.log_culling_dhw_effects = log_culling_dhw_effects; %Note: will only work with 16 size classes
end

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


% Updates from Tina (July 2023) - now have fixed target reef list. Control at T, then P, then N.
load('New_regions_TS.mat') %this has been updated for new GBRMPA 2023 PR list and includes target reefs now
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

%% Prioritisation options
META.COTS_cull_fixed_reeflist = 0; % Specifies whether boat order for reef visitation is fixed (1) or randomised (0) at every time step
META.min_control_cover = 20; % minimum coral cover needed for control to still happen when high COTS in f_make_ReefList_TS
META.max_COTS = 3; % Specifies max COTS per tow above which control wouldn't happen as too many.

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
META.boatProperties.totalInidvDives(i)=META.boatProperties.boatDays(i)*META.diven*META.boatProperties.divers(i)*META.control_effort_allocation;%tdives of individual divers; note that this assumes divers fromt eh same boat can all be on different sites on the same reef durign the same dive which is probably unrealistic

end

META.max_dives_per_site = 300; %For high effort reefs, specify stopping rule threshold number of dives at cull site level - for hours convert * (40/60) = 200 hours

META.max_dives_per_reef = 3000; % For high effort reefs, specify stopping rule threshold number of dives at reef level - for hours convert * (40/60) = 2000 hours

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
