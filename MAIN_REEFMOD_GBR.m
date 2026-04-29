%__________________________________________________________________________
% 
% REEFMOD-GBR MAIN SCRIPT, version 7.3
% Yves-Marie Bozec, y.bozec@uq.edu.au (Mar 2025 - Nov 2025)
%
% Extensive revision of CoTS control (f_COTS_control_NEW) with correct conversion CoTS density <-> per tow and correct
% calculation of culled density relative to the ecological threshold. Also revised entirely the computation of priority
% reef lists (makeReefList_NEW) with contracted and optimsed code. 
% CoTS control records now stored in 'RECORD' instead of 'RESULT'.
%
% Added size-specific background mortality for corals (mortality decreasing with size). 
%
% New calibration of coral stock-recruitment (beta parameter) under new assumption for the default fertilisation success
% FS = 0.4 when Allee effects is OFF (consistent with the maximum FS when Allee effects is ON).
% Allee effects in coral fertilisation is OFF by default.
% 
% Refined the selection of dispersal model for coral and CoTS (settings_CONNECTIVITY_NEW_2) with a dedicated CONNECT
% structure array for coral and another for CoTS. Connectivity matrices now averaged over multiple spawning events 
% of the same year to optimise storage and standardise across dispersal models.
% Recommended settings (default): 
% -> GBR1 for coral (but can be swapped back to GBR4 without modifying BH_alpha/BH_beta). 
% -> GBR4 for CoTS; GBR1 is available but not recommended. GBRLUP also available but requires new parametrisation of
% CoTS mortality, BH_alpha and BH_beta (ongoing work by Suki).
% Also refined the tracking of sequence of connectivity matrices (in META.connectivity and RECORD). 
% 
% Outputs now compiled and formatted in a dedicated script (FORMAT_OUTPUTS) which is run at the end of MAIN_REEFMOD.
% Outputs can be produced either for each seasonal step (6 months) or yearly (choice to be set in format_extract).
%__________________________________________________________________________
clear

SaveDir ='';

NB_SIMULATIONS = 1; % Number of repeated runs

% NB_TIME_STEPS has to be an even number. 
% Always run the hindcast before future projections (initialisation = winter 2007)
% Example: 32 (hindcast 2008-2023) + 154 (forecast 2024-2100)

% Historic cyclones and DHW currently available until 2024 (inclusive). 
NB_TIME_STEPS = 34+52; % HINDCAST: summer 2008 to winter 2024 (34 steps)
% NB_TIME_STEPS = 34+152; % HINDCAST+FORECAST summer 2008 - winter 2100

% Set the format for output extraction and saving
format_extract = 'short' ; % annual outputs (every year)
% format_extract = 'long' ; % seasonal outputs (every 6 months) (native time resolution)

% Set the name for the output file. Will add suffix 's' for 'short' or 'l' for 'long'
OutputName = 'testing_dynamic'; options = [1 1 1 1 0 1 0.3 0]; % see list of options below

%% select the Global Circulation Model for climate change projection (CMIP-6)
GCM = 1; % 1=CNRM-ESM2-1, 2=EC-Earth3-Veg, 3=IPSL-CM6A-LR, 4=MRI-ESM2-0, 5=UKESM1-0-LL, ...
% 6=GFDL-ESM4, 7=MIROC-ES2L, 8=MPI-ESM1-2-HR, 9=MIROC6, 10=NorESM2-LM

SSP = 3; % 1=SSP1-1.9, 2=SSP1-2.6, 3=SSP2-4.5, 4=SSP3-7.0, 5=SSP5-8.5
% Note SSP1-1.9 is not available for GFDL-ESM4, MPI-ESM1-2-HR and NorESM2-LM

%% --------------------------------------------------------------------------------
% Climate forecasts from CMIP6 global circulation models
All_GCMs = ["CNRM-ESM2-1" ; "EC-Earth3-Veg" ; "IPSL-CM6A-LR" ; "MRI-ESM2-0" ; "UKESM1-0-LL" ; ...
    "GFDL-ESM4" ; "MIROC-ES2L" ; "MPI-ESM1-2-HR" ; "MIROC6" ; "NorESM2-LM" ];
All_SSPs = ["119" ; "126" ; "245" ; "370" ; "585" ];

OPTIONS.GCM = All_GCMs(GCM);
OPTIONS.SSP = All_SSPs(SSP);

% Stressor options: yes(1)/no(0)
OPTIONS.doing_cyclones = options(1);
OPTIONS.doing_bleaching = options(2) ;
OPTIONS.doing_COTS = options(3);
OPTIONS.doing_WQ = options(4);
OPTIONS.doing_restoration = options(5) ;
OPTIONS.doing_COTS_control= options(6);
OPTIONS.heritability_HT = options(7); % heritability of heat tolerance for all groups (inheritance model of HT phenotypes)
OPTIONS.allee_effect = options(8);

OPTIONS.doing_size_frequency = 1; % for tracking population size structure (incl. juveniles)
OPTIONS.doing_genetics = 0 ; % for running the genetic model (Bozec & Mumby 2019) - needs revision
OPTIONS.genetic_parms = [ 1 1 1.5 ]; % sigma cold/sigma hot/esd (Bozec & Mumby 2019) - useless if doing_genetics = 0

%% --------------------------------------------------------------------------------
% Below options are used to force simulations with specific starting conditions - keep empty if of no use
OPTIONS.init_coral_cover = []; %0.01*ones(1,6); % as proportional cover (vector of 6 values)
OPTIONS.init_sand_cover = []; % as proportional cover
OPTIONS.init_rubble_cover = []; % as proportional cover
OPTIONS.ssc = []; %in mg/L

%% --------------------------------------------------------------------------------
% RESTORATION INTERVENTIONS (high level controls; set parameters in settings_RESTORATION)

% 1) Outplanting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of outplanting (set to 0 for not doing this intervention)
RESTORATION.doing_coral_outplanting = uint8(zeros(1,NB_TIME_STEPS)); % if zero, no coral deployment at time step for all reefs
RESTORATION.doing_coral_outplanting(1,38:2:46) = 1 ; % set to 1 to indicate outplanting: first in summer 2026 and last in summer 2030 (BCA)
% RESTORATION.doing_coral_outplanting(1,[38:2:46 78:2:86 118:2:126 158:2:166] ) = 1 ; % outplanting starts 2026-2030, 2046-2050, 2066-2070, 2086-2090

RESTORATION.total_nb_outplants = 1e6; % Nb of outplants available for the GBR at each deployment (= time step). Set to Inf if density/nb_reefs are the drivers
RESTORATION.outplanted_density = 6.8;  % only in the case of fixed density of outplants (ignores RESTORATION.total_nb_outplants?)
RESTORATION.MEAN_HT_outplants = 0*[1 1 1 1 1 1];% as +°C-week (DHW) compared to the contemporary mean of the group

% 2) Larval enrichment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of larval enrichment (set to 0 for not doing this intervention)
RESTORATION.doing_larval_enrichment = uint8(zeros(1,NB_TIME_STEPS));
% RESTORATION.doing_larval_enrichment(1,38:2:46) = 1 ; % set to 1 to do outplanting
% RESTORATION.doing_larval_enrichment(1,2:2:30) = 1 ;

RESTORATION.nb_reefs_enriched = 10 ; % max number of reefs where larval enrichment is undertaken at each time step
% Will be picked from a priority list generated in f_generate_priority_list OR specified in SETTINGS_RESTORATION)
% If 0, enrichment cannot happen. For the counterfactual, set to Inf with total_nb_larvae = 0 for ghost deployment
RESTORATION.total_nb_larvae = Inf; % Max number of 'larvae' (ie, 1 yr old corals) available at each time step. Set to Inf if unlimited.
RESTORATION.MEAN_HT_larvae = 2*[1 1 1 1 1 1];% as +°C-week (DHW) compared to the contemporary mean of the group

% 3) Rubble stabilisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the restoration effort: number of reefs where rubble is stabilised at each time step
RESTORATION.nb_reefs_stabilised = 0 ; % (if 0 rubble stabilisation cannot happen)
% Set the timing of intervention (if 1 intervention is deployed at step t, if 0 no intervention at t)
RESTORATION.doing_rubble_stabilisation = uint8(zeros(1,NB_TIME_STEPS));
% RESTORATION.doing_rubble_stabilisation(1,38:1:end) = 1 ; % set to 1 to do outplanting: first in summer 2026 and last in summer 2030

% 4) Solar Radiation Management  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fogging for RRAP 2022 intervention simulations
RESTORATION.doing_fogging = uint8(zeros(1,NB_TIME_STEPS)); %if zero don't do fogging at time step
% RESTORATION.doing_fogging(1,32:2:end) = 1 ; % set to 1 to do fogging (only in summer)

RESTORATION.nb_reefs_fogged = []; % number of reefs if random deployment, otherwise =[] but need to listthe  reef IDs
RESTORATION.fogged_reef_ID = [695 697 969 970]; % using 5 fogging units (22.5 km2) covering 23 km2 (equivalent 2D reef areas)
% 695-Moore; 697-Elford; 698-Briggs ; 969-Milln; 970-Thetford
% Leave empty otherwise []

% Cloud brightening (=cooling) from Bozec & Mumby (2019) (needs revision)
RESTORATION.doing_cooling = uint8(zeros(1,NB_TIME_STEPS)); %if zero don't do cooling at time step
% RESTORATION.doing_fogging(1,32:2:end) = 1 ; % set to 1 to do cooling (only in summer)

%% --------------------------------------------------------------------------------
%% Create output file name
if NB_TIME_STEPS > 34
    OPTIONS.OutputFileName = [SaveDir format_extract(1) OutputName  '_SSP' char(OPTIONS.SSP) '_' char(All_GCMs(GCM)) '.mat'];
    % Display the running scenario
    ['Running scenario ' OutputName '_herit' num2str(OPTIONS.heritability_HT) '_SSP' char(OPTIONS.SSP) '_' char(All_GCMs(GCM)) ' .....']
else
    OPTIONS.OutputFileName = [SaveDir format_extract(1) OutputName  '_HINDCAST.mat' ];
    % Display the running scenario
    ['Running scenario ' OutputName '_herit' num2str(OPTIONS.heritability_HT) '_HINDCAST'  ' .....' ]
end

%% --------------------------------------------------------------------------------
%% RUN REEFMOD
OUTPUTS = struct('REEF', [],'RESULT', [],'RECORD', []);
TEMP_META = struct('META', []);

% parfor run_id = 1:NB_SIMULATIONS
for run_id = 1:NB_SIMULATIONS

    run_id

    tic
    [meta, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, NB_TIME_STEPS, run_id);
    toc

    OUTPUTS(run_id).RESULT = RESULT ;
    OUTPUTS(run_id).RECORD = RECORD ;
    OUTPUTS(run_id).REEF = REEF ;
    TEMP_META(run_id).META = meta ;
end

META = TEMP_META(1).META; % Keep only one META because common to all simulations
clear TEMP_META ADAPT run_id meta REEF RESULT RECORD GCM GCM_list options SaveDir OutputName RCP RCP_list

%% --------------------------------------------------------------------------------
%% FORMAT, OPTIMISE & SAVE
FORMAT_OUTPUTS % new script that extract all outputs for each run

% Possibility here to delete specific outputs before saving
clear NB_TIME_STEPS ...
    reef_shelter_volume_absolute ...
    selection_diff ... % selection differential = difference in mean HT between start/end of the year
    COTS_densities ... % CoTS density-at-age
    COTS_predicted_densities ... % predicted CoTS density-at-age before erased by observation (hindcast)
    macroUprightFleshy macroEncrustFleshy macroTurf % algal covers (only worth keeping if grazing < 1)

save (OPTIONS.OutputFileName)
