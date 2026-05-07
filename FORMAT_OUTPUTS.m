%______________________________________________________________________________________________
% 
% FORMAT REEFMOD OUTPUTS
% Yves-Marie Bozec, y.bozec@uq.edu.au (Nov 2025)
%
% Gathers the code previously in MAIN_REEFMOD_GBR below the running model command.
% Now provides 2 different formats for outputs (to be set in MAIN_REEFMOD_GBR):
% format_extract = 'short', which compiles model outputs on a yearly basis
% format_extract = 'extended', which compiles model outputs on a 6-month basis
%______________________________________________________________________________________________

load('GBR_REEF_POLYGONS_2026.mat') % new habitat areas based on geomorphic map

% Scale according to chosen format
switch format_extract

    case 'long' % for outputs every 6 months (native time resolution)
        scale = 0;
        TIME = NB_TIME_STEPS;

    case 'short' % for annual outputs (yearly averaged outputs)
        scale = 1;
        TIME = (NB_TIME_STEPS)/2;
end

% cols = 2:2:NB_TIME_STEPS; % column indices for averaging
cols = 2:(1+scale):(NB_TIME_STEPS+1-scale); % column indices for averaging


% Initial step is end (winter) of 2007 (=2007.5), first step is summer 2008 (=2008.0), last is end (winter) of 2017 (=2020.5)
start_year = 2007.5 ;
YEARS =  start_year + (0:NB_TIME_STEPS)/2 ;
YEARS = YEARS(2:2:end) ;

% Extract the defintiion of the simulated reefs
myGBR_REEFS = GBR_REEFS(META.reef_ID,[1 3:7 10:16]);

%% Memory allocation
% 1) coral outputs
coral_cover_per_taxa = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types,'single');
coral_larval_supply = coral_cover_per_taxa;
total_shelter_volume_per_taxa = coral_cover_per_taxa;
coral_HT_mean = nan(NB_SIMULATIONS, META.nb_reefs, TIME+1,META.nb_coral_types,'single');
coral_HT_var = coral_HT_mean;
selection_diff = coral_HT_mean;

nb_coral_offspring = coral_cover_per_taxa;
nb_coral_recruit = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types,'uint16');

if OPTIONS.doing_size_frequency == 1
    % Bin edges for juvenile bins (cm): 1  3  5 (2 size classes)
    % Adolescents (cm):  5  9  13  17 (3 size classes)
    % Adults (cm): 17 31 45 59 73 87 101 (6 size classes)
    nb_coral_juv = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types, 2, 'uint16') ;
    nb_coral_adol = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types, 3, 'uint16') ;
    nb_coral_adult = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types, 6, 'uint16') ;
end

% 2) Coral loss for each stressor (need to add initial step)
coral_cover_lost_bleaching = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, META.nb_coral_types, 'single');
coral_cover_lost_cyclones = coral_cover_lost_bleaching;
coral_cover_lost_COTS = coral_cover_lost_bleaching;
coral_mortality_bleaching = coral_cover_lost_bleaching;

% 3) Stress records
record_applied_cyclones = zeros(NB_SIMULATIONS, META.nb_reefs, TIME,'uint8');
record_applied_DHWs = zeros(NB_SIMULATIONS, META.nb_reefs, TIME,'single');
record_bleaching_category = zeros(NB_SIMULATIONS, META.nb_reefs, TIME,'single');

% 4) Connectivity sequence records (NOW IN META.connectivity for a single run but need to be collected for all runs)
% Note index values were generated every time steps, yet spawning only happens in summer, so selection starts with the first value then every 2 values
record_spawning_chronology_CORAL = zeros(NB_SIMULATIONS, TIME,'uint8'); % sequence of connectivity years for corals (depends on dispersal model)
record_spawning_chronology_COTS = zeros(NB_SIMULATIONS, TIME,'uint8'); % sequence of connectivity years for COTS (depends on dispersal model)

% 5) Restoration records (TO BE TESTED)
if OPTIONS.doing_restoration==1

    if sum(sum(OUTPUTS(1).RECORD.outplanted_reefs))> 1
        % Total nb of outplants per reef per time step
        record_total_outplants_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, TIME, META.nb_coral_types,'uint16');
        record_outplanted_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
    end

    if sum(sum(OUTPUTS(1).RECORD.enriched_reefs))> 1
        record_total_larvae_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, TIME, META.nb_coral_types,'uint16');
        record_enriched_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
    end

    % record_rubble_pct2D_stabilised = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'single');
    % record_stabilised_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
    % record_fogged_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');

    coral_cover_per_taxa_restored_sites = zeros(NB_SIMULATIONS,META.nb_reefs,TIME+1,META.nb_coral_types,'single');

else
    clear RESTORATION
end

% 6) Other variables
rubble = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
nongrazable = zeros(NB_SIMULATIONS, META.nb_reefs,'uint8');
macroTurf = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
macroEncrustFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');
macroUprightFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1,'uint8');

% 7) CoTS outputs
if OPTIONS.doing_COTS == 1
    COTS_mantatow = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, 'single');
    COTS_densities = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, 16, 'single'); % 16 age classes
    COTS_settler_density = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, 'uint16');
    COTS_larval_supply = COTS_mantatow;
    COTS_larval_output = COTS_mantatow;

    % If recording CoTS densities predicted by the model before being replaced by observations (only during hindcast) 
    % COTS_predicted_densities = zeros(NB_SIMULATIONS, META.nb_reefs, TIME+1, 16, 'single'); % 16 age classes
end

%% Populate outputs
% Anonymous function to average over time for the dimensional outputs
average_time_3D = @(X, scale, cols) cat(2, X(:,1,:), (X(:,cols,:) + X(:,cols+mod(scale,2),:))/2);
average_time_4D = @(X, scale, cols) cat(2, X(:,1,:,:), (X(:,cols,:,:) + X(:,cols+mod(scale,2),:,:))/2);
sum_time_3D = @(X, scale, cols) cat(2, X(:,1,:), (X(:,cols,:) + X(:,cols+mod(scale,2),:)));
sum_time_4D = @(X, scale, cols) cat(2, X(:,1,:,:), (X(:,cols,:,:) + X(:,cols+mod(scale,2),:,:)));
% cols = 2:scale:TIME; % column indices for averaging

for simul = 1:NB_SIMULATIONS

    % 1) coral outputs
    tmp_coral_cover_per_taxa = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D));
    tmp_coral_larval_supply = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_larval_supply)); % nb of incoming larvae per unit of reef area (400m2)
    tmp_nb_coral_recruit = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_settler_count));
    tmp_nb_coral_offspring = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_total_fecundity)); % nb of larvae produced per unit of reef area (400m2)
    tmp_total_shelter_volume_per_taxa = squeeze(cat(4,OUTPUTS(simul).RESULT.total_shelter_volume_per_taxa));
    tmp_coral_HT_mean = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_HT_mean));
    tmp_coral_HT_var = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_HT_var));
    tmp_selection_diff = squeeze(cat(4,OUTPUTS(simul).RESULT.selection_diff));
  
    coral_cover_per_taxa(simul,:,:,:) = average_time_3D(tmp_coral_cover_per_taxa, scale, cols);
    coral_larval_supply(simul,:,:,:) = average_time_3D(tmp_coral_larval_supply, scale, cols);
    nb_coral_recruit(simul,:,:,:) = average_time_3D(tmp_nb_coral_recruit, scale, cols);
    nb_coral_offspring(simul,:,:,:) = average_time_3D(tmp_nb_coral_offspring, scale, cols);
    total_shelter_volume_per_taxa(simul,:,:,:) = average_time_3D(tmp_total_shelter_volume_per_taxa, scale, cols);
    coral_HT_mean(simul,:,:,:) = average_time_3D(tmp_coral_HT_mean, scale, cols);
    coral_HT_var(simul,:,:,:) = average_time_3D(tmp_coral_HT_var, scale, cols);
    selection_diff(simul,:,:,:) = average_time_3D(tmp_selection_diff, scale, cols);

    if OPTIONS.doing_size_frequency == 1
        tmp_nb_coral_juv = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_juv_count(:,:,:,:)));
        tmp_nb_coral_adol = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adol_count(:,:,:,:)));
        tmp_nb_coral_adult = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adult_count(:,:,:,:)));

        nb_coral_juv(simul,:,:,:,:)= average_time_4D(tmp_nb_coral_juv, scale, cols);
        nb_coral_adol(simul,:,:,:,:) = average_time_4D(tmp_nb_coral_adol, scale, cols);
        nb_coral_adult(simul,:,:,:,:) = average_time_4D(tmp_nb_coral_adult, scale, cols);
    end

    % 2) Coral cover loss following stressors
    tmp_coral_cover_lost_bleaching = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_bleaching);
    tmp_coral_cover_lost_cyclones = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_cyclones);
    tmp_coral_cover_lost_COTS = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_COTS);
    tmp_coral_mortality_bleaching = squeeze(OUTPUTS(simul).RECORD.bleaching_mortality);

    coral_cover_lost_bleaching(simul,:,2:end,:) = tmp_coral_cover_lost_bleaching(:,1:(scale+1):end,:);
    coral_cover_lost_cyclones(simul,:,2:end,:) = tmp_coral_cover_lost_cyclones(:,1:(scale+1):end,:);
    coral_cover_lost_COTS(simul,:,2:end,:) = tmp_coral_cover_lost_COTS(:,1:(scale+1):end,:);
    coral_mortality_bleaching(simul,:,2:end,:) = tmp_coral_mortality_bleaching(:,1:(scale+1):end,:);

    % 3) Stress records
    tmp_record_applied_cyclones = squeeze(OUTPUTS(simul).RECORD.hurricane_events);
    tmp_record_applied_DHWs = squeeze(OUTPUTS(simul).RECORD.applied_DHWs);
    % Suki March 2026 add GBRMPA bleaching categories
    tmp_record_bleaching_category = squeeze(OUTPUTS(simul).RECORD.bleaching_category);

    record_applied_cyclones(simul,:,:) = tmp_record_applied_cyclones(:,1:(scale+1):end);
    record_applied_DHWs(simul,:,:) = tmp_record_applied_DHWs(:,1:(scale+1):end);
    record_bleaching_category(simul,:,:) = tmp_record_bleaching_category(:,1:(scale+1):end);

   % 4) Connectivity sequence records  
    record_spawning_chronology_CORAL(simul,:) =  cell2mat(OUTPUTS(simul).RECORD.spawning_chronology_CORAL(2,1:(scale+1):end));
    record_spawning_chronology_COTS(simul,:) =  cell2mat(OUTPUTS(simul).RECORD.spawning_chronology_COTS(2,1:(scale+1):end));

    % 5) Restoration records
    if  OPTIONS.doing_restoration==1

        if sum(sum(OUTPUTS(simul).RECORD.outplanted_reefs))> 1
            record_total_outplants_deployed(simul,:,1:end,:) = squeeze(OUTPUTS(simul).RECORD.total_outplanted(:,1:(scale+1):end,:));
            record_outplanted_reefs(simul,:,1:end) = squeeze(OUTPUTS(simul).RECORD.outplanted_reefs(:,1:(scale+1):end));
        end

        if sum(sum(OUTPUTS(simul).RECORD.enriched_reefs))> 1
            record_total_larvae_deployed(simul,:,1:end,:) = squeeze(OUTPUTS(simul).RECORD.total_enriched(:,1:(scale+1):end,:));
            record_enriched_reefs(simul,:,1:end) = squeeze(OUTPUTS(simul).RECORD.enriched_reefs(:,1:(scale+1):end));
        end

        % record_rubble_pct2D_stabilised(simul,:,2:end) = squeeze(OUTPUTS(simul).RECORD.rubble_cover_pct2D_stabilised(:,1:end));
        % record_stabilised_reefs(simul,:,2:end) = squeeze(OUTPUTS(simul).RECORD.stabilised_reefs(:,1:(scale+1):end));
        % record_fogged_reefs(simul,:,2:end) = RECORD.fogged_reefs;

        tmp_coral_cover_per_taxa_restored_sites = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D_restored_sites));
        coral_cover_per_taxa_restored_sites(simul,:,:,:) = average_time_4D(tmp_coral_cover_per_taxa_restored_sites, scale, cols);
    end

    % 6) Other variables
    nongrazable(simul,:) = squeeze(cat(4,OUTPUTS(simul).REEF.nongrazable_substratum));

    tmp_rubble = squeeze(OUTPUTS(simul).RESULT.rubble_cover_pct2D);
    tmp_macroUprightFleshy = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,2)));
    tmp_macroEncrustFleshy = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,3)));
    tmp_macroTurf = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,4)));

    rubble(simul,:,:)= average_time_3D(tmp_rubble, scale, cols);
    macroUprightFleshy(simul,:,:)= average_time_4D(tmp_macroUprightFleshy, scale, cols);
    macroEncrustFleshy(simul,:,:)= average_time_4D(tmp_macroEncrustFleshy, scale, cols);
    macroTurf(simul,:,:)= average_time_4D(tmp_macroTurf, scale, cols);

    % 7) CoTS outputs
    if OPTIONS.doing_COTS == 1
        
        tmp_COTS_densities = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities); % Density for 400m2 of each age class
        %tmp_COTS_predicted_densities = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities_predicted); % Density for 400m2 (predicted before erasure by obs)
        tmp_COTS_settler_density = squeeze(OUTPUTS(simul).RESULT.COTS_settler_densities); % Density for 400m2
        tmp_COTS_larval_supply = squeeze(OUTPUTS(simul).RESULT.COTS_larval_supply); % Density for 400m2
        tmp_COTS_larval_output = squeeze(OUTPUTS(simul).RESULT.COTS_larval_output); % Density for 400m2

        COTS_densities(simul,:,:,:) = average_time_4D(tmp_COTS_densities, scale, cols);
        %COTS_predicted_densities(simul,:,:,:) = average_time_4D(tmp_COTS_predicted_densities, scale, cols);
        COTS_settler_density(simul,:,:,:) = average_time_3D(tmp_COTS_settler_density, scale, cols);
        COTS_larval_supply(simul,:,:,:) = average_time_3D(tmp_COTS_larval_supply, scale, cols);
        COTS_larval_output(simul,:,:,:) = average_time_3D(tmp_COTS_larval_output, scale, cols);

        % Estime CoTS per maqnta tow, assuming 0.6 CoTS per grid ~ 0.22 CoTS per tow
        % (0.22 per tow is equivalent to 1500 COTS per km2 (Moran & De'ath 92), so that 1 COTS per grid (400m2) is equivalent to 0.22*2500/1500
        % Suki May 2026: use new non-linear conversion by YM
        for t=1:size(COTS_densities,3)
            TMP = squeeze(COTS_densities(simul,:,t,:)).*META.COTS_detectability(ones(1, META.nb_reefs),:);
            COTS_mantatow(simul,:,t) = f_convert_CoTS_density_2_tow(sum(TMP(:,META.COTS_adult_min_age:end),2), META.total_area_cm2);
        end
    end

    % COTS control (not available yearly but for every 6 month step)
    if OPTIONS.doing_COTS_control == 1

        if NB_TIME_STEPS > META.COTS_control_start

            record_COTS_control(simul).control_records = OUTPUTS(simul).RECORD.control_records;
            % regional_effort_allocation is stored per timestep in record_COTS_control(simul).control_records(t).regional_effort_allocation
            % New way of storing control records, for each simulation (separetely) at each time step (separetely), including:
            % - culled_reef_ID: list of reef ID that were culled (the last reef ID may have been partially culled)
            % - culled_density: total density (per 400 m2) of CoTS adults killed per culled reef
            % - nb_dives: number of dives simulated at each culled reef
            % - nb_culled_sites: number of sites culled at each culled reef
            % - nb_culled_reefs: total number of culled reefs
            % - nb_visited_reefs: total number of visited reefs (including culled reefs)
        end
    end
end

% Express shelter volume as total across all coral group in dm3 per m2
reef_shelter_volume_absolute = sum(total_shelter_volume_per_taxa,4)/(META.total_area_cm2/1e4);

% Express shelter volume relative to maximum possible value (ie, 1 tabular coral of the size of the cell in every cell)
total_shelter_volume_max = META.grid_x_count*META.grid_y_count*exp(-8.32 + 1.50*log(META.cell_area_cm2));
reef_shelter_volume_relative = sum(total_shelter_volume_per_taxa,4)/total_shelter_volume_max;
reef_shelter_volume_relative(reef_shelter_volume_relative>1)=1;
reef_shelter_volume_relative(reef_shelter_volume_relative<0)=0;

%% Further optimisation of storage
% New (08/2021): only record COTS densities by yearly classes to reduce output size
% COTS_densities = COTS_densities0(:,:,:,1:2:end)+COTS_densities0(:,:,:,2:2:end);
% 09/2021: now just sum across all juveniles and across all adults
if OPTIONS.doing_COTS == 1
    COTS_juv_densities = sum(COTS_densities(:,:,:,1:(META.COTS_adult_min_age-1)),4);
    COTS_adult_densities = sum(COTS_densities(:,:,:,META.COTS_adult_min_age:end),4);
end

% (11/2025): only record total nb of corals, not nb of corals persize class
nb_coral_juv = uint16(sum(nb_coral_juv,5));
nb_coral_adol = uint16(sum(nb_coral_adol,5));
nb_coral_adult = uint16(sum(nb_coral_adult,5));

% (11/2025): save coral cover lost as integrer values
coral_cover_lost_bleaching = uint8(coral_cover_lost_bleaching);
coral_cover_lost_cyclones = uint8(coral_cover_lost_cyclones);
coral_cover_lost_COTS = uint8(coral_cover_lost_COTS);

%% SAVE OUTPUTS
% delete all the tmp objects
clearvars -regexp ^tmp

clear ans cols scale s t TMP start_year...
    format_extract GBR_REEFS OUTPUTS simul NB_SIMULATIONS  ...
    All_GCMs All_SSPs GCM SSP TIME ...
    total_shelter_volume_max total_shelter_volume_per_taxa ...
    total_shelter_volume_max ...
    average_time_3D average_time_4D sum_time_3D sum_time_4D
