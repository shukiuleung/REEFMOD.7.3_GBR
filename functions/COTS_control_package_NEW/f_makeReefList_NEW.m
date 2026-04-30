% ----------------------------------------------------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, Mar 2025.
%
% Creates a list of ordered reefs for prioritising CoTS culling under different scenarios of control management.
% Optimised and revised version of the previous scripts:
% - f_makeReefList from Karlo Hock (REEFMOD-GBR.6.3, Sep 2019)
% - f_makeReefListCCS from Carolina Castro-Sanguino (REEFMOD-GBR.6.6, Mar 2022)
% - f_makeReefListTS from Tina Skinner (REEFMOD-GBR.6.8, May 2023)
% ----------------------------------------------------------------------------------------------------------------------
function [full_list_ID, unvisited, monitor_only_IDs, relocation_dist_km] = f_makeReefList_NEW(META, current_COTS_densities, current_reef_ET, COTS_densities_per_site,...
    total_coral_pct2D, COTS_larval_output, CONNECT_CORAL, t, last_reef_COTScontrolled, ...
    bleaching_category, severely_bleached, control_records_prev, ...
    current_COTS_per_tow)

% Tina 07/2023: now have fixed target reef list in 'New_regions_TS.mat', updated for new GBRMPA 2023 PR list.
% Control at target (T), then priority (P), then nonpriority reef (N), as specified in 'reef_type'.
% Region as FN, N, C, S -> no 'out' anymore (as in GBR_REEFS.AREA_DESCR), ie outside of the Marine Park, meaning that
% all reefs in the Top North are now included for control.

% YM 09/2025: 'New_regions_TS.mat' now replaces 'GBRMPAprioritylist.mat'. This canonical list is loaded in settings_COTS_CONTROL
% and stored in META.COTS_cull_reeflist with other reef characteristics ('Reef_ID', 'reef_type', 'Region', 'nb_sites', 'AIMS_sector', 'GreenZone').
% From this, run-specific lists of T, P and N are created in settings_COTS_CONTROL:
% META.COTS_cull_reeflist_targetRUN > META.COTS_cull_reeflist_priorityRUN > META.COTS_cull_reeflist_nonpriorityRUN
% Each list is randomised per run, ensuring reefs are visited in a different order. Lists are also filtered by the specified regions of intervention
% (META.COTS_cull_region), eg, only include Far North reefs if 'FN' was selected, allowing every strategy to be simulated within specific regions.
% ------- WHAT THIS SCRIPT DOES:
% 1) Set GBRMPA priorities as initially generated for the run, or further randomise for this time step (META.COTS_cull_fixed_reeflist).
% 2) Create the full list of reef ID based on the specified strategy (META.COTS_reefs2cull_strat)
% 3) Finally place on top of list the last reef controlled at the previous step, in case it was only paritially controlled

unvisited = []; % initialize the list of unvisited reefs, to keep track of which reefs in the priority list are not visited at each time step due to bleaching
monitor_only_IDs = []; % initialize list of reefs to be monitored only (not culled) due to regional bleaching
relocation_dist_km = struct('reef_ID', [], 'original_reef_ID', [], 'dist_km', []); % relocation travel distances: candidate reef ID, original (bleaching cat>3) reef ID, distance in km

%% Allow for permutation of reef prioritisation
if META.COTS_cull_fixed_reeflist == 1 % if the reef prioritisation list is fixed over time
    target_ID = META.COTS_cull_reeflist_targetRUN.Reef_ID;
    priority_ID = META.COTS_cull_reeflist_priorityRUN.Reef_ID;
    nonpriority_ID = META.COTS_cull_reeflist_nonpriorityRUN.Reef_ID;
else % otherwise shuffle the lists every time step
    target_ID = META.COTS_cull_reeflist_targetRUN.Reef_ID(randperm(length(META.COTS_cull_reeflist_targetRUN.Reef_ID)));
    priority_ID = META.COTS_cull_reeflist_priorityRUN.Reef_ID(randperm(length(META.COTS_cull_reeflist_priorityRUN.Reef_ID)));
    nonpriority_ID = META.COTS_cull_reeflist_nonpriorityRUN.Reef_ID(randperm(length(META.COTS_cull_reeflist_nonpriorityRUN.Reef_ID)));
end

%% Apply the specified prioritisation strategy
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
% 19 - Adaptive Control during bleaching = default (WIP)

switch META.COTS_reefs2cull_strat

    case 1 % GBRMPA strategy that goes to Target reefs first, then Priority reefs, then Non Priority reefs
        full_list_ID = vertcat(target_ID, priority_ID, nonpriority_ID); % just catenate the 3 lists

    case 9  % Outbreak front: GBRMPA strategy that goes to target reefs first, then also goes to 0.5' lat (~50 km)
        % from target reefs with outbreaks - whole GBR.

        %first, find all target reefs with outbreaking sites
        target_outbreaks=[];

        for rfs = 1:length(target_ID)
            this_reef_ID = target_ID(rfs);
            this_reef_COTS_densities_per_site = COTS_densities_per_site{this_reef_ID,1}; % Extract COTS density for all sites
            % Calculates manta tow equivalent (number of coTS per tow) for all sites
            this_reef_COTS_per_tow_per_site = (0.22/0.6)*sum(this_reef_COTS_densities_per_site(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
            % Find sites with manta tows above ecological threshold (ET)
            sites_over_ET = find(this_reef_COTS_per_tow_per_site > current_reef_ET(this_reef_ID));

            if sum(sites_over_ET)>0 % if at least 1 site is above ET
                target_outbreaks = vertcat(target_outbreaks, target_ID(rfs)); % add this reef ID to the list
            end
        end

        %then, find all reefs that are within 0.5° latitude N-S from outbreak target reefs
        close_reefs=[]; % let's gather reefs within 0.5° lat around each target % YM 09/25: quite large window, can give >3,400 reefs

        if ~isempty(target_outbreaks)

            for pob = 1:size(target_outbreaks,1)
                this_lat = META.reef_lat(find(META.reef_ID == target_outbreaks(pob)));
                % Calculate the absolute latitude difference with all reefs
                abs_lat_diff = abs(this_lat - META.reef_lat);
                % Create a logical index for reefs within 0.5 degrees
                within_range = abs_lat_diff < 0.5;
                % Use logical indexing to exclude those not within range
                this_np = find(within_range);
                if ~isempty(this_np)
                    close_reefs = [close_reefs; META.reef_ID(this_np)];
                end
            end
        end

        close_reefs = unique(close_reefs);%clean it up
        close_reefs = close_reefs(~ismember(close_reefs, target_ID)); %remove the reefs that are already in target list
        other_reefs = META.reef_ID(~ismember(META.reef_ID, vertcat(target_ID, close_reefs))); % find the remaining reefs
        full_list_ID = vertcat(target_ID, close_reefs, other_reefs); % update priority list with length = length(META.reef_ID)

    case {10, 11}  % Outbreak front: look for the AIMS sector with the highest CoTS density
        % 10: search for all reefs within a sector
        % 11: search only among priority reefs within a sector
        SectorIndices = unique(META.COTS_cull_reeflist.AIMS_sector, 'sorted'); % Number of AIMS sectors represented in the reefs in simulation
        sorted_reef_ID = cell(length(SectorIndices), 1); % Initialize a cell array to store the sorted indices for each sector

        % First calculate the overall COTS density for each sector
        sectorCOTSDensity = zeros(length(SectorIndices), 1);

        for i = 1:length(SectorIndices)
            sector = SectorIndices(i);
            reef_ID_sector = META.COTS_cull_reeflist.Reef_ID(find(META.COTS_cull_reeflist.AIMS_sector == sector));
            sectorCOTSDensity(sector) = sum(current_COTS_densities(reef_ID_sector,META.COTS_adult_min_age:end), 'all');
        end
        [~, sortedSectorIndices] = sort(sectorCOTSDensity, 'descend'); % Sort sectors based on overall COTS density

        % Then loop through each sorted sector and order reefs
        switch META.COTS_reefs2cull_strat
            case 10 ; ReefList = vertcat(target_ID, priority_ID);
            case 11 ; ReefList = vertcat(target_ID, priority_ID); % search only within the priority list
        end

        % Then loop through each sorted sector and order priority reefs
        for i = 1:length(sortedSectorIndices)
            sector = sortedSectorIndices(i);
            reef_ID_sector = intersect(META.COTS_cull_reeflist.Reef_ID(find(META.COTS_cull_reeflist.AIMS_sector == sector)), ReefList); % Find indices of reefs in the current sector that match the list
            [~, sortedTimestepIndices] = sort(sum(current_COTS_densities(reef_ID_sector,META.COTS_adult_min_age:end),2), 'descend'); % Sort the indices based on COTS densities for each timestep
            sorted_reef_ID{i} = reef_ID_sector(sortedTimestepIndices);  % Save the sorted indices for the current sector
        end

        full_list_ID = cat(1, sorted_reef_ID{:});     % Combine the sorted indices for all sectors into a column vector

    case {12, 13} % Effort sink: at each timestep, still make priority reef list as before but don't include reefs with lots of COTS
        % 12: no more than 3 CoTS per tow (META.max_COTS).
        % 13: no more than 3 CoTS per tow (META.max_COTS), otherwise with total coral cover more than 20% (META.min_control_cover)
        priority_list_tmp0 = vertcat(target_ID, priority_ID, nonpriority_ID); % each list might be shuffled at every time step, but doesn't matter here.

        COTS_per_tow = (0.22/0.6)*sum(current_COTS_densities(priority_list_tmp0,META.COTS_adult_min_age:end),2); % Convert into COTS per tow and sort in priority order

        switch META.COTS_reefs2cull_strat
            case 12
                I = find(COTS_per_tow <= META.max_COTS);
            case 13
                p_coral_cover = total_coral_pct2D(priority_list_tmp0); % total coral cover at current timestep
                I = find((COTS_per_tow <= META.max_COTS) | (COTS_per_tow > META.max_COTS & p_coral_cover > META.min_control_cover));
        end
        priority_list_tmp1 = priority_list_tmp0(I); % only keeps reef ID with more than 3 CoTS per tow

        J = priority_list_tmp0(~ismember(priority_list_tmp0, priority_list_tmp1)); % finds remaining reef ID
        full_list_ID = vertcat(priority_list_tmp1, J); % add the list of non priority reefs underneath (preserves first prioritisation)

    case {14, 15}  % Protection status - At each timestep, still make priority reef list etc but only include green zone reefs.
        priority_list_tmp0 = vertcat(target_ID, priority_ID, nonpriority_ID); % each list might be shuffled at every time step, but doesn't matter here.
        priority_list_GreenZone = (META.COTS_cull_reeflist.GreenZone(priority_list_tmp0)==1);

        switch META.COTS_reefs2cull_strat
            case 14; priority_list_tmp1 = priority_list_tmp0(priority_list_GreenZone==1); % only keeps reef ID within Green Zones
            case 15; priority_list_tmp1 = priority_list_tmp0(priority_list_GreenZone==0); % only keeps reef ID within Blue Zones
        end

        J = priority_list_tmp0(~ismember(priority_list_tmp0, priority_list_tmp1)); % finds remaining reef ID
        full_list_ID = vertcat(priority_list_tmp1, J); % add the list of non priority reefs underneath (preserves first prioritisation)

    case { 16, 17, 18}  % Weighting with CoTS connec (potential as source for CoTS larvae) and coral cover.
        % Where similar or in same range, weight to green or blue zone.
        % #16: priority for green zones; #17: priority for blue zones; #18: no preferential weighting relative to zoning
        priority_list_tmp0 = vertcat(target_ID, priority_ID); % each list might be shuffled at every time step, but doesn't matter here.
        nonpriority_list_tmp0 = nonpriority_ID; % each list might be shuffled at every time step, but doesn't matter here.

        % First, extract total coral cover
        p_coral_cover = total_coral_pct2D(priority_list_tmp0); % Extract total coral cover at current timestep and sort by prioritisation
        np_coral_cover = total_coral_pct2D(nonpriority_list_tmp0); % Extract total coral cover at current timestep and sort by prioritisation

        % Now get larval output at the time step before. Larvae only in summer, so summer and winter values were
        % combined here
        p_larval_output = COTS_larval_output(priority_list_tmp0);
        np_larval_output = COTS_larval_output(nonpriority_list_tmp0);

        % Combine then normalise to get a score
        p_score = zscore((p_coral_cover+1).*(p_larval_output+1)); % adding 1 to avoid 0 as minimum
        np_score = zscore((np_coral_cover+1).*(np_larval_output+1)); % adding 1 to avoid 0 as minimum

        % Now sort the list based on the combined score
        [~, p_sorted_indices] = sort(p_score, 'descend');
        priority_list_tmp1 = priority_list_tmp0(p_sorted_indices);

        [~, np_sorted_indices] = sort(np_score, 'descend');
        nonpriority_list_tmp1 = nonpriority_list_tmp0(np_sorted_indices);

        % Prioritise following zoning status
        switch META.COTS_reefs2cull_strat
            case 16  % Prioritize reefs with META.GreenZone == 1 in case of ties
                I = find(META.COTS_cull_reeflist.GreenZone(priority_list_tmp1)==1);
                J = find(META.COTS_cull_reeflist.GreenZone(nonpriority_list_tmp1)==1);

            case 17
                I = find(META.COTS_cull_reeflist.GreenZone(priority_list_tmp1)==0);
                J = find(META.COTS_cull_reeflist.GreenZone(nonpriority_list_tmp1)==0);

            case 18
                I = []; J = [];
        end

        priority_list_tmp2 = priority_list_tmp1(I); % New top of priority list
        nonpriority_list_tmp2 = nonpriority_list_tmp1(J); % New top of nonpriority list

        % Combine all lists - this preserves the order priority > nonpriority while re-ordering for zoning within each
        full_list_ID = vertcat(priority_list_tmp2, priority_list_tmp1(~ismember(priority_list_tmp1, priority_list_tmp2)),...
            nonpriority_list_tmp2, nonpriority_list_tmp1(~ismember(nonpriority_list_tmp1, nonpriority_list_tmp2))); % New list
    case 19 % control - maximise benefits from Tina's NESP results - WIP
         full_list_ID = vertcat(target_ID, priority_ID);

    case { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 } % Suki March 2026: Adaptive control scenarios and alternatives
        % 20 = default adaptive scenario (cat 1 cull, cat 2 +monitor, cat 3 +culling, cat 4 relocate 250km, cat 5 relocate 1000km)
        % 21 = shifted +1 cat (cat 1-2 cull, cat 3 +monitor, cat 4 +culling, cat 5 relocate 1000km)
        % 22 = shifted -1 cat (cat 1 +monitor, cat 2 +culling, cat 3-4 relocate 250km, cat 5 relocate 1000km)
        % 23 = alternative regional bleaching threshold: skip region when >=60% of reefs are severely bleached
        % 24 = alternative regional bleaching threshold: skip region when >=30% of reefs are severely bleached
        % 25 = alternative remaining workable region threshold: expand effort to non-priority reefs even if there is 1 available workable region (i.e., all_skip requires 4 regions all to be skipped, here we require only 3 regions. if there are at least 3 skip regions). e.g., FN, C, and S are skip regions, in these regions, we expand effort to cull non-priority reefs even though N is workable.
        % 26 = alternative remaining workable region threshold: do not redistribute effort to closest workable region, always expand effort. effort always stay in same region with no travelling.
        % 27 = alternative knowledge level: CoTS risk + predicted manta
        % 28 = alternative knowledge level: CoTS benefits + predicted manta
        % 29 = adaptive - Maximise benefits

        % Case-specific bleaching thresholds for reef-level decisions (cases 21 and 22 shift the default case-20 thresholds by +/-1 category)
        % Cases 23 and 24 use the same reef-level thresholds as case 20 but change the regional bleaching threshold.
        switch META.COTS_reefs2cull_strat
            case 20; cull_threshold = 3; cand_threshold = 3; radius_by_cat = [0 0 0 250 1000]; regional_threshold = 0.9; % cat 4->250km, cat 5->1000km; skip region when >=90% severely bleached
            case 21; cull_threshold = 4; cand_threshold = 4; radius_by_cat = [0 0 0 0 1000];   regional_threshold = 0.9; % cat 5->1000km only
            case 22; cull_threshold = 2; cand_threshold = 2; radius_by_cat = [0 0 250 250 1000]; regional_threshold = 0.9; % cat 3-4->250km, cat 5->1000km
            case 23; cull_threshold = 3; cand_threshold = 3; radius_by_cat = [0 0 0 250 1000]; regional_threshold = 0.6; % same as case 20 but skip region when >=60% severely bleached
            case 24; cull_threshold = 3; cand_threshold = 3; radius_by_cat = [0 0 0 250 1000]; regional_threshold = 0.3; % same as case 20 but skip region when >=30% severely bleached
            otherwise; cull_threshold = 3; cand_threshold = 3; radius_by_cat = [0 0 0 250 1000]; regional_threshold = 0.9; % default for cases 25-29, same as case 20
        end

        % Filter all reef lists to workable reefs only: within GBRMP, within 250 km of port, and shallow
        workable_reefIDs = META.workable.ReefID(META.workable.GBRMP == 1 & META.workable.Within250kmOfPort == 1 & META.workable.Shallow == 1);

        priority_list_tmp0 = vertcat(target_ID, priority_ID);  % each list might be shuffled at every time step, but doesn't matter here.
        priority_list_tmp0 = priority_list_tmp0(ismember(priority_list_tmp0, workable_reefIDs)); % only keep workable reefs
        % default: only select the 500 priority reefs

        % if doing maximise benefit scenarios - sort reef list first
        %switch
        %    case 29
        %        priority_list_tmp0 = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, priority_list_tmp0, control_records_prev, current_COTS_per_tow, 4, 0); % refugia + BZ
        %end

        %% Regional bleaching check point
        % Map reef IDs in reeflist to position indices, then assign bleaching info
        [~, reef_pos] = ismember(META.COTS_cull_reeflist.Reef_ID, META.reef_ID);
        META.COTS_cull_reeflist.bleaching_category = bleaching_category(reef_pos);
        META.COTS_cull_reeflist.severely_bleached = severely_bleached(reef_pos);

        % sort non-priority by CoTS risk and bleaching category
        % expand to non-priority reefs if necessary
        nonpriority_tmp0 = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, nonpriority_ID(ismember(nonpriority_ID, workable_reefIDs)), control_records_prev, current_COTS_per_tow, 6, 0);

        % Calculate proportion of reefs by region with severely bleached reef = 1
        regional_bleaching = zeros(length(META.COTS_cull_region), 1);
        
        for r = 1:length(META.COTS_cull_region)
            current_region = META.COTS_cull_region{r};
            region_reefs = META.COTS_cull_reeflist(ismember(META.COTS_cull_reeflist.Region, current_region), :);

            if height(region_reefs) > 0
                num_severely_bleached = sum(region_reefs.severely_bleached == 1);
                regional_bleaching(r,1) = num_severely_bleached / height(region_reefs);
            end
        end
        
        % Identify regions exceeding the regional bleaching threshold (proportion of severely bleached reefs)
        skip_region_idx = find(regional_bleaching >= regional_threshold); % numeric indices into META.COTS_cull_region
        skip_region = META.COTS_cull_region(skip_region_idx); % region labels (e.g., 'FN', 'N', 'C', 'S')

        % Identify reefs in skip regions - they will be monitored only (not culled)
        skip_priority_reefs = [];
        for s = 1:length(skip_region)
            in_skip = META.COTS_cull_reeflist.Reef_ID(ismember(META.COTS_cull_reeflist.Region, skip_region{s}));
            skip_priority_reefs = vertcat(skip_priority_reefs, priority_list_tmp0(ismember(priority_list_tmp0, in_skip)));
        end

        % Non-skip-region priority reefs go through the reef-level bleaching check as before
        priority_list_tmp1 = priority_list_tmp0(~ismember(priority_list_tmp0, skip_priority_reefs));

        %% Reef level check point

        priority_list_tmp2 = []; % initialize the list of priority reefs to keep after bleaching check

        for reef = 1:length(priority_list_tmp1) % for every reef within the priority list

            % Skip if this reef is already in priority_list_tmp2 (e.g., added as a candidate from a previous iteration)
            if ismember(priority_list_tmp1(reef), priority_list_tmp2)
                continue
            end

            I_reef = find(META.reef_ID == priority_list_tmp1(reef)); % position index for this reef
            if bleaching_category(I_reef) <= cull_threshold
                priority_list_tmp2 = vertcat(priority_list_tmp2, priority_list_tmp1(reef)); % keep reef in list - further actions are set in f_COTS_control_NEW

            else         
                % search for nearby reefs with lower bleaching category to relocate efforts
                search_radius = radius_by_cat(bleaching_category(I_reef));

                distances = META.distance_matrix(:, priority_list_tmp1(reef)); % extract distances between this reef and all other reef
                reefs_within_radius = find(distances <= search_radius & distances > 0); % look for reefs within search radius that is not itself

                % Ensure intersection with valid reef IDs
                reefs_within_radius = intersect(reefs_within_radius, META.reef_ID);

                if isempty(reefs_within_radius)
                 continue % skip this reef - no replacement is found
                end

                [~, rwr_pos] = ismember(reefs_within_radius, META.reef_ID); % position indices for reefs within radius
                candidates = reefs_within_radius(bleaching_category(rwr_pos) <= cand_threshold); % look for reefs within radius that have bleaching category <= cand_threshold
                % Reef is only workable when it is within the GBRMP, within 250 km of a port, and shallow (shallows point of the reef <10 m depth).
                candidates = candidates(ismember(candidates, workable_reefIDs)); % only keep candidates that are workable reefs

                % Exclude reefs from skip regions as relocation candidates (they will be monitored, not culled)
                if ~isempty(skip_region)
                    skip_region_all_reefs = META.COTS_cull_reeflist.Reef_ID(ismember(META.COTS_cull_reeflist.Region, skip_region));
                    candidates = candidates(~ismember(candidates, skip_region_all_reefs));
                end

                [~, cand_pos] = ismember(candidates, META.reef_ID); % position indices for candidates
                if ~isempty(control_records_prev)
                    % Filter candidates: keep borderline-category reefs only if previously culled
                    candidates = candidates((bleaching_category(cand_pos) ~= cand_threshold) | ((bleaching_category(cand_pos) == cand_threshold) & ismember(candidates, control_records_prev.culled_reef_ID)));
                else
                    candidates = candidates(bleaching_category(cand_pos) < cand_threshold);
                end
                
                % if candidate already in the priority list, remove it from the list to avoid duplicates
                candidates = candidates(~ismember(candidates, priority_list_tmp2));

                if isempty(candidates)
                    continue % no valid candidate found after filtering, skip this reef
                end

                % Create scoring system to select the best candidate reef,
                % see f_scoring - Suki March 2026
                % Normalise to between 0 and 1, then weight by importance of each criteria (to be defined) - See Tina's code. 
                % - Manta tow COTS density (lower is better)
                % - Culled COTS density (lower is better) or benefits of culling (higher is better)
                % - Reef priority status (T > P > N)
                % - bleaching category (priority to 1, then 2, then 3)
                % - distance to the priority reef (closer is better)
                switch META.COTS_reefs2cull_strat
                    case 27 % score based on CoTS risk (density) + predicted manta tow density 
                        sorted_candidates = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, candidates, control_records_prev, current_COTS_per_tow, 2, priority_list_tmp1(reef));
                    case 28 % score based on CoTS benefits (culled density) + predicted manta tow density
                        sorted_candidates = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, candidates, control_records_prev, current_COTS_per_tow, 3, priority_list_tmp1(reef));
                    otherwise % default
                        sorted_candidates = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, candidates, control_records_prev, current_COTS_per_tow, 1, priority_list_tmp1(reef));
                end

                priority_list_tmp2 = vertcat(priority_list_tmp2, sorted_candidates(1)); % add highest-scored candidate reef to the list

                % Record relocation distance: original reef (bleaching cat>3) -> selected candidate
                relocation_dist_km.reef_ID         = [relocation_dist_km.reef_ID;         sorted_candidates(1)];
                relocation_dist_km.original_reef_ID = [relocation_dist_km.original_reef_ID; priority_list_tmp1(reef)];
                relocation_dist_km.dist_km          = [relocation_dist_km.dist_km;          META.distance_matrix(sorted_candidates(1), priority_list_tmp1(reef))];

            end
        end

        % Remove any skip-region priority reefs added as relocation candidates to avoid duplication
        skip_priority_reefs = skip_priority_reefs(~ismember(skip_priority_reefs, priority_list_tmp2));

        % Non-priority reefs to cull: filter to bleaching category <= cull_threshold and not already in priority list (tmp2)
        [~, np_pos] = ismember(nonpriority_tmp0, META.reef_ID);
        nonpriority_tmp1 = nonpriority_tmp0(bleaching_category(np_pos) <= cull_threshold & ~ismember(nonpriority_tmp0, priority_list_tmp2));

        % Build monitor_only_IDs: skip-region priority reefs are always monitor-only.
        % If the 'expand' condition is met: also include non-priority reefs with bleaching category > 3
        % (they are surveyed but cannot be culled efficiently), and no effort is redistributed to other regions.
        % Case 20 (default): expand only when ALL regions are skip.
        % Case 25: expand when at least 3 regions are skip (>=3 unworkable -> no redistribution).
        % Case 26: always expand (never redistribute, regardless of how many regions are skip).
        switch META.COTS_reefs2cull_strat
            case 25
                all_skip = length(skip_region_idx) >= 3;
            case 26
                all_skip = ~isempty(skip_region_idx);
            otherwise
                all_skip = length(skip_region_idx) == length(META.COTS_cull_region);
        end
        monitor_only_IDs = skip_priority_reefs;
        if all_skip
            nonpriority_cat_gt3 = nonpriority_tmp0(bleaching_category(np_pos) > cull_threshold & ~ismember(nonpriority_tmp0, priority_list_tmp2));
            monitor_only_IDs = vertcat(monitor_only_IDs, nonpriority_cat_gt3);
        end

        % Build full reef list: skip-region priority (monitor-only) first to protect their global budget,
        % then workable priority reefs (including relocation candidates), then cullable non-priority,
        % then (if all-skip) monitor-only non-priority reefs.
        % skip_priority_reefs must come first: they are cheap to monitor (1 dive/site) but could be
        % starved of global dives if relocation candidates (which may be large non-priority reefs) go first.
        full_list_ID = vertcat(skip_priority_reefs, priority_list_tmp2, nonpriority_tmp1);
        if all_skip
            full_list_ID = vertcat(full_list_ID, nonpriority_cat_gt3);
        end

        % unvisited: priority reefs skipped due to individual reef bleaching (cat 4/5, locally relocated)
        % Note: skip-region priority reefs are no longer unvisited - they are monitored (in monitor_only_IDs)
        skip_reef_local = priority_list_tmp1(~ismember(priority_list_tmp1, priority_list_tmp2));
        [~, idx_local] = ismember(skip_reef_local, META.COTS_cull_reeflist.Reef_ID);
        [~, region_local] = ismember(META.COTS_cull_reeflist.Region(idx_local), META.COTS_cull_region);

        if ~isempty(skip_reef_local)
            unvisited = table(skip_reef_local, region_local, ...
                'VariableNames', {'Reef_ID', 'Region'});
        else
            unvisited = table(zeros(0,1), zeros(0,1), ...
                'VariableNames', {'Reef_ID', 'Region'});
        end
        % Columns: Reef_ID, Reason (1=individual reef bleaching/local relocation), Region

end

%% Force to start with the last controlled reef, since it may not have been completely culled
if last_reef_COTScontrolled ~= 0 % would be 0 if no reefs have been visited before (ie, before control starts)
    New_toplist = full_list_ID(full_list_ID == last_reef_COTScontrolled,:);
    full_list_ID(full_list_ID == last_reef_COTScontrolled,:)=[]; % delete the last controlled reef from the list, wherever it is
    full_list_ID = [New_toplist ; full_list_ID]; % Put the last controlled reef on top of the list -> will be visited first
end