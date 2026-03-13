% ----------------------------------------------------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, Mar 2025.
%
% Creates a list of ordered reefs for prioritising CoTS culling under different scenarios of control management.
% Optimised and revised version of the previous scripts:
% - f_makeReefList from Karlo Hock (REEFMOD-GBR.6.3, Sep 2019)
% - f_makeReefListCCS from Carolina Castro-Sanguino (REEFMOD-GBR.6.6, Mar 2022)
% - f_makeReefListTS from Tina Skinner (REEFMOD-GBR.6.8, May 2023)
% ----------------------------------------------------------------------------------------------------------------------
function [full_list_ID, unvisited] = f_makeReefList_NEW(META, current_COTS_densities, current_reef_ET, COTS_densities_per_site,...
    total_coral_pct2D, COTS_larval_output, last_reef_COTScontrolled, all_DHWs)

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

% Initialize unvisited as empty - only populated for case 19
unvisited = [];

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
        close_reefs = setdiff(close_reefs, target_ID); %remove the reefs that are already in target list
        other_reefs = setdiff(META.reef_ID,vertcat(target_ID, close_reefs)); % find the remaining reefs
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
            case 10 ; ReefList = vertcat(target_ID, priority_ID, nonpriority_ID);
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

        J = setdiff(priority_list_tmp0,priority_list_tmp1); % finds remaining reef ID
        full_list_ID = vertcat(priority_list_tmp1, J); % add the list of non priority reefs underneath (preserves first prioritisation)

    case {14, 15}  % Protection status - At each timestep, still make priority reef list etc but only include green zone reefs.
        priority_list_tmp0 = vertcat(target_ID, priority_ID, nonpriority_ID); % each list might be shuffled at every time step, but doesn't matter here.
        priority_list_GreenZone = (META.COTS_cull_reeflist.GreenZone(priority_list_tmp0)==1);

        switch META.COTS_reefs2cull_strat
            case 14; priority_list_tmp1 = priority_list_tmp0(priority_list_GreenZone==1); % only keeps reef ID within Green Zones
            case 15; priority_list_tmp1 = priority_list_tmp0(priority_list_GreenZone==0); % only keeps reef ID within Blue Zones
        end

        J = setdiff(priority_list_tmp0,priority_list_tmp1); % finds remaining reef ID
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
        full_list_ID = vertcat(priority_list_tmp2, setdiff(priority_list_tmp2, priority_list_tmp1),...
            nonpriority_list_tmp2, setdiff(nonpriority_list_tmp2, nonpriority_list_tmp1)); % New list

    case 19 % relocation only Suki Feb 2026
         priority_list_tmp0 = vertcat(target_ID, priority_ID);  % each list might be shuffled at every time step, but doesn't matter here. 
         % Here we assume priority >>> non-priority. so even though a
         % p reef exceeds the DHW threshold, it still has a higher
         % priority than a n reef. 

         % sort visitation order based on distance from port
         p_distance_score = META.distance_port(priority_list_tmp0);
         [~, p_distance_sorted_indices] = sort(p_distance_score, 'ascend'); 
         priority_list_tmp1 = priority_list_tmp0(p_distance_sorted_indices);

         n_distance_score = META.distance_port(nonpriority_ID);
         [~, n_distance_sorted_indices] = sort(n_distance_score, 'ascend'); 
         nonpriority_list_tmp0 = nonpriority_ID(n_distance_sorted_indices);

         % set DHW threshold and search radius
         DHW_threshold = 2;
         search_radius = 50;
         
         % create all_DHWs table with reef ID for indexing
         all_DHWs_ID = zeros(META.nb_reefs,2);
         all_DHWs_ID(:,1) = META.reef_ID;
         all_DHWs_ID(:,2) = all_DHWs(META.reef_ID)';

         % isolate reefs not exceeding DHW threshold
         non_exceeding_reefs = all_DHWs_ID(all_DHWs_ID(:,2) <= DHW_threshold); 

         priority_list_tmp2 = [];

         for reef = 1:length(priority_list_tmp1) % for every reef within the priority list

             priority_reef_id = priority_list_tmp1(reef);

             if all_DHWs(priority_reef_id) <= DHW_threshold % if this reef is below DHW threshold
                 priority_list_tmp2(end+1) = priority_reef_id; % keep the reef in the list
                 continue
             end
             
             % otherwise find replacement (relocate efforts)

             % extract distances between this reef and all other reef
             distances = META.distance_matrix(:, priority_list_tmp1(reef)); 

             reefs_within_radius = find(distances <= search_radius & distances > 0); % look for reefs within search radius that is not itself
                 
             if isempty(reefs_within_radius)
                 continue % skip this reef - no replacement is found
             end

             % sort by distance (ascending)              
             [~, sort_idx] = sort(distances(reefs_within_radius));
             reefs_within_radius = reefs_within_radius(sort_idx);

             % apply criteria for suitable reefs = within radius, below DHW threshold, is priority reef
             suitable_reefs = reefs_within_radius(ismember(reefs_within_radius, non_exceeding_reefs) & ...
                 ismember(reefs_within_radius, priority_list_tmp1));

             if isempty(suitable_reefs)
                 continue % skip this reef - no replacement is found
             end

             % add closest suitable reef instead
             priority_list_tmp2(end+1) = suitable_reefs(1);
         end

         priority_list_tmp2 = unique(priority_list_tmp2,'stable')';
         unvisited = setdiff(priority_list_tmp1,priority_list_tmp2);

         % Combine all lists - this preserves the order priority > nonpriority
         % in this case priority + DHW above threshold is still better than
         % any non-priority reef
        full_list_ID = vertcat(priority_list_tmp2, unvisited, nonpriority_list_tmp0);

end

%% Force to start with the last controlled reef, since it may not have been completely culled
if last_reef_COTScontrolled ~= 0 % would be 0 if no reefs have been visited before (ie, before control starts)
    New_toplist = full_list_ID(full_list_ID == last_reef_COTScontrolled,:);
    full_list_ID(full_list_ID == last_reef_COTScontrolled,:)=[]; % delete the last controlled reef from the list, wherever it is
    full_list_ID = [New_toplist ; full_list_ID]; % Put the last controlled reef on top of the list -> will be visited first
end