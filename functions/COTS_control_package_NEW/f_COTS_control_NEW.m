% ----------------------------------------------------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, Mar 2025.
%
% Implements CSIRO-like control effort and triggers under various scenarios of control management.
%
% Revised, corrected and optimised from code originally written by Karlo Hock:
% - f_COTS_control_K (REEFMOD-GBR.6.3, Sep 2019)
% - f_COTS_control_CSIRO (REEFMOD-GBR.6.6, Sep 2021)
% ----------------------------------------------------------------------------------------------------------------------

function [COTS_all_densities, control_records, last_reef_COTScontrolled] = ...
    f_COTS_control_NEW(META, COTS_all_densities, last_reef_COTScontrolled, total_coral_pct2D, COTS_larval_output, CONNECT_CORAL, t, ...
    all_DHWs, bleaching_category, severely_bleached, control_records_prev)

%% TEMP FOR DEBUGGING (see also last line for recording last_reef_COTScontrolled in RESULT)
% function [RESULT, control_records, last_reef_COTScontrolled] = f_COTS_control_NEW(META, RESULT, t)
% 
% if t==META.COTS_control_start%set this to zero if this is the first time COTS control is run
%     RESULT.last_reef_COTScontrolled=0;%no reefs controlled before, so set this to zero for makeReefList file
% end
% t=t+1;
% COTS_all_densities = reshape(RESULT.COTS_all_densities(:,t,1:META.COTS_maximum_age),META.nb_reefs,META.COTS_maximum_age);
% % COTS_all_densities = reshape(RESULT,META.nb_reefs,META.COTS_maximum_age);

%%

% COTS_all_densities: current CoTS densities for all classes at the end of the time step - will be updated at the end
% (after culling). CoTS densities right before culling will be recorded in 'COTS_records'
% last_reef_COTScontrolled is zero before CoTS control starts, ie, no reefs controlled before (used by makeReefList)
current_COTS_densities = squeeze(COTS_all_densities); % density of COTS at every age for every reef

% Sum across subadults and adults after adjusting for imperfect
% detectability - YM May 2026
current_COTS_total_densities = sum(current_COTS_densities(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
% Estimate the equivalent CoTS per tow
current_COTS_per_tow = f_convert_CoTS_density_2_tow(current_COTS_total_densities, META.total_area_cm2);

% thisboatorder = META.boatProperties.fixedOrder(1); % each boat to visit reefs in fixed (1) or random (0) order but only considering first boat for now
timestep_dives = sum(META.boatProperties.totalInidvDives); % full dive quota for timestep
if META.use_regional_dives == 1
    region_dives = round(META.regional_effort_allocation * timestep_dives); % dive quota per region for timestep
    % note there could be rounding residuals but the number is negligible 
end

% Travel cost conversion: 8 nm/hr steaming speed, 12 hr/day steaming available, 16 dive-hr lost per steaming day, 45 min per dive
% km -> nm -> steaming hours -> steaming days -> dive-hours lost -> dives lost
km_to_dives = @(d) ceil(d ./ 1.852 ./ 8 ./ 12 .* 16 ./ 0.75);

%% Set the ecological threshold (ET) for each reef
% - dynamically for a current level of coral cover; used in both reef-level and site-level calculations since coral tracked on a reef-level
% current_COTS_ET = f_calculate_COTS_ET(RESULT,t, META);
% - fixed value independent of coral cover as in Castro et al (2023). Deliberate choices based on a few reasons,
% including that GBRMPA did not always implement a dynamic ET during the control program. ET = 0.075 CoTS per tow (Fletcher et al. 2021)
current_reef_ET = repmat(0.075,[size(current_COTS_densities,1) 1]); % CANNOT BE CHANGED AS THE EQUATION OF EFFORT REQUIRED TO CULL TO ET ONLY WORKS FOR ET=0.075
% current_COTS_ET = repmat(0,[size(current_COTS,1) 1]); % JUST TO COMPARE WITH PREVIOUS RUNS BUT WRONG, ONLY WORKS FOR ET=0.075

%% Redistribute reef-level COTS density to individual control sites based on randomly generated proportions.
% Gives a cell array with, for every reef, a matrix of CoTS density for all size classes at every site of the reef.
% Also keeps track of generated proportions per site. YM 03/2025: corrected to avoid negative values.
[COTS_densities_per_site, stored_betarndN] = f_redistrib_COTS(current_COTS_densities, META.COTS_cull_reeflist.nb_sites);

%% Make a list of reefs that determines the order in which they are visited
% Only do this once per timestep, then go down the list checking the individual triggers until out of quota
% current_COTS_ET = current_reef_ET;
% current_COTS = current_COTS_densities;
% [~, criteriaS, global_trigger, RESULT] = f_makeReefList_TS(META, RESULT, t, 0, current_COTS, thisboatorder, current_COTS_ET, COTS_densities_per_site);
% ReefList = criteriaS.criteria(:,1); % TEMP - from now only use the ordered list of reef index
[ReefList, unvisited, monitor_only_IDs, relocation_dist_km, skip_region_idx_out, regional_bleaching_out] = f_makeReefList_NEW(META, current_COTS_densities, current_reef_ET, COTS_densities_per_site, total_coral_pct2D, COTS_larval_output, ...
    CONNECT_CORAL, t, last_reef_COTScontrolled, ...
    bleaching_category, severely_bleached, control_records_prev, ...
    current_COTS_per_tow);

% ReefList must be a selection of META.reef_ID sorted in decreasing order of priority for control
% Keep track record of key control variables
control_records=struct('culled_reef_ID',[], 'culled_density_reef',[], 'nb_culled_COTS',[],...
    'nb_dives',[], 'nb_culled_sites',[], 'nb_culled_reefs',[], 'nb_visited_reefs',[], ...
    'nb_control_sites' , [], 'nb_unvisited_reefs', height(unvisited), ...
    'additional_cull_ID', [], 'nb_additional_cull', 0, ...
    'skip_region_monitor_ID', [], 'skip_region_monitor_dives', [], ...  % reefs monitored due to regional bleaching trigger
    'skip_local_monitor_ID', [], 'skip_local_monitor_dives', [], ...    % priority reefs monitored because individually too hot (no valid relocation candidate)
    'additional_monitor_ID', [], 'additional_monitor_dives', [], ...
    'total_monitor_dives', 0, ...
    'nb_monitored_reefs', 0, 'remaining_dives', 0, 'remaining_region_dives', [], ...
    'visited_reef_ID', [], 'visited_manta_tow', [], ...
    'total_expended_dives', 0, ...
    'region_break_reef_ID', [], 'region_break_region', [], 'region_break_remaining', [], ...
    'global_break_reef_ID', [], 'global_break_remaining', [], ...
    'regional_effort_allocation', [], ... % [FN N C S] proportions used this timestep - for tracking and debugging
    'site_skip_count', 0, 'site_skip_reef_ID', [], ... % count of skip-then-find events per timestep
    'site_skip_first_site', [], 'site_skip_success_site', [], ... % first site that exceeded budget, and the site that fit within budget
    'port_travel_dist_km', [], 'port_travel_dives', [], 'port_travel_from_region', [], 'port_travel_to_region', [], ... % case 1: inter-port travel when redistributing effort across regions
    'reef_travel_dist_km', [], 'reef_travel_dives', [], 'reef_travel_from', [], 'reef_travel_to', [], 'total_travel_dives', 0, ... % case 2: inter-reef travel between relocation reefs and original bleaching cat>3 reefs
    'skipped_regions', skip_region_idx_out, ... % numeric indices into META.COTS_cull_region of regions that exceeded the bleaching threshold this timestep ([] for non-adaptive strategies)
    'regional_bleaching_pct', regional_bleaching_out); % proportion of severely bleached reefs per region at this timestep ([] for non-adaptive strategies)
% Suki Apr 2026: site_skip_* fields:
% These track the situation where a cull site on a reef requires more dives
% than the remaining budget (global or regional), so it is skipped, and a
% cheaper site on the SAME reef is found and successfully culled instead.
%   site_skip_count:        total number of times this happened in the timestep
%   site_skip_reef_ID:      reef where each event occurred (vector, one entry per event)
%   site_skip_first_site:   the site index that was first skipped (too expensive)
%   site_skip_success_site: the site index that was culled instead (affordable)
% Multiple events can occur on the same reef if several expensive sites are
% skipped before each cheaper site is found. If all remaining sites on a
% reef are too expensive (no cheaper site found), NO event is recorded.

visited_reefs = 0; % counter to keep track of how many reefs were culled
n = 1; % Start with the first reef on the list
remaining_dives = timestep_dives; %full dive quota at the beginning of a timestep
if META.use_regional_dives == 1
    remaining_region_dives = region_dives; % full regional dive quota at the beginning of a timestep
else
    % Global-only mode: set regional budgets to inf so every regional check
    % (remaining_region_dives < X) evaluates to false throughout the loop,
    % effectively disabling all regional constraints without guarding each site.
    remaining_region_dives = inf(length(META.COTS_cull_region), 1);
end

% Bleaching response thresholds for adaptive strategies (20 = default, 21 = shift +1 cat, 22 = shift -1 cat)
is_adaptive = ismember(META.COTS_reefs2cull_strat, [20 21 22 23 24 25 26 27 28 29]); % cases 23-29 use same reef level thresholds as 20
if is_adaptive
    switch META.COTS_reefs2cull_strat
        case 20; add_monitor_cat = 2; add_cull_cat = 3;
        case 21; add_monitor_cat = 3; add_cull_cat = 4;
        case 22; add_monitor_cat = 1; add_cull_cat = 2;
        otherwise; add_monitor_cat = 2; add_cull_cat = 3; % default for cases 23-26, same as case 20
    end
end

%% Suki Apr 2026 - regional effort redistribution for skip regions (strategy 19)
% Skip-region priority reefs are monitored (not culled) - handled using monitor_only_IDs.
% After accounting for expected monitoring cost, any surplus regional budget is redistributed
% to the closest workable (non-skip) adjacent region following proximity rules:
%   FN(1) -> N(2) only;  N(2) -> C(3), S(4), FN(1);  C(3) -> S(4), N(2) only;  S(4) -> C(3), N(2) only
% If no workable adjacent region exists (e.g., FN when N is also a skip region), the surplus
% stays in that region and is consumed by non-priority culling/monitoring within it.
if META.use_regional_dives && is_adaptive && ~isempty(monitor_only_IDs)
    proximity_order = {2; [3 4 1]; [4 2]; [3 2]}; 

    % Identify skip-region priority reefs (type T or P, AND in a skip region) within the monitor_only list.
    % Excludes individually-hot priority reefs (skip_reef_local) which are also T/P in monitor_only_IDs
    % but should NOT trigger effort redistribution - their regional budget stays in their own region.
    [~, mo_rows] = ismember(monitor_only_IDs, META.COTS_cull_reeflist.Reef_ID);
    mo_types = META.COTS_cull_reeflist.reef_type(mo_rows);
    [~, mo_region_tmp] = ismember(META.COTS_cull_reeflist.Region(mo_rows), META.COTS_cull_region);
    priority_monitor_IDs = monitor_only_IDs(~strcmp(mo_types, 'N') & ismember(mo_region_tmp, skip_region_idx_out)); % T or P in skip regions only

    if ~isempty(priority_monitor_IDs)
        [~, pm_rows] = ismember(priority_monitor_IDs, META.COTS_cull_reeflist.Reef_ID);
        [~, pm_region_idx] = ismember(META.COTS_cull_reeflist.Region(pm_rows), META.COTS_cull_region);
        skip_regions_idx = unique(pm_region_idx);
        % Determine whether to expand effort within skip regions (no redistribution) or redistribute to workable regions.
        % Case 20 (default): expand only when ALL regions are skip.
        % Case 25: expand when at least 3 regions are skip.
        % Case 26: always expand (never redistribute).
        switch META.COTS_reefs2cull_strat
            case 25; all_skip = length(skip_regions_idx) >= 3;
            case 26; all_skip = ~isempty(skip_regions_idx);
            otherwise; all_skip = length(skip_regions_idx) == length(META.COTS_cull_region);
        end

        if ~all_skip
            % For each skip region: compute monitoring cost for its priority reefs,
            % then redistribute surplus to the closest workable adjacent region.
            for si = 1:length(skip_regions_idx)
                r = skip_regions_idx(si);
                reefs_this_skip = priority_monitor_IDs(pm_region_idx == r);
                [~, rs_rows] = ismember(reefs_this_skip, META.COTS_cull_reeflist.Reef_ID);
                monitoring_cost_r = sum(META.COTS_cull_reeflist.nb_sites(rs_rows)); % 1 monitoring dive per cull site
                surplus = max(0, remaining_region_dives(r) - monitoring_cost_r);

                if surplus > 0
                    adj_regions = proximity_order{r};
                    for c = 1:length(adj_regions)
                        if ~ismember(adj_regions(c), skip_regions_idx)
                            % Case 1: compute inter-port travel cost (skip region port -> receiving region port)
                            % before transferring surplus to the receiving region
                            port_row = META.port_distances.IN_REGION == r & META.port_distances.NEAR_REGION == adj_regions(c);
                            if ~any(port_row)
                                continue % no valid port route to this adjacent region - try the next one in proximity order
                            end
                            port_dist_km = META.port_distances.NEAR_DIST_KM(port_row);
                            port_travel_dv = km_to_dives(port_dist_km);
                            net_surplus = max(0, surplus - port_travel_dv); % surplus after travel cost
                            remaining_region_dives(adj_regions(c)) = remaining_region_dives(adj_regions(c)) + net_surplus;
                            remaining_region_dives(r) = remaining_region_dives(r) - surplus; % debit full surplus from skip region
                            remaining_dives = max(0, remaining_dives - port_travel_dv); % deduct travel cost from global pool
                            % Record inter-port travel
                            control_records.port_travel_dist_km    = [control_records.port_travel_dist_km;    port_dist_km];
                            control_records.port_travel_dives       = [control_records.port_travel_dives;      port_travel_dv];
                            control_records.port_travel_from_region = [control_records.port_travel_from_region; r];
                            control_records.port_travel_to_region   = [control_records.port_travel_to_region;  adj_regions(c)];
                            break
                        end
                    end
                    % If no workable adjacent region found, surplus stays for non-priority use in region r 
                end
            end
        end
        % If all regions are skip: no redistribution - all budgets used for monitoring + non-priority use within own region
    end
end

%% Check CoTS density of this reef

while remaining_dives > 0 && n <= length(ReefList) %while there are dives remaining and the list not exhausted
    % General description from KH: get first reef ID from the list, check if reef-level trigger is used, then use it;
    % check if site-level trigger is valid at any sites. Go only to those sites, recheck after every run to see if trigger
    % still valid, if not, increase the current_reef by one until out of quota
   
    % Suki Apr 2026: If all regions are exhausted, break the global loop
    if all(remaining_region_dives <= 0)
        remaining_dives = 0;
        break
    end

    % Check if reef level trigger is used; if yes, check whether this reef satisifies the criterion; if not, advance to the next reef on list
    if META.COTS_reef_trigger==1 % YET TO BE TESTED

        treat_this_reef = 0;

        while treat_this_reef==0

            this_reef_ID = ReefList(n);
            I = find(META.reef_ID == this_reef_ID); % locate this reef in the reef definition list

            if META.doing_COTS_heat_stress_detectibility == 1
                current_COTS_per_tow(I) = current_COTS_per_tow(I).*exp(-0.25 * all_DHWs(I)); % Apply detectability reduction due to heat stress if applicable (Cook et al. in review)
            end
            
            if current_COTS_per_tow(I) < 0.22 %% If number of CoTS per tow on the reef is above outbreak threshold (Moran and De'ath 1992)

                treat_this_reef = 1; % OK let's do the culling

            else % go to the next reef
                n = n + 1;

                if n > length(ReefList) % break the control while loop if all reefs have been visited
                    break
                end
            end
        end
    end

    % For the selected reef
    this_reef_ID = ReefList(n); % (needs to be done again if condition above wasn't met)
    I = find(META.reef_ID == this_reef_ID); % locate this reef in the reef definition list

    % Suki March 2026
    % Look up this reef's region and check regional dive budget
    this_reef_region = META.COTS_cull_reeflist.Region(META.COTS_cull_reeflist.Reef_ID == this_reef_ID);
    [~, region_idx] = ismember(this_reef_region, META.COTS_cull_region);

    if remaining_region_dives(region_idx) <= 0
        n = n + 1;
        continue % skip this reef - no dives remaining for its region
    end

    %% Monitor-only handling for adaptive strategies (20/21/22)
    % Covers: (a) skip-region priority reefs, (b) non-priority reefs above culling threshold when all regions are skip.
    % These reefs are surveyed (monitoring dives) but NOT culled.
    if is_adaptive && ismember(this_reef_ID, monitor_only_IDs)
        monitor_dives_mo = ceil(META.COTS_cull_reeflist.nb_sites(I) * 1); % 1 monitoring dive per cull site

        % Check global budget
        if remaining_dives < monitor_dives_mo
            n = n + 1;
            continue
        end

        % Check regional budget
        if remaining_region_dives(region_idx) < monitor_dives_mo
            n = n + 1;
            continue % no regional budget available for monitoring
        end

        % Budget confirmed: debit global and regional pools
        remaining_dives = remaining_dives - monitor_dives_mo;
        remaining_region_dives(region_idx) = remaining_region_dives(region_idx) - monitor_dives_mo;

        % Record monitoring - route to the correct tracking field:
        %   skip-region: reef's region exceeded the regional bleaching threshold
        %   local-skip:  reef was individually too hot with no valid relocation candidate
        if ismember(region_idx, skip_region_idx_out)
            control_records.skip_region_monitor_ID    = [control_records.skip_region_monitor_ID;    this_reef_ID];
            control_records.skip_region_monitor_dives = [control_records.skip_region_monitor_dives; monitor_dives_mo];
        else
            control_records.skip_local_monitor_ID    = [control_records.skip_local_monitor_ID;    this_reef_ID];
            control_records.skip_local_monitor_dives = [control_records.skip_local_monitor_dives; monitor_dives_mo];
        end
        control_records.total_monitor_dives = sum(control_records.skip_region_monitor_dives) + sum(control_records.skip_local_monitor_dives) + sum(control_records.additional_monitor_dives);
        control_records.nb_monitored_reefs  = control_records.nb_monitored_reefs + 1;
        control_records.nb_control_sites    = [control_records.nb_control_sites; META.COTS_cull_reeflist.nb_sites(I)];

        % Record as visited (manta tow)
        % Sum detected densities across age classes FIRST, then convert (correct: f(sum) not sum(f), nonlinear conversion)
        tmp_COTS_per_site = COTS_densities_per_site{I,1};
        tmp_detected_total = sum(tmp_COTS_per_site(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end), 2);
        tmp_COTS_per_tow   = f_convert_CoTS_density_2_tow(tmp_detected_total, META.total_area_cm2);
        if META.doing_COTS_heat_stress_detectibility == 1
            tmp_COTS_per_tow = tmp_COTS_per_tow .* exp(-0.25 * all_DHWs(I));
        end

        control_records.visited_reef_ID   = [control_records.visited_reef_ID; this_reef_ID];
        control_records.visited_manta_tow = [control_records.visited_manta_tow; sum(tmp_COTS_per_tow)];

        visited_reefs = visited_reefs + 1;
        n = n + 1;
        continue % monitoring done - skip culling for this reef
    end

    %% Extract CoTS density for this reef
    % Extract COTS density for all sites
    this_reef_COTS_densities_per_site = COTS_densities_per_site{I,1};
    
    % Keep per-age-class detected densities — needed later for the weighted culling efficiency under heat stress
    this_reef_COTS_densities_per_site_detected = this_reef_COTS_densities_per_site(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end);
    % Sum detected densities across ALL adult age classes first (correct order: f(sum) not sum(f), because the Moran & De'ath conversion is nonlinear)
    this_reef_detected_total_per_site = sum(this_reef_COTS_densities_per_site_detected, 2); % [n_sites x 1] aggregate detected density per site
    % Estimate the equivalent CoTS per tow — nonlinear conversion applied ONCE to the aggregate
    this_reef_COTS_per_tow_per_site = f_convert_CoTS_density_2_tow(this_reef_detected_total_per_site, META.total_area_cm2); % [n_sites x 1]


    %% Additional monitoring effort + travel costs for bleaching category 2 and relocation reefs (effort debited before culling)
    % Suki March 2026 - adaptive control during bleaching
    % At the reef level, if the bleaching is mild, culling operations
    % continue but with addition effort for monitoring
    additional_monitoring = 0; % default no additional monitoring

    if is_adaptive
        monitor_dives = 0; % initialise monitor dives

        if bleaching_category(I) == add_monitor_cat
            additional_monitoring = 1;
            monitor_dives = monitor_dives + ceil(META.COTS_cull_reeflist.nb_sites(I)*1); % additional no. of dives * no. of cull sites
            % assume that we need additional monitor effort equivalent to 1
            % extra dive per cull site to reach 100% confidence that all CoTS
            % are detected
        end

        % Case 2: look up relocation travel cost if this reef was selected as a relocation candidate in f_makeReefList_NEW.
        % Computed here (after monitor_dives) so both upfront costs can be checked and debited together.
        reloc_idx = find(relocation_dist_km.reef_ID == this_reef_ID, 1);
        reloc_travel_dv = 0; % default: no relocation travel cost
        if ~isempty(reloc_idx)
            reloc_travel_dv = km_to_dives(relocation_dist_km.dist_km(reloc_idx));
        end

        % Suki Apr 2026: Check global and regional budget for all upfront costs together:
        % relocation travel (case 2) + additional monitoring (bleaching cat 2).
        % Checking together prevents debiting travel then skipping the reef for insufficient monitoring budget.
        upfront_dives = reloc_travel_dv + monitor_dives;
        if upfront_dives > 0
            % Check global budget first
            if remaining_dives < upfront_dives
                n = n + 1;
                continue % skip this reef - not enough global dives for upfront costs
            end

            % Check regional budget (own region only - no reallocation)
            if remaining_region_dives(region_idx) < upfront_dives
                n = n + 1;
                continue % skip this reef - not enough regional dives for upfront costs
            end

            % Budget confirmed: subtract all upfront costs from global and regional totals
            remaining_dives = remaining_dives - upfront_dives;
            remaining_region_dives(region_idx) = remaining_region_dives(region_idx) - upfront_dives;

            % Record relocation travel (from = original bleaching cat>3 reef; to = this relocation candidate)
            if reloc_travel_dv > 0
                control_records.reef_travel_dist_km = [control_records.reef_travel_dist_km; relocation_dist_km.dist_km(reloc_idx)];
                control_records.reef_travel_dives   = [control_records.reef_travel_dives;   reloc_travel_dv];
                control_records.reef_travel_from    = [control_records.reef_travel_from;    relocation_dist_km.original_reef_ID(reloc_idx)];
                control_records.reef_travel_to      = [control_records.reef_travel_to;      this_reef_ID];
            end
        end
    end

    %% MANTA TOW COTS DENSITY AND SITE-LEVEL TRIGGER

    if META.doing_COTS_heat_stress_detectibility == 1 && additional_monitoring == 0 % Cook et al. in review
        % Apply scalar heat stress factor to the aggregate MTC (correct: scalar multiplier on already-aggregated value)
        this_reef_COTS_per_tow_per_site_reduced = this_reef_COTS_per_tow_per_site .* exp(-0.25 * all_DHWs(I));
        %Scalar unit. Proportion of CoTS per increase in DHW (Ranges from 1-0.02 as a decay function)

        % Find sites with manta tows (reduced detectability due to DHW) above ecological threshold (ET)
        sites_over_ET = find(this_reef_COTS_per_tow_per_site_reduced > current_reef_ET(I));

    else
        % Find sites with manta tows above ecological threshold (ET)
        % assume that with additional monitoring we can detect CoTS as if
        % conditions were normal

        sites_over_ET = find(this_reef_COTS_per_tow_per_site > current_reef_ET(I));

    end

    % Suki Apr 2026: Sort sites_over_ET by descending CoTS density for optimal effort allocation
    if META.use_site_skip_find == 1
        if ~isempty(sites_over_ET)
            % this_reef_COTS_per_tow_per_site is [n_sites x 1] after the aggregate-then-convert fix — no extra sum needed
            [~, sort_order] = sort(this_reef_COTS_per_tow_per_site(sites_over_ET), 'descend');
            sites_over_ET = sites_over_ET(sort_order);
        end
    end

    %% PERFORM COTS CULLING ON SITES ABOVE ET WHEN REMAINING DIVES AVAILABLE

    % We go only to those sites
    if ~isempty(sites_over_ET) && remaining_dives > 0 % if there are sites to treat, and dives remaining, do control

        last_reef_COTScontrolled = this_reef_ID; % set this reef as the last controlled one for the next step; set here because we don't know when we run out of dives
        record_site_post_densities = this_reef_COTS_densities_per_site; % init

        site = 1 ; % iterator of culled sites
        record_nb_dives = 0 ; % track the number of dives to control all sites of a reef

        sites_culled = 0 ; % Suki Apr 2026: track number of actually culled sites
        global_budget_hit = false ; % Suki Apr 2026: flag indicating if global budget was hit
        region_budget_hit = false(size(remaining_region_dives)); % Suki Apr 2026: flag indicating if regional budget was hit for this reef
        first_skipped_site = 0; 

        while site <= length(sites_over_ET) && remaining_dives > 0 % go through all sites above ET

            ctrl_site = sites_over_ET(site) ;
            COTS_per_tow_this_site = this_reef_COTS_per_tow_per_site(ctrl_site); % scalar MTC for this site (aggregate-converted)

            % Pre-compute age-class proportions and per-class heat stress culling effects for this site.
            % These are used in both the budget (additional_culling) and culling factor (culling_factor_heat) sections.
            % age_props: proportion of each adult age class in the detectable density at this site.
            % culling_dhw_effects: age-specific culling kill-rate multiplier under current DHW.
            detected_this_site        = this_reef_COTS_densities_per_site_detected(ctrl_site, :); % [1 x n_adult_ages]
            total_detected_this_site  = this_reef_detected_total_per_site(ctrl_site);              % scalar
            if total_detected_this_site > 0
                age_props = detected_this_site / total_detected_this_site; % [1 x n_adult_ages], sums to 1
            else
                age_props = ones(1, numel(detected_this_site)) / numel(detected_this_site);
            end
            
            if META.doing_COTS_heat_stress_detectibility == 1
                culling_dhw_effects = exp(META.log_culling_dhw_effects(META.COTS_adult_min_age:end)' * all_DHWs(I)); % [n_adult_ages x 1]
            end

            additional_culling = 0;
            if is_adaptive && bleaching_category(I) == add_cull_cat
                % calculate adjusted CoTS per tow for this site under heat stress.
                % The weighted correction term represents the extra MTC needed to achieve the same
                % kill as without heat stress, averaged across age classes by their proportion.
                % caveat: assuming CoTS can ultimately be culled to ET when in reality that might not be the case
                additional_culling = 1;
                COTS_per_tow_this_site_sizeclass = COTS_per_tow_this_site * age_props;
                adjusted_COTS_per_tow_this_site = COTS_per_tow_this_site_sizeclass + current_reef_ET(I).*(1./culling_dhw_effects - 1);
                control_dives = ceil(4.18 * (sum(adjusted_COTS_per_tow_this_site)/0.015)^0.667);

            else                
                    % First need to check if we've got enough remaining dives for this site
                    % Number of control dives required for culling to this threshold
                    % Karlo: here 8 is the standard number of divers, I use team dives per site, using different boats would require a rewrite
                    % Ceil to round, as the whole team will finish the dive on the same site
                    control_dives = (ceil((4.18)*(COTS_per_tow_this_site/0.015)^0.667)); %indiv dives
                    % YM: this_site_tow/0.015 converts number of CoTS per tow into a density of CoTS per hectare
        
                    % control_dives=ceil(166.7*(COTS_per_tow_this_site/0.015)^0.665); % Fletcher et al. 2021, divided by 40 assuming 40 min bottom time?
                    % which is equivalent to: control_dives=ceil(exp(5.116+0.665*log(COTS_per_tow_this_site/0.015))); % Fletcher et al. 2021, divided by 40 assuming 40 min bottom time?
                    % Consider expressing all CoTS densities per hectare and convert for coral onsumption?
            end
  
            % Check global and regional budgets before committing to this site.
            % Use global_insufficient and region_insufficient for the skip decision
            % Global_budget_hit and region_budget_hit are only used for outer-loop break. 
            % Once budget is hit at any site, we stop processing the next reef. this will be the last reef we work on (globally or for this region)
            global_insufficient = remaining_dives < control_dives;
            region_insufficient = remaining_region_dives(region_idx) < control_dives;

            if global_insufficient
                global_budget_hit = true;
            end
            if region_insufficient
                region_budget_hit(region_idx) = true;
            end

            if global_insufficient || region_insufficient
                if META.use_site_skip_find
                    if first_skipped_site == 0
                        first_skipped_site = ctrl_site;
                    end
                    site = site + 1;
                    continue % skip this site, try a cheaper one
                else
                    break % original: stop processing sites on this reef
                end
            end

            % subtract dives from global and regional totals when there is still enough budget
            remaining_region_dives(region_idx) = remaining_region_dives(region_idx) - control_dives;
            remaining_dives = remaining_dives - control_dives;

            % If allowed to cull (May 2026): account for non-linearity in conversion CoTS density <-> CoTS per tow
            COTS_per_tow_culled_this_site = COTS_per_tow_this_site - current_reef_ET(I);

            if META.doing_COTS_heat_stress_detectibility == 1 && additional_culling == 0
                COTS_per_tow_culled_this_site = sum(COTS_per_tow_this_site * age_props .* culling_dhw_effects) - current_reef_ET(I);
            end

            % Clamp to zero: heat-stress-reduced detectibility can push effective CoTS per tow below ET,
            % yielding a negative value that causes complex arithmetic in f_convert_CoTS_tow_2_density
            COTS_per_tow_culled_this_site = max(0, COTS_per_tow_culled_this_site);

            % Convert the loss of CoTS per tow into loss of CoTS density - this time use a deterministic prediction (ie, NB_tows = 0)
            % to avoid negatives down the line
            COTS_density_culled_this_site = f_convert_CoTS_tow_2_density(COTS_per_tow_culled_this_site, 0, META.total_area_cm2);

            % always cap to maximum 1 just in case culled > total detected
            culling_factor = min(1, COTS_density_culled_this_site/this_reef_detected_total_per_site(ctrl_site));

            record_site_post_densities(ctrl_site,META.COTS_adult_min_age:end) = (this_reef_COTS_densities_per_site(ctrl_site,META.COTS_adult_min_age:end))*(1 - culling_factor);
            
            sites_culled = sites_culled + 1 ; % Suki Apr 2026: increment actual cull count

            % record skip-then-find event (a site was too expensive, but a cheaper site on the same reef was found)
            if first_skipped_site > 0
                control_records.site_skip_count = control_records.site_skip_count + 1;
                control_records.site_skip_reef_ID = [control_records.site_skip_reef_ID; this_reef_ID];
                control_records.site_skip_first_site = [control_records.site_skip_first_site; first_skipped_site];
                control_records.site_skip_success_site = [control_records.site_skip_success_site; ctrl_site];
                first_skipped_site = 0; % reset for potential next skip-then-find on same reef
            end

            record_nb_dives = record_nb_dives + control_dives;
            site = site + 1 ;
        end

        % Suki Apr 2026: Zero out any regions that hit their budget during site processing
        % so the reef-level check skips future reefs in those regions
        hit_regions = find(region_budget_hit);
        for hr = 1:length(hit_regions)
            control_records.region_break_reef_ID = [control_records.region_break_reef_ID; this_reef_ID];        % reef where the regional budget ran out mid-culling
            control_records.region_break_region = [control_records.region_break_region; hit_regions(hr)];        % which region (1=FN, 2=N, 3=C, 4=S) ran out of dives
            control_records.region_break_remaining = [control_records.region_break_remaining; remaining_region_dives(hit_regions(hr))]; % leftover regional dives at the point of running out (diagnostic)
        end
        remaining_dives = remaining_dives - sum(remaining_region_dives(region_budget_hit)); % deduct stranded regional dives from global pool to keep pools consistent
        remaining_region_dives(region_budget_hit) = 0; % zero out exhausted regions so future reefs in those regions are skipped

        % Suki Apr 2026: Only record reef as culled if at least one site was actually controlled
        if sites_culled > 0

            %% After culling all sites above ET, we record the results and update the COTS density of the reef

            % track records for the reef
            control_records.culled_reef_ID = [control_records.culled_reef_ID ; this_reef_ID];                    % ID of every reef where at least one site was culled this timestep
            control_records.nb_culled_sites = [control_records.nb_culled_sites ; sites_culled];                  % number of sites actually culled on this reef (may be < total sites above ET if budget ran out)
            record_reef_post_densities = sum(record_site_post_densities,1)/size(record_site_post_densities,1);   % mean CoTS density per site after culling (averaged across all sites of the reef)
            record_reef_density_culled = sum(squeeze(COTS_all_densities(I, 1, META.COTS_adult_min_age:end)) - record_reef_post_densities(META.COTS_adult_min_age:end)'); % total adult CoTS removed: pre-cull minus post-cull density, summed across age classes
            control_records.culled_density_reef = [control_records.culled_density_reef ; record_reef_density_culled]; % CoTS removed on this reef (adults per 400 m2, one value per culled reef)
            control_records.nb_dives = [control_records.nb_dives ; record_nb_dives];                            % total dives spent culling this reef (summed across all culled sites)

            % Record new density of adults for that reef
            COTS_all_densities(I, META.COTS_adult_min_age:end) = record_reef_post_densities(META.COTS_adult_min_age:end); % update the model state: CoTS density on this reef is now the post-cull value

            % track additional culling effort (strategy 19, bleaching category 3: extra dives needed due to heat-stress-reduced kill rate)
            if additional_culling == 1
                control_records.additional_cull_ID = [control_records.additional_cull_ID; this_reef_ID];        % ID of reefs where extra dives were used to compensate for heat-stress reduced culling efficiency
                control_records.nb_additional_cull = control_records.nb_additional_cull + 1;                    % running count of such reefs this timestep
            end
        end
    end

    %% Suki March 2026: track effort and monitoring at the reef level (recorded for every visited reef, culled or not)
    control_records.nb_control_sites = [control_records.nb_control_sites; META.COTS_cull_reeflist.nb_sites(I)]; % total number of cull sites on this reef (from the reef list, independent of how many were actually culled)
    if is_adaptive && additional_monitoring == 1
        % additional_monitoring == 1 only for bleaching category add_monitor_cat reefs: culling proceeds but extra survey dives are added
        control_records.additional_monitor_ID    = [control_records.additional_monitor_ID; this_reef_ID];     % ID of category-2 reefs where extra monitoring dives were added on top of culling
        control_records.nb_monitored_reefs = control_records.nb_monitored_reefs + 1;                            % running count of reefs with extra monitoring (category 2 only here; monitor-only reefs are counted in the monitor-only block above)
        control_records.additional_monitor_dives  = [control_records.additional_monitor_dives; monitor_dives];  % extra monitoring dives spent on this reef (1 dive per cull site)
        control_records.total_monitor_dives = sum(control_records.skip_region_monitor_dives) + sum(control_records.skip_local_monitor_dives) + sum(control_records.additional_monitor_dives); % total monitoring dives spent across all monitored reefs this timestep
    end

    control_records.visited_reef_ID = [control_records.visited_reef_ID; this_reef_ID];                          % ID of every reef where a manta tow was performed (includes culled, monitor-only, and category-2 reefs)
    if META.doing_COTS_heat_stress_detectibility == 1 && additional_monitoring == 0 
        control_records.visited_manta_tow = [control_records.visited_manta_tow; sum(this_reef_COTS_per_tow_per_site_reduced, 'all')]; % observed CoTS per tow summed across all sites, with detectability reduced by heat stress (Cook et al.) - what the survey team actually sees
    else
        control_records.visited_manta_tow = [control_records.visited_manta_tow; sum(this_reef_COTS_per_tow_per_site, 'all')];         % observed CoTS per tow summed across all sites, full detectability (no heat stress, or additional monitoring overrides reduction)
    end
    %% Move on to the next reef or exit if out of dives or list exhausted
    visited_reefs = visited_reefs + 1; % up the count of controlled reefs
    n = n + 1; % go to the next reef

    % Suki Apr 2026: If global budget was hit during site processing, this is the last reef.
    % Don't travel to another reef - break the global loop.
    if exist('global_budget_hit', 'var') && global_budget_hit
        control_records.global_break_reef_ID = [control_records.global_break_reef_ID; this_reef_ID];            % reef where the global dive budget was exhausted (last reef visited)
        control_records.global_break_remaining = [control_records.global_break_remaining; remaining_dives];     % leftover global dives at that point (diagnostic; should be ~0)
        remaining_dives = 0; % set remaining dives to 0 to exit the main loop
        break
    end
end

%% End-of-timestep summaries
control_records.remaining_dives = remaining_dives;                                                              % global dives left unspent at end of timestep (normally 0 - budget fully used)
if META.use_regional_dives
    control_records.remaining_region_dives = remaining_region_dives;                                            % per-region dives left unspent at end of timestep (normally 0)
else
    control_records.remaining_region_dives = [];                                                                % not tracked in global-only mode
end
control_records.nb_culled_reefs = length(control_records.culled_reef_ID);                                       % total number of reefs where culling occurred (at least one site culled)
control_records.nb_visited_reefs = visited_reefs;                                                               % total reefs visited (manta tow performed), including culled, monitor-only, and category-2 reefs
control_records.nb_culled_COTS = uint64(sum(control_records.culled_density_reef.*control_records.nb_culled_sites)*1e5/400); % total CoTS removed across the GBR this timestep (1 site = 1e5 m2 = 10 ha, Skinner et al. 2024)
control_records.total_travel_dives   = sum(control_records.reef_travel_dives) + sum(control_records.port_travel_dives); % total dives lost to travel: relocation (case 2) + inter-port redistribution (case 1)
control_records.total_expended_dives = sum(control_records.nb_dives) + control_records.total_monitor_dives;    % all dives spent this timestep: culling dives + monitoring dives combined


% %% TEMP FOR DEBUGGING (see also last line for recording last_reef_COTScontrolled in RESULT)
% % REEF.last_reef_COTScontrolled = last_reef_COTScontrolled;
% % Remember t was incremented above by 1 so need to overwrite t, not t+1
% RESULT.COTS_all_densities(control_records.culled_reef_ID,t,META.COTS_adult_min_age:end) = COTS_all_densities(control_records.culled_reef_ID, META.COTS_adult_min_age:end);
% RESULT.COTS_culled_density(control_records.culled_reef_ID,t-1) = control_records.culled_density;
% RESULT.COTS_culled_reefs(1:META.nb_reefs,t-1) = ismember(META.reef_ID,control_records.culled_reef_ID);
% RESULT.COTS_control_remaining_dives(1,t-1)= remaining_dives;
