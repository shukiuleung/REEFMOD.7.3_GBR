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
    f_COTS_control_NEW(META, COTS_all_densities, last_reef_COTScontrolled, total_coral_pct2D, COTS_larval_output, all_DHWs)

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
current_COTS_per_tow = (0.22/0.6)*sum(current_COTS_densities(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);

% thisboatorder = META.boatProperties.fixedOrder(1); % each boat to visit reefs in fixed (1) or random (0) order but only considering first boat for now
timestep_dives = sum(META.boatProperties.totalInidvDives); % full dive quota for timestep

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
ReefList = f_makeReefList_NEW(META, current_COTS_densities, current_reef_ET, COTS_densities_per_site, total_coral_pct2D, COTS_larval_output, last_reef_COTScontrolled, all_DHWs);

% ReefList must be a selection of META.reef_ID sorted in decreasing order of priority for control
% Keep track record of key control variables
control_records=struct('culled_reef_ID',[], 'culled_density_reef',[], 'culled_density_total',[],...
    'nb_dives',[],'nb_culled_sites',[], 'nb_culled_reefs',[], 'nb_visited_reefs',[], ...
    'nb_control_sites' , [], 'nb_unvisited_reefs', []);

visited_reefs = 0; % counter to keep track of how many reefs were culled
n = 1; % Start with the first reef on the list
remaining_dives = timestep_dives; %full dive quota at the beginning of a timestep

while remaining_dives > 0 && n <= length(ReefList) %while there are dives remaining and the list not exhausted
    % General description from KH: get first reef ID from the list, check if reef-level trigger is used, then use it;
    % check if site-level trigger is valid at any sites. Go only to those sites, recheck after every run to see if trigger
    % still valid, if not, increase the current_reef by one until out of quota

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

    % Extract COTS density for all sites
    this_reef_COTS_densities_per_site = COTS_densities_per_site{I,1};

    this_reef_COTS_per_tow_per_site = (0.22/0.6).*(this_reef_COTS_densities_per_site(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end));


    if META.doing_COTS_heat_stress_detectibility == 1 % Cook et al. in review

        this_reef_COTS_per_tow_per_site_reduced = (0.22/0.6).*(this_reef_COTS_densities_per_site(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end)).*exp(-0.25 * all_DHWs(I));
        %Scalar unit. Proportion of CoTS per increase in DHW (Ranges from 1-0.02 as a decay function)

        % Find sites with manta tows (reduced detectability due to DHW) above ecological threshold (ET)
        sites_over_ET = find(sum(this_reef_COTS_per_tow_per_site_reduced,2) > current_reef_ET(I));

    else
        % Find sites with manta tows above ecological threshold (ET)

        sites_over_ET = find(sum(this_reef_COTS_per_tow_per_site,2) > current_reef_ET(I));

    end


    % We go only to those sites
    if ~isempty(sites_over_ET) && remaining_dives>0 % if there are sites to treat, and dives remaining, do control

        last_reef_COTScontrolled = this_reef_ID; % set this reef as the last controlled one for the next step; set here because we don't know when we run out of dives
        record_site_post_densities = this_reef_COTS_densities_per_site;

        site = 1 ; % iterator of culled sites
        record_nb_dives = 0 ; % track the number of dives to control all sites of a reef

        while site <= length(sites_over_ET) && remaining_dives>0 % go through all sites above ET

            ctrl_site = sites_over_ET(site) ;
            COTS_per_tow_this_site = this_reef_COTS_per_tow_per_site(ctrl_site, :); % original detectability

            % First need to check if we've got enough remaining dives for this site
            % Number of control dives required for culling to this threshold
            % Karlo: here 8 is the standard number of divers, I use team dives per site, using different boats would require a rewrite
            % Ceil to round, as the whole team will finish the dive on the same site
            control_dives=(ceil((4.18)*(sum(COTS_per_tow_this_site)/0.015)^0.667)); %indiv dives
            % YM: this_site_tow/0.015 converts number of CoTS per tow into a density of CoTS per hectare

            % control_dives=ceil(166.7*(COTS_per_tow_this_site/0.015)^0.665); % Fletcher et al. 2021, divided by 40 assuming 40 min bottom time?
            % which is equivalent to: control_dives=ceil(exp(5.116+0.665*log(COTS_per_tow_this_site/0.015))); % Fletcher et al. 2021, divided by 40 assuming 40 min bottom time?
            % Consider expressing all CoTS densities per hectare and convert for coral onsumption?

            remaining_dives = remaining_dives-control_dives; % subtract perfomed dives from total (can be negative)

            if remaining_dives < 0 % break if not enough dives available for this site

                remaining_dives = 0 ; % set to 0 to break the global while loop after this while loop
                break % stop control on this reef (breaks this while loop)
            end

            if META.doing_COTS_heat_stress_detectibility == 1

                % If allowed to cull, apply heat stress detectability reduction to calculate culling factor = fewer CoTS are removed than required
                culling_factor = sum(COTS_per_tow_this_site .* (exp(META.log_culling_dhw_effects(META.COTS_adult_min_age:end)' * all_DHWs(I))))/current_reef_ET(I); % this is the reduction needed to take CoTS per tow to the ET
            
            else

                % If not applying heat stress detectability reduction, use standard culling factor
                culling_factor = sum(COTS_per_tow_this_site)/current_reef_ET(I); % this is the reduction needed to take CoTS per tow to the ET

            end            
            
            record_site_post_densities(ctrl_site,META.COTS_adult_min_age:end) = this_reef_COTS_densities_per_site(ctrl_site,META.COTS_adult_min_age:end)./culling_factor;
            
            record_nb_dives = record_nb_dives + control_dives;
            site = site + 1 ;
        end

        % track records for the reef
        control_records.culled_reef_ID = [control_records.culled_reef_ID ; this_reef_ID];
        control_records.nb_culled_sites = [control_records.nb_culled_sites ; site-1]; % number of sites is site-1, whether we break before or cull all sites
        record_reef_post_densities = sum(record_site_post_densities,1)/size(record_site_post_densities,1); % density after culling: sum across sites and divide by number of sites
        record_reef_density_culled = sum(COTS_all_densities(I, META.COTS_adult_min_age:end) - record_reef_post_densities(META.COTS_adult_min_age:end));
        control_records.culled_density_reef = [control_records.culled_density_reef ; record_reef_density_culled]; % as total number of adults per 400 m2
        control_records.nb_dives = [control_records.nb_dives ; record_nb_dives];

        % Suki: track records related to survey effort
        control_records.nb_control_sites = [control_records.nb_control_sites; META.COTS_cull_reeflist.nb_sites(this_reef_ID)];

        % Record new density of adults for that reef
        COTS_all_densities(I, META.COTS_adult_min_age:end) = record_reef_post_densities(META.COTS_adult_min_age:end);
    end

    visited_reefs = visited_reefs + 1;%up the count of controlled reefs
    n = n + 1; % go to the next reef
end

% % control_records.remaining_dives = remaining_dives; % record the remaining control effort when simulating less then 3,806 reefs - useless because always 0!!
control_records.nb_culled_reefs = length(control_records.culled_reef_ID);
control_records.nb_visited_reefs = visited_reefs;
control_records.culled_density_total = sum(control_records.culled_density_reef.*control_records.nb_culled_sites)*1e5/400;
% This is the total number of CoTS killed, ie the summed culled density across the culled sites of each reef, across all reefs
% This assumes 1 site ~ 1e5 m2 = 0.1 km2 (10 hectares, Skinner et al. 2024)

% %% TEMP FOR DEBUGGING (see also last line for recording last_reef_COTScontrolled in RESULT)
% % REEF.last_reef_COTScontrolled = last_reef_COTScontrolled;
% % Remember t was incremented above by 1 so need to overwrite t, not t+1
% RESULT.COTS_all_densities(control_records.culled_reef_ID,t,META.COTS_adult_min_age:end) = COTS_all_densities(control_records.culled_reef_ID, META.COTS_adult_min_age:end);
% RESULT.COTS_culled_density(control_records.culled_reef_ID,t-1) = control_records.culled_density;
% RESULT.COTS_culled_reefs(1:META.nb_reefs,t-1) = ismember(META.reef_ID,control_records.culled_reef_ID);
% RESULT.COTS_control_remaining_dives(1,t-1)= remaining_dives;
