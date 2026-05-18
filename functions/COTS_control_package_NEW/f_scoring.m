% ----------------------------------------------------------------------------------------------------------------------
% Suki Leung, MSEL, March 2026 
% Based on Tina's f_makeReefList_NEW cases 25-30 that omputes a combined CoTS-risk / key-source-reef (KSR) score
% CCIP scenario rqeuire ranking by different criteria:
% Default (1):
% - previous cull record
% - previous manta tow record
% - Reef priority (T > P > N)
% - bleaching category
% - distance from relocation reefs
%
% Perfect knowledge A (2):
% - replace previous cull record with CoTS risk rank
% - replace manta tow record with CoTS_manta_tow
%
% Perfect knowledge B (3):
% - replace previous cull record with CoTS benefits
% - replace manta tow record with CoTS_manta_tow
%
% New reef priority (4,5):
% - replace previous cull record with CoTS risk rank
% - replace reef priority with KSR score
% - include bonuses for refugia and hotspots
% - scenario 4 = refugia + BZ, scenario 5 = hotspot + BZ
%
% Non-priority reef scoring (6):
% - score non-priority reefs by CoTS risk and bleaching category (to help identify which non-priority reefs might be worth monitoring or culling if resources allow)
% 
% INPUTS
%   META - model metadata struct (uses: connectivity, KSR_reefsizes, COTS_cull_reeflist.AIMS_sector, cots_risk_rank, refugia, hotspot, COTS_cull_reeflist.GreenZone, COTS_reefs2cull_strat)
%   CONNECT_CORAL - struct array of connectivity matrices for KSR
%   t - current timestep index
%   total_coral_pct2D - total coral cover per reef (row vector, length = META.nb_reefs)
%   priority_list - ordered reef IDs to score (combined target + priority + nonpriority)
%   scenario - integer (1-4) indicating which scenario to run (see above)
%
% OUTPUT
%   sorted_ID - reef IDs sorted from highest to lowest combined score
% ----------------------------------------------------------------------------------------------------------------------
function [sorted_ID] = f_scoring(META, CONNECT_CORAL, t, total_coral_pct2D, priority_list, control_records_prev, current_COTS_per_tow, scenario, relocate_reef)

[~, pl_rows] = ismember(priority_list, META.COTS_cull_reeflist.Reef_ID); % row indices into META.COTS_cull_reeflist (should be equivalent to full reef list in simulation)

%% Initialise empty control_records_prev if this is the first timestep (no prior CoTS control)
if isempty(control_records_prev)
    control_records_prev = struct( ...
        'culled_reef_ID',      zeros(META.nb_reefs, 1), ...
        'culled_density_reef', zeros(META.nb_reefs, 1), ...
        'visited_reef_ID',     zeros(META.nb_reefs, 1), ...
        'visited_manta_tow',   zeros(META.nb_reefs, 1));
end

%% ============================================================
%  SCORE ALL REEFS (full domain), THEN FILTER TO PRIORITY LIST
%  ============================================================

%% Previous cull record
culled_density_prev = zeros(META.nb_reefs, 1);
for i = 1:META.nb_reefs
    reef_ID = META.COTS_cull_reeflist.Reef_ID(i);
    idx = find(control_records_prev.culled_reef_ID == reef_ID, 1);
    if ~isempty(idx)
        culled_density_prev(i) = control_records_prev.culled_density_reef(idx);
    end
end

culled_rank_n = normalize(culled_density_prev, 'range'); % normalise to 0-1

%% Previous manta tow record
manta_tow_prev = zeros(META.nb_reefs, 1);
for i = 1:META.nb_reefs
    reef_ID = META.COTS_cull_reeflist.Reef_ID(i);
    idx = find(control_records_prev.visited_reef_ID == reef_ID, 1);
    if ~isempty(idx)
        manta_tow_prev(i) = control_records_prev.visited_manta_tow(idx);
    end
end

manta_tow_rank_n = normalize(manta_tow_prev, 'range'); % normalise to 0-1

%% Reef priority (static)
reef_priority = string(META.COTS_cull_reeflist.reef_type); % T > P > N
reef_priority_rank = zeros(size(reef_priority));
reef_priority_rank(reef_priority == 'T') = 3;
reef_priority_rank(reef_priority == 'P') = 2;
reef_priority_rank(reef_priority == 'N') = 1;
reef_priority_rank = normalize(reef_priority_rank, 'range'); % normalise to 0-1

%% Bleaching category
bleaching_category = META.COTS_cull_reeflist.bleaching_category; % 1-5
bleaching_category_rank = 1 - normalize(bleaching_category, 'range'); % normalise to 0-1 - inverted because lower category is better

%% Distance from target reef for relocation
% only applies if relocate_reef is specified (non-zero), otherwise all distances are set to 0 and won't contribute to the score
% relocate_reef is a single reef ID - original target reef that is cat 4 or above and needs to be relocated
if relocate_reef ~= 0
    distances = META.distance_matrix(relocate_reef, :); % extract distances to all reefs
    distance_rank_n = (1 - normalize(distances(:), 'range')); % normalise to 0-1 - inverted because closer is better
else
    distance_rank_n = zeros(META.nb_reefs, 1); % if no relocation reef, set all distances to 0
end

%% Predicted CoTS manta tow 
cots_per_tow_rank_n = normalize(current_COTS_per_tow(:), 'range'); % normalise to 0-1, higher CoTS = higher priority

%% CoTS benefits
% Each year maps to two timesteps (e.g. 2019 -> t=23 and t=24), but the
% benefits table only stores odd t values (23, 25, 27...). Round t down to
% the nearest stored timestep before lookup.
t_benefits = t - mod(t - 23, 2); % t=23->23, t=24->23, t=25->25, t=26->25, etc.
cots_benefits = META.cots_benefits(META.cots_benefits.t == t_benefits, :);
cots_benefits_rank_n = normalize(cots_benefits.meanDelta, 'range');

%% CoTS risk rank (static, from counterfactual)
cots_risk_rank = META.cots_risk_rank;

%% Key Source Reef (KSR) score via connectivity-weighted larval enrichment
con_yr_id = META.connectivity.CORAL_sequence{2,t};
con = full(CONNECT_CORAL(con_yr_id).LINK_STRENGTH);

cover = total_coral_pct2D';  % need coral cover
critic = 20;                 % reefs below 20% coral cover need larval input
alpha  = 1;

[row, col] = size(con);
out = zeros(row, 13);
region = META.COTS_cull_reeflist.AIMS_sector;
out(:,12) = region(:,1);

% Weights based on current cover (need = how much a reef requires input)
need = zeros(1, row);
for i = 1:row
    if cover(1,i) < critic
        need(1,i) = (critic - cover(1,i)) ./ critic;
    end
end

simul_KSR_reefsizes = META.KSR_reefsizes(META.reef_ID);

% create weighted contributions to connectivity
larvalinputs     = (con' .* cover .* simul_KSR_reefsizes)'; % multiples each connectivity source coral cover and reef size
%larvalinput_reef=sum(larvalinputs,2); % vector of total inputs to each reef
larvalinput_reef = sum(larvalinputs, 1);
larvalflux_reef  = larvalinput_reef ./ simul_KSR_reefsizes; % convert to a flux based on size

for i = 1:row % calc enrichment done by each source
    enrichunweight = 0;
    enrichweight   = 0;
    for j = 1:col
        if con(i,j) > 0
            inputflux     = larvalinputs(i,j) ./ simul_KSR_reefsizes(1,j);  % convert to a flux
            propinputflux = inputflux ./ larvalflux_reef(1,j); % express input as proprotion of total flux to that reef
            enrichunweight = enrichunweight + propinputflux; %create total
            enrichweight   = enrichweight   + (propinputflux .* need(1,j)); % weight by need at the sink
        end
    end
    out(i,1) = i;
    out(i,2) = enrichunweight;
    out(i,3) = enrichweight;
    out(i,4) = need(1,i);
end

 % identify number of reefs each source contacts and then the weighted
 % number which ignores those that don't need it
for i = 1:row
    numreefs       = 0;
    numreefsneeded = 0;
    for j = 1:col
        if larvalinputs(i,j) > 0 && j ~= i
            numreefs = numreefs + 1;
            if cover(1,j) <= critic
                numreefsneeded = numreefsneeded + 1;
            end
        end
    end
    out(i,5) = numreefs;
    out(i,6) = numreefsneeded;
end

% convert the larval enrichment and numreefsneeded into z scores
out(:,7) = normalize(out(:,3), 'range');
out(:,8) = normalize(out(:,6), 'range');
out(:,9) = out(:,7) + (alpha .* out(:,8));

% Re-scale by AIMS region
regions = unique(region);
for i = 1:length(regions)
    dat = out(out(:,12) == regions(i,1), :); %extract those cells
    dat(:,10) = normalize(dat(:,3), 'range');
    dat(:,11) = normalize(dat(:,6), 'range');
    dat(:,13) = dat(:,10) + (alpha .* dat(:,11));
    for j = 1:size(dat, 1)
        out(dat(j,1), :) = dat(j,:); % puts them back in correct place
    end
end

% KSR rank (rank 1 = most important source)
KSR_score = out(:,13);  %KSR final regionalised score
[~, order] = sort(KSR_score, 'descend');
KSR_rank = zeros(size(KSR_score));
KSR_rank(order) = 1:numel(KSR_score);

%% Combine CoTS risk and KSR into a single score (both normalised, inverted so 1 = most important)
cots_risk_rank_n = 1 - normalize(cots_risk_rank, 'range'); 
KSR_rank_n       = 1 - normalize(KSR_rank,       'range');  

%% Final ranks for scoring
culled_rank_n         = culled_rank_n(pl_rows);
manta_tow_rank_n      = manta_tow_rank_n(pl_rows);
reef_priority_rank    = reef_priority_rank(pl_rows);
bleaching_category_rank = bleaching_category_rank(pl_rows);
distance_rank_n       = distance_rank_n(pl_rows);
cots_per_tow_rank_n   = cots_per_tow_rank_n(pl_rows);
cots_benefits_rank_n  = cots_benefits_rank_n(pl_rows);
cots_risk_rank_n      = cots_risk_rank_n(pl_rows);
KSR_rank_n            = KSR_rank_n(pl_rows);

clearvars -except META priority_list pl_rows scenario culled_rank_n manta_tow_rank_n reef_priority_rank bleaching_category_rank distance_rank_n cots_per_tow_rank_n cots_benefits_rank_n cots_risk_rank_n KSR_rank_n

%% weights of criteria - can be changed easily
w_cull_record = 0.4;
w_manta_record = 0.2;
w_priority = 0.15;
w_bleaching_category = 0.15;
w_distance = 0.15;
w_predicted_manta = 0.2;
w_cots_benefits = 0.4;
w_cots_risk = 0.4;
w_KSR = 0.5;

%% Prepare masks for refugia and hotspots (built over full domain, then extracted to priority_list)
refugia_mask = false(3806, 1);
hotspot_mask = false(3806, 1);
refugia_mask(META.refugia) = true;
hotspot_mask(META.hotspot) = true;

refugia_mask = refugia_mask(priority_list);  % extract to priority_list
hotspot_mask = hotspot_mask(priority_list);  % extract to priority_list

% green/blue already sized to priority_list (from META.COTS_cull_reeflist)
green_mask = META.COTS_cull_reeflist.GreenZone(:) == 1;
blue_mask  = ~green_mask;

green_mask = green_mask(pl_rows);
blue_mask = blue_mask(pl_rows);

%% Add bonuses according to strategy
bonus_thermal = 0.1;  % easy to tune
bonus_zoning  = 0.1;  % easy to tune

switch scenario
    case 1   % default AC
        score = culled_rank_n * w_cull_record + manta_tow_rank_n * w_manta_record + reef_priority_rank * w_priority + bleaching_category_rank * w_bleaching_category + distance_rank_n * w_distance;
    case 2   % perfect knowledge a - case 27 in f_makeReefList_NEW
        score = cots_risk_rank_n * w_cots_risk + cots_per_tow_rank_n * w_predicted_manta + reef_priority_rank * w_priority + bleaching_category_rank * w_bleaching_category + distance_rank_n * w_distance;
    case 3   % perfect knowledge b - case 28 in f_makeReefList_NEW
        score = cots_benefits_rank_n * w_cots_benefits + cots_per_tow_rank_n * w_predicted_manta + reef_priority_rank * w_priority + bleaching_category_rank * w_bleaching_category + distance_rank_n * w_distance;
    case 4   % refugia + BZ
        w_cots_risk = 0.5;
        w_KSR  = 0.5;
        bonus_thermal = 0.1;
        bonus_zoning  = 0.1;
        score = w_cots_risk*cots_risk_rank_n + w_KSR*KSR_rank_n + bonus_thermal * refugia_mask + bonus_zoning * blue_mask;

    case 5  % hotspot + BZ
        w_cots_risk = 0.5;
        w_KSR  = 0.5;
        bonus_thermal = 0.1;
        bonus_zoning  = 0.1;
        score = w_cots_risk*cots_risk_rank_n + w_KSR*KSR_rank_n + bonus_thermal * hotspot_mask + bonus_zoning * blue_mask;
    case 6 % scoring of non-priority reefs
        score = cots_risk_rank_n * 0.8 + bleaching_category_rank * 0.2;
end

%% Sort priority list by score and return
[~, idx] = sort(score, 'descend'); % score is already indexed to priority_list, sort in descending order
sorted_ID = priority_list(idx);
