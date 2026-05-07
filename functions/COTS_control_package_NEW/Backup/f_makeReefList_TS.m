function [sorted_indices, criteria, global_trigger, RESULT] = f_makeReefList_TS(META, RESULT, t, connmetrics, remaining_COTS, thisboatorder, current_COTS_ET, COTS_site_densities)
%March 2023: updated control by Tina to fit new scenarios and changes to control strategy
%Only site-trigger for GBR, North-Central, Central, and Central-South New
%regions updated here so follow that.
%Applying new Habitat maps and re-defined GBR regions

% Updates from Tina (July 2023) - now have fixed target reef list. Control at T, then P, then N.
load('New_regions_TS_SL.mat') %this has been updated for new GBRMPA 2023 PR list and includes target reefs now
%new_regions=cell2table(nregions, 'VariableNames', {'Index' %'Priority''region'}); %read in data as a table, so just change column names so the same
new_regions = nregions;
new_regions.Properties.VariableNames = {'Index' 'Reef_type' 'region'}; 

criteria=[];
others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
idx=META.cntrl_reefID; %% Reefs within MPA

% Also need a CF with NO number.
% Regional = 1 = GBR; 2 = FN; 3 = FNN; 4 = N; 5 = NC; 6 = C; 7 = CS; 8 = S.
% Outbreak front = OF LAT = 9; OF SEC all reefs = 10; OF SEC PR = 11; 
% Effort sink = high COTS = 12; Effort sink = high COTS with coral cover threshold = 13; high dive total threshold and /cull site in separate function, change in f_runmodel (L1265/1266) + GBR wide strategy 1.
% Protection status = NotBlueZones = 14; BlueZones = 15.
% Dynamic scenario = 18. Dynamic + weighting for GZ = 16, dynamic + weighting for BZ = 17.

switch META.COTS_reefs2cull_strat
    case 1 %% No regional strategy: GBR-wide
        % YM: tentative changes below to initialise randomised priority and
        % non-priority lists and keep them unchanged for the rest of the simulation (Caro's comparison with CoCoNet).
        % Requires turning 'COTS fixed list' from 0 to 1 in settings_COTS_CONTROL to avoid any further randomisation (thisboatorder below).
        
        %change which file is read below to change focal region for different strategies
        if t==(META.COTS_control_start+1)%load only once to save time
            gbrmpalist=load('GBRMPAprioritylist.mat'); %Tina - this has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            %RESULT.COTS_priority_list(1).list=priority;
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO

            is_target = ismember(new_regions.Reef_type, 'T'); %Find all target reefs
            target = new_regions.Index(is_target, :);     %Select all target reefs
            RESULT.target_list = target; 
            
            remaining_priority = setdiff(priority, target); %Remove the new target reefs from PR list to see which are remaining
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_remaining_priority = remaining_priority(randperm(length(remaining_priority))); %Randomly shuffle the remaining PR list
            new_priority_list = vertcat(shuffled_target, shuffled_remaining_priority); %Join both together so still have the full list of 500 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
         end    
            
        priority = RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        %NOTE: create a strategy without nonprioty reefs if we exclusively want
        %to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 1 here


    case 2  % Spatial strategy to visit only the FN region
        
        if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'FN'}); %% Find all reefs for Far North
           focus = new_regions(is_focus, :); %Take only ReefIDs from selected region
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all TR
           target = focus(is_target, :);   %Select only TR from selected region
           target = double(table2array(target(:,1)));  %Take TR ReefIDs
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from selected region
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from selected region
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP - actually doesn't return anything so fine to leave for now
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the region
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

       %end of case 2 here
    
    
    case 3  % Spatial strategy to visit only the FN and N region
        
        if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'FN', 'N'}); %% Find all reefs for Far North and North
           focus = new_regions(is_focus, :); %Take only ReefIDs from selected region
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all TR
           target = focus(is_target, :);   %Select only TR from selected region
           target = double(table2array(target(:,1)));  %Take TR ReefIDs
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from selected region
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from selected region
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP - actually doesn't return anything so fine to leave for now
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the region
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

       %end of case 3 here    
    
     case 4  % Spatial strategy to visit only the N region
        
        if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'N'}); %% Find all reefs for North
           focus = new_regions(is_focus, :); %Take only ReefIDs from selected region
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all TR
           target = focus(is_target, :);   %Select only TR from selected region
           target = double(table2array(target(:,1)));  %Take TR ReefIDs
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from selected region
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from selected region
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP - actually doesn't return anything so fine to leave for now
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the region
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

       %end of case 4 here    
    
    
    case 5  % Spatial strategy to visit only the NC region   
    
     % Spatial strategy to visit only the North and Central sectors. 
        
        if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'C', 'N'}); %% Find all reefs for North and Central
           focus = new_regions(is_focus, :); %Take only ReefIDs from North and Central
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all PR
           target = focus(is_target, :);   %Select only PR from North and Central
           target = double(table2array(target(:,1)));  %Take PR ReefIDs: 107 TR in NC
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from North and Central
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs: 41 PR in NC

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from North and Central
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs: 327 non PR reefs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the North and Central region (3331).
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of 148 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

         %end of case 5 here
        
        
 case 6  % Spatial strategy to visit Central
        
        if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'C'}); %% Find all reefs for Central
           focus = new_regions(is_focus, :); %Take only ReefIDs from Central
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all PR
           target = focus(is_target, :);   %Select only PR from Central
           target = double(table2array(target(:,1)));  %Take PR ReefIDs:
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from Central
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs: 

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from Central
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs:
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the North and Central region (3331).
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of 148 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

         %end of case 6 here
        
 case 7  % Spatial strategy to visit only Central & South.

           if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'C', 'S'}); %% Find all reefs for CS
           focus = new_regions(is_focus, :); %Take only ReefIDs from North and Central
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all PR
           target = focus(is_target, :);   %Select only PR from CS
           target = double(table2array(target(:,1)));  %Take PR ReefIDs: 
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from CS
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs:

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from CS
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the North and Central region (3331).
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of 148 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
      
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 7 here

 case 8  % Spatial strategy to visit only South.

           if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = ismember(new_regions.region, {'S'}); %% Find all reefs for South
           focus = new_regions(is_focus, :); %Take only ReefIDs from North and Central
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all PR
           target = focus(is_target, :);   %Select only PR from North and Central
           target = double(table2array(target(:,1)));  %Take PR ReefIDs: 107 TR in NC
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from North and Central
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs: 41 PR in NC

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from North and Central
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs: 327 non PR reefs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the North and Central region (3331).
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of 148 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
      
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 8 here

    case 9  % Outbreak front: GBRMPA strategy that goes to target reefs first, then also goes to 0.5' lat (~50 km) from target reefs with outbreaks - whole GBR.
      if t>=(META.COTS_control_start+1)%load only once to save time
            gbrmpalist=load('GBRMPAprioritylist.mat'); %Tina - this has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            %RESULT.COTS_priority_list(1).list=priority;
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
 %          nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO

            is_target = ismember(new_regions.Reef_type, 'T'); %Find all target reefs
            target = new_regions.Index(is_target, :);     %Select all target reefs
            RESULT.target_list = target; 
            
            remaining_priority = setdiff(priority, target); %Remove the new target reefs from PR list to see which are remaining
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_remaining_priority = remaining_priority(randperm(length(remaining_priority))); %Randomly shuffle the remaining PR list
            new_priority_list = vertcat(shuffled_target, shuffled_remaining_priority); %Join both together so still have the full list of 500 reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
           
           priority = RESULT.COTS_priority_list(1).list;
           nonpriority=RESULT.COTS_nonpriority_list(1).list;
        
     %first, find all target reefs with sites 
        target_outbreaks=[];
        for rfs=1:size(RESULT.target_list,1)
            this_reef=RESULT.target_list(rfs);
            this_reef_sites=COTS_site_densities{this_reef,1};
            this_sites_tow=zeros(size(this_reef_sites,1),1);
            this_sites_over_ET=zeros(size(this_reef_sites,1),1);
            for sts=1:size(this_reef_sites,1)%convert COTS at sites to COTS per tow, also check which ones are over ET
                this_sites_tow(sts,1)=sum(this_reef_sites(sts,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
                if this_sites_tow(sts,1)>current_COTS_ET(this_reef,1)
                    this_sites_over_ET(sts,1)=1;
                end
            end
            sites_over_ET=find(this_sites_over_ET);
            if ~isempty(sites_over_ET)
                target_outbreaks=vertcat(target_outbreaks,this_reef);
            end
        end
        
        %then, find all reefs that are within 0.5' latitude N-S from outbreak target reefs
        %Update-- exclude reefs outside MPA
        META.reef_lat(others)=[];
        
        close_lat=[];
        if ~isempty(target_outbreaks)
            for pob=1:size(target_outbreaks,1)
                this_pri_ob=target_outbreaks(pob,1);
                this_lat=META.reef_lat(this_pri_ob,1);
                % Calculate the absolute latitude difference for all reefs
                abs_lat_diff = abs(this_lat - META.reef_lat);
                % Create a logical index for reefs within 0.5 degrees
                within_range = abs_lat_diff < 0.5;
                 % Use logical indexing to exclude those not within range
                 this_np = find(within_range);
                %this_np=find(abs(this_lat-META.reef_lat(:,1))<0.5);  %Change to 0.5 as 1 degree latitude is 100km, 0.5 is 50 km.
                if ~isempty(this_np)
                   %close_lat= vertcat(close_lat,this_np);
                   close_lat = [close_lat; this_np];
                end
            end
        end

        close_lat=unique(close_lat);%clean it up
        close_lat=setdiff(close_lat,target); %remove the reefs that are already in target list if they were within 1' of another target reef with outbreak
        others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
        close_lat=setdiff(close_lat,others); %%remove reefs outside MPA-- control doesnt operate there
        
        otherreefs=transpose(1:3806);
        pnp=vertcat(target, close_lat);
        otherreefs(pnp)=[];
        otherreefs=setdiff(otherreefs,others);%%remove reefs outside MPA-- control doesnt operate ther
      end

        RESULT.COTS_nonpriority_list(1).list= nonpriority;
        RESULT.COTS_dynamic_nonpriority_list(t).npl=close_lat;
        
        RESULT.COTS_N2consider_reefs=size(target,1)+size(close_lat,1);
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(target,close_lat, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,target)
                    target(find(target==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,target,close_lat, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,close_lat)
                    close_lat(find(close_lat==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,target,close_lat, otherreefs);
                else
                    sorted_indices_all=vertcat(target,close_lat, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(target(randperm(length(target))),close_lat(randperm(length(close_lat))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,target)
                    target(find(target==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,target(randperm(length(target))),close_lat(randperm(length(close_lat))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,close_lat)
                    close_lat(find(close_lat==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,target(randperm(length(target))),close_lat(randperm(length(close_lat))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(target(randperm(length(target))),close_lat(randperm(length(close_lat))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
            
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
  
        %         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

        %%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 9 here

 case 10 % Outbreak front - At each timestep look for the sector (of the 11) that has the highest density of COTS on ALL reefs, start control there, then remaining.
       if t>=(META.COTS_control_start+1)%load only once to save time
            
           numReefs = size(RESULT.COTS_adult_densities, 1); %Number of reefs
           numSectors = max(META.AIMS_sector);  % Number of sectors
           sortedIndices = cell(numSectors, 1); % Initialize a cell array to store the sorted indices for each sector

           % First calculate the overall COTS density for each sector
           sectorCOTSDensity = zeros(numSectors, 1);  
           for sector = 1:numSectors
                sectorIndices = find(META.AIMS_sector == sector);
                sectorCOTSDensity(sector) = sum(RESULT.COTS_adult_densities(sectorIndices, t), 'all');
           end
        
           [~, sortedSectorIndices] = sort(sectorCOTSDensity, 'descend'); % Sort sectors based on overall COTS density
         
            % Then loop through each sorted sector and order reefs
            for i = 1:numSectors
                sector = sortedSectorIndices(i);
                sectorIndices = find(META.AIMS_sector == sector);  % Find indices of reefs in the current sector
                [~, sortedTimestepIndices] = sort(RESULT.COTS_adult_densities(sectorIndices, t), 'descend'); % Sort the indices based on COTS densities for each timestep
                sortedIndices{i} = sectorIndices(sortedTimestepIndices);  % Save the sorted indices for the current sector
            end
        
            allSortedIndices = cat(1, sortedIndices{:});     % Combine the sorted indices for all sectors into a column vector
       end

        RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        priority = allSortedIndices;  %Set sorted indices as priority to keep below code, but really no priority or non-priority.
        nonpriority=transpose(1:3806); 
        nonpriority(priority)=[];  %Again, no P or NP. Just to keep code below.
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list. Keep this as 1.
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
                
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 10 here

case 11 % Outbreak front - At each timestep look for the sector (of the 11) that has the highest proportion of PRIORITY reefs with outbreaks. 
        % Prioritise effort there first, then address the remaining sectors. New order each timestep.

       if t>=(META.COTS_control_start+1)%load only once to save time

            gbrmpalist=load('GBRMPAprioritylist.mat'); % This has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            nonpriority=setdiff(1:3806,priority)'; %remove priority
    
            numReefs = size(RESULT.COTS_adult_densities, 1); % Number of reefs
            numSectors = max(META.AIMS_sector);  % Number of sectors
            sortedIndices = cell(numSectors, 1); % Initialize a cell array to store the sorted indices for each sector
            
            % First calculate the overall COTS density for each sector considering only priority reefs
            sectorCOTSDensity = zeros(numSectors, 1);
        
        for sector = 1:numSectors
            sectorIndices = intersect(find(META.AIMS_sector == sector), priority);
            sectorCOTSDensity(sector) = sum(RESULT.COTS_adult_densities(sectorIndices, t), 'all');
        end
        
        [~, sortedSectorIndices] = sort(sectorCOTSDensity, 'descend'); % Sort sectors based on overall COTS density
        
        % Then loop through each sorted sector and order priority reefs
        for i = 1:numSectors
            sector = sortedSectorIndices(i);
            sectorIndices = intersect(find(META.AIMS_sector == sector), priority); % Find indices of priority reefs in the current sector
            [~, sortedTimestepIndices] = sort(RESULT.COTS_adult_densities(sectorIndices, t), 'descend'); % Sort the indices based on COTS densities for each timestep
            sortedIndices{i} = sectorIndices(sortedTimestepIndices);  % Save the sorted indices for the current sector
        end
        
        allSortedIndices = cat(1, sortedIndices{:}); % Combine the sorted indices for all sectors into a column vector
      
       end

        RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        priority = allSortedIndices;  %Priority now ordered by descending COTS density, first by sectors, then reefs within sectors. Only P, not NP.
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list. Keep this as 1.
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
                
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

        %end of case 11 here

case 12  % Effort sink - At each timestep, still make priority reef list same as before but don't include reefs that have COTS > 3 per tow at that timestep. 
        if t>=(META.COTS_control_start+1) %load 
           gbrmpalist=load('GBRMPAprioritylist.mat'); % This has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO

            is_target = ismember(new_regions.Reef_type, 'T'); %Find all target reefs
            target = new_regions.Index(is_target, :);     %Select all target reefs
            RESULT.target_list = target; 
            
            remaining_priority = setdiff(priority, target); %Remove the new target reefs from PR list to see which are remaining
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_remaining_priority = remaining_priority(randperm(length(remaining_priority))); %Randomly shuffle the remaining PR list
            new_priority_list = vertcat(shuffled_target, shuffled_remaining_priority); %Join both together so still have the full list of 500 reefs
            priority = new_priority_list; % Keep same naming for continuity

            %Tina: for all reefs, next want to remove all > 1.1 adult COTS per grid which is >3 per tow.
            COTS_mt = (0.22/0.6)*RESULT.COTS_total_perceived_density(priority, t); %Convert to COTS per tow and find all for priority at that ts
            priority_indices = find(COTS_mt <= META.max_COTS);  % Get the indices of reefs with COTS < META.max_COTS
                if isempty(priority_indices)==0
                    priority_filtered = priority(priority_indices); % Remove reefs with COTS META.max_COTS from the priority list
                    else priority_filtered = priority;
                end
            COTS_rm_p = setdiff(priority, priority_filtered);
            priority = priority_filtered; % Keep same naming for continuity
            RESULT.COTS_priority_list = priority; % Update the priority list with the new list

            %Tina: for all non-priority reefs, next want to remove all > 1.1 adult COTS per grid which is >3 per tow.
            COTS_mt_n = (0.22/0.6)*RESULT.COTS_total_perceived_density(nonpriority, t); %Convert to COTS per tow and find all for priority at that ts
            nonpriority_indices = find(COTS_mt_n <= META.max_COTS);  % Get the indices of reefs with COTS < META.max_COTS
                if isempty(nonpriority_indices)==0
                    nonpriority_filtered = nonpriority(nonpriority_indices); % Remove reefs with coral cover <= min_cover from the priority list
                    else nonpriority_filtered = nonpriority;
                 end
           COTS_rm_np = setdiff(nonpriority, nonpriority_filtered);
           nonpriority = nonpriority_filtered; % Keep same naming for continuity
           RESULT.COTS_nonpriority_list = nonpriority;

           COTS_rm = vertcat(COTS_rm_p, COTS_rm_np);
           RESULT.COTS_rm_list(t).high_effort_reefs = COTS_rm; % Store target list
       end   

        priority=RESULT.COTS_priority_list;
        nonpriority=RESULT.COTS_nonpriority_list;
        %otherreefs=RESULT.COTS_otherreefs_list(1).list;
        %RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        %NOTE: create a strategy without nonprioty reefs if we exclusively want to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else% here removed randommness and keep on with list which might now contain different reefs if some removed because of high COTS
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority,nonpriority);
                %sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
                %sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
               
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 12 here
 
case 13  % Effort sink with coral cover minimum- At each timestep, still make priority reef list same as before but don't include reefs that have COTS > 3 per tow at that timestep unless they have > threshold coral cover.
        if t>=(META.COTS_control_start+1) %load 
            gbrmpalist=load('GBRMPAprioritylist.mat'); % This has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO

            is_target = ismember(new_regions.Reef_type, 'T'); %Find all target reefs
            target = new_regions.Index(is_target, :);     %Select all target reefs
            RESULT.target_list = target; 
            
            remaining_priority = setdiff(priority, target); %Remove the new target reefs from PR list to see which are remaining
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_remaining_priority = remaining_priority(randperm(length(remaining_priority))); %Randomly shuffle the remaining PR list
            new_priority_list = vertcat(shuffled_target, shuffled_remaining_priority); %Join both together so still have the full list of 500 reefs
            priority = new_priority_list; % Keep same naming for continuity

            %Tina: for all reefs, next want to remove all > 1.1 adult COTS 
            %per grid which is >3 per tow UNLESS coral cover is > 10%? So
            %COTS <= COTS_mt META.max_COTS & > META.min_control_cover
            priority_cover = sum(squeeze(RESULT.coral_pct2D(priority, t, :)), 2); %Find coral cover at current timestep
            COTS_mt = (0.22/0.6)*RESULT.COTS_total_perceived_density(priority, t); %Convert to COTS per tow and find all for priority at that ts

            %priority_indices = find(priority_cover > META.min_control_cover); % Get the indices for minimum coral cover, as here, don't care if high COTS when coral cover high enough. BUT this causes problems in future when coral cver not enough anymore

            priority_indices = find((COTS_mt <= META.max_COTS) | (COTS_mt > META.max_COTS & priority_cover > META.min_control_cover));

            if isempty(priority_indices)==0
                 priority_filtered = priority(priority_indices); % Remove reefs with COTS META.max_COTS from the priority list
                 else priority_filtered = priority;
            end

            COTS_rm_p = setdiff(priority, priority_filtered);
            priority = priority_filtered; % Keep same naming for continuity
            RESULT.COTS_priority_list = priority; % Update the priority list with the new list

            %Tina: for all non-priority reefs, next want to remove all > 1.1 adult COTS per grid which is >3 per tow.
            nonpriority_cover = sum(squeeze(RESULT.coral_pct2D(nonpriority, t, :)), 2); %Find coral cover at current timestep
            COTS_mt_n = (0.22/0.6)*RESULT.COTS_total_perceived_density(nonpriority, t); %Convert to COTS per tow and find all for priority at that ts
           
            % nonpriority_indices = find(nonpriority_cover > META.min_control_cover);  % Get the indices of reefs with COTS < META.max_COTS
            
            nonpriority_indices = find((COTS_mt_n <= META.max_COTS) | (COTS_mt_n > META.max_COTS & nonpriority_cover > META.min_control_cover));

            if isempty(nonpriority_indices)==0
                nonpriority_filtered = nonpriority(nonpriority_indices); % Remove reefs with coral cover <= min_cover from the priority list
                else nonpriority_filtered = nonpriority;
             end

           COTS_rm_np = setdiff(nonpriority, nonpriority_filtered);
           nonpriority = nonpriority_filtered; % Keep same naming for continuity
           RESULT.COTS_nonpriority_list = nonpriority;

           COTS_rm = vertcat(COTS_rm_p, COTS_rm_np);
           RESULT.COTS_rm_list(t).high_effort_reefs = COTS_rm; % Store target list
           end   

        priority=RESULT.COTS_priority_list;
        nonpriority=RESULT.COTS_nonpriority_list;
        %otherreefs=RESULT.COTS_otherreefs_list(1).list;
        %RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        %NOTE: create a strategy without nonprioty reefs if we exclusively want to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else% here removed randommness and keep on with list which might now contain different reefs if some removed because of high COTS
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority,nonpriority);
                %sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
                %sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
               
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 13 here

   case 14  % Protection status - At each timestep, still make priority reef list etc but only include green zone reefs.
        
      if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = (META.GreenZone == 1); % Find all reefs that are green zones
           focus = new_regions(is_focus, :); % Take only rows where GreenZone is 1
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all TR
           target = focus(is_target, :);   %Select only TR from selected region
           target = double(table2array(target(:,1)));  %Take TR ReefIDs
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from selected region
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from selected region
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP - actually doesn't return anything so fine to leave for now
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the region
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

       %end of case 14 here

   case 15 % Protection status - At each timestep, still make priority reef list etc but don't include any green zone reefs.

           if t==(META.COTS_control_start+1)%load only once to save time
           is_focus = (META.GreenZone == 0); % Find all reefs that are green zones
           focus = new_regions(is_focus, :); % Take only rows where GreenZone is 0
           
           is_target = ismember(focus.Reef_type, 'T'); %Find all TR
           target = focus(is_target, :);   %Select only TR from selected region
           target = double(table2array(target(:,1)));  %Take TR ReefIDs
           
           is_priority = ismember(focus.Reef_type, 'P'); %Find all PR
           priority = focus(is_priority, :);   %Select only PR from selected region
           priority = double(table2array(priority(:,1)));  %Take PR ReefIDs

           is_nonpriority = ismember(focus.Reef_type, 'N'); %Find all non PR
           nonpriority = focus(is_nonpriority, :);     %Select all non PR from selected region
           nonpriority = double(table2array(nonpriority(:,1)));  %Take non PR ReefIDs
                           
           others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP - actually doesn't return anything so fine to leave for now
           otherreefs=transpose(1:3806);
           otherreefs(vertcat(target, priority, nonpriority))=[]; 
           otherreefs=setdiff(otherreefs,others); %%These are all the reefs not in the region
            
            shuffled_target = target(randperm(length(target))); %Randomly shuffle the target list
            shuffled_priority = priority(randperm(length(priority))); %Randomly shuffle the PR list
            new_priority_list = vertcat(shuffled_target, shuffled_priority); %Join both together so still have the full list of reefs
            priority = new_priority_list; % Keep same naming for continuity
            RESULT.COTS_priority_list(1).list = priority; % Update the priority list with the new list
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
         end    

        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
       
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

       %end of case 15 here


 case 16 %% GreenZone weighting with CoTS connec and coral cover. 
            % i) At each timestep, focus on PR still. First, remove reefs with < 20 coral cover. 
            % ii) For remaining reefs, give score based on 1) coral, and 2) connec. This gives general order. 
            % iii) Where similar or in same range, weight to green zone. Next scenario will be blue zone first. 

           if t >= (META.COTS_control_start + 1); % 
            gbrmpalist=load('GBRMPAprioritylist.mat'); %Tina - this has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            %RESULT.COTS_priority_list(1).list=priority;
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
            
            %Going to have to do this whole process with priority reefs first, then NP. 
            %First, exclude any priority with < META.MIN Control cover. 
            p_coral_cover = sum(squeeze(RESULT.coral_pct2D(priority, t, :)), 2); %Get CC at PR at that timestep
            p_valid_indices = find(p_coral_cover > META.min_control_cover); %Find indices with coral cover over minimum threshold
           
            %Remove reefs with > minimum coral cover 
            if isempty(p_valid_indices)== 0
                       valid_priority = priority(p_valid_indices); 
            else      %if no reefs > min coral cover                      
                       valid_priority = priority;  %keep whole list
            end
         
            % Get larval output at the time step before. Larvae only in summer, so combine summer and winter values
            p_larval_output = sum(RESULT.COTS_larval_output(p_valid_indices, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            p_normalized_coral_cover = zscore(p_coral_cover(p_valid_indices));
            p_normalized_larval_output = zscore(p_larval_output);
            
            % Calculate the combined score
            p_score = p_normalized_coral_cover + p_normalized_larval_output;
            
            % Prioritize reefs with META.GreenZone == 1 in case of ties
            green_zones = META.GreenZone(p_valid_indices);
            green_zone_indices = find(green_zones == 1);

            % Check if there are valid green zone indices
            if ~isempty(green_zone_indices)
                % Increment score for reefs in green zone
                p_score(green_zone_indices) = p_score(green_zone_indices) + 0.01;
            end
            
            % Now sort the valid reefs based on the combined score
            [sorted_score, sorted_indices] = sort(p_score, 'descend');
            p_sorted_reefs = valid_priority(sorted_indices);

            %Now need to add in priority reefs removed for not having enough CC. Q: do we prefer to control those next, or non PR
            %with > minimum coral? For now, lets say all PR first.

            % Find the remaining priority reefs
            remaining_priority = setdiff(priority, valid_priority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_priority) == 0
                % Calculate the score for the remaining reefs only on COTS larvae, as we know not enough CC
                larval_output_remaining = sum(RESULT.COTS_larval_output(remaining_priority, (t-2):(t-1)),2);
                normalized_larval_output_remaining = zscore(larval_output_remaining);
                score_remaining = normalized_larval_output_remaining;
                
                % Prioritize reefs with META.GreenZone == 1 in case of ties
                green_zones = META.GreenZone(remaining_priority);
                green_zone_indices = find(green_zones == 1);
                score_remaining(green_zone_indices) = score_remaining(green_zone_indices) + 0.01;
            
                % Sort the remaining reefs based on the larval output
                [sorted_score_remaining, sorted_indices_remaining] = sort(score_remaining, 'descend');
                sorted_reefs_remaining = remaining_priority(sorted_indices_remaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_priority_reefs = [p_sorted_reefs; sorted_reefs_remaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_priority_reefs = p_sorted_reefs;
            end

            priority = sorted_priority_reefs; % Keep same naming for continuity
            RESULT.COTS_priority_list = priority; % Update the priority list with the new list

            %Now NonPR reefs, same process for all.
            np_coral_cover = sum(squeeze(RESULT.coral_pct2D(nonpriority, t, :)), 2); %Get CC at NPR
            np_valid_indices = find(np_coral_cover > META.min_control_cover);
           
            if isempty(np_valid_indices)== 0
                       valid_nonpriority = nonpriority(np_valid_indices); % Remove reefs with < MIN coral cover 
            else      %if not reefs > min coral cover                      
                       valid_nonpriority = nonpriority;  %keep whole list
            end
         
            % Get larval output at the time step before, tho here 2 before as no data. Check with YM about timings of steps.
            np_larval_output = sum(RESULT.COTS_larval_output(np_valid_indices, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            np_normalized_coral_cover = zscore(np_coral_cover(np_valid_indices));
            np_normalized_larval_output = zscore(np_larval_output);
            
            % Calculate the combined score
            np_score = np_normalized_coral_cover + np_normalized_larval_output;
            
            % Prioritize reefs with META.GreenZone == 1 in case of ties
            green_zones = META.GreenZone(valid_nonpriority);
            green_zone_indices = find(green_zones == 1);

            % Check if there are valid green zone indices
            if ~isempty(green_zone_indices)
                % Increment score for reefs in green zone
                np_score(green_zone_indices) = np_score(green_zone_indices) + 0.01;
            end
            
            % Now sort the valid reefs based on the combined score
            [sorted_score, sorted_indices] = sort(np_score, 'descend');
            np_sorted_reefs = valid_nonpriority(sorted_indices);

            %Now need to add in nonpriority reefs removed for not having enough CC.

            % Find the remaining nonpriority reefs
            remaining_nonpriority = setdiff(nonpriority, valid_nonpriority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_nonpriority) == 0
                % Calculate the score for the remaining reefs only 
                non_larval_output_remaining = sum(RESULT.COTS_larval_output(remaining_nonpriority, (t-2):(t-1)),2);
                non_normalized_larval_output_remaining = zscore(non_larval_output_remaining);
                non_score_remaining = non_normalized_larval_output_remaining;
                
                % Prioritize reefs with META.GreenZone == 1 in case of ties
                green_zones = META.GreenZone(remaining_nonpriority);
                green_zone_indices = find(green_zones == 1);
                non_score_remaining(green_zone_indices) = non_score_remaining(green_zone_indices) + 0.01;
            
                % Sort the remaining reefs based on the larval output
                [sorted_score_remaining, sorted_indices_remaining] = sort(non_score_remaining, 'descend');
                non_sorted_reefs_remaining = remaining_nonpriority(sorted_indices_remaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_nonpriority_reefs = [np_sorted_reefs; non_sorted_reefs_remaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_nonpriority_reefs = np_sorted_reefs;
            end

            nonpriority = sorted_nonpriority_reefs; % Keep same naming for continuity
            RESULT.COTS_nonpriority_list = nonpriority; % Update the priority list with the new list

        end  
  
        priority=RESULT.COTS_priority_list;
        nonpriority=RESULT.COTS_nonpriority_list;
        %otherreefs=RESULT.COTS_otherreefs_list(1).list;
        %RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        %NOTE: create a strategy without nonprioty reefs if we exclusively want to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else% here removed randommness and keep on with list which might now contain different reefs if some removed because of high COTS
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority,nonpriority);
                %sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
                %sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
               
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

%end of case 16 here


    case 17  % GreenZone weighting with CoTS connec and coral cover, same as case 16 above. Except this time, preferential weighting to bluezones.
            if t >= (META.COTS_control_start + 1); % 
            gbrmpalist=load('GBRMPAprioritylist.mat'); %Tina - this has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            %RESULT.COTS_priority_list(1).list=priority;
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
            
            %Going to have to do this whole process with priority reefs first, then NP. 
            %First, exclude any priority with < META.MIN Control cover. 
            p_coral_cover = sum(squeeze(RESULT.coral_pct2D(priority, t, :)), 2); %Get CC at PR
            p_valid_indices = find(p_coral_cover > META.min_control_cover);
           
            if isempty(p_valid_indices)== 0
                       valid_priority = priority(p_valid_indices); % Remove reefs with < MIN coral cover 
            else      %if not reefs > min coral cover                      
                       valid_priority = priority;  %keep whole list
            end
         
            % Get larval output at the time step before, tho here 2 before as no data. Check with YM about timings of steps.
            p_larval_output = sum(RESULT.COTS_larval_output(p_valid_indices, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            p_normalized_coral_cover = zscore(p_coral_cover(p_valid_indices));
            p_normalized_larval_output = zscore(p_larval_output);
            
            % Calculate the combined score
            p_score = p_normalized_coral_cover + p_normalized_larval_output;
            
            % Prioritize reefs with META.GreenZone == 0 in case of ties
            blue_zones = META.GreenZone(p_valid_indices);
            blue_zone_indices = find(blue_zones == 0);
                        
            if ~isempty(blue_zone_indices)
                % Increment score for reefs in blue zone
                p_score(blue_zone_indices) = p_score(blue_zone_indices) + 0.01;
            end

            % Now sort the valid reefs based on the combined score
            [sorted_score, sorted_indices] = sort(p_score, 'descend');
            p_sorted_reefs = valid_priority(sorted_indices);

            %Now need to add in priority reefs removed for not having enough CC. Q: do we prefer to control those next, or non PR
            %with > minimum coral? For now, lets say all PR first.

            % Find the remaining priority reefs
            remaining_priority = setdiff(priority, valid_priority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_priority) == 0
                % Calculate the score for the remaining reefs only on COTS
                larval_output_remaining = sum(RESULT.COTS_larval_output(remaining_priority, (t-2):(t-1)),2);
                normalized_larval_output_remaining = zscore(larval_output_remaining);
                score_remaining = normalized_larval_output_remaining;
                
                % Prioritize reefs with META.GreenZone == 0 in case of ties
                blue_zones = META.GreenZone(remaining_priority);
                blue_zone_indices = find(blue_zones == 0);
                score_remaining(blue_zone_indices) = score_remaining(blue_zone_indices) + 0.01;
            
                % Sort the remaining reefs based on the larval output
                [sorted_score_remaining, sorted_indices_remaining] = sort(score_remaining, 'descend');
                sorted_reefs_remaining = remaining_priority(sorted_indices_remaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_priority_reefs = [p_sorted_reefs; sorted_reefs_remaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_priority_reefs = p_sorted_reefs;
            end

            priority = sorted_priority_reefs; % Keep same naming for continuity
            RESULT.COTS_priority_list = priority; % Update the priority list with the new list

            %Now NonPR reefs, same process for all.
            np_coral_cover = sum(squeeze(RESULT.coral_pct2D(nonpriority, t, :)), 2); %Get CC at NPR
            np_valid_indices = find(np_coral_cover > META.min_control_cover);
           
            if isempty(np_valid_indices)== 0
                       valid_nonpriority = nonpriority(np_valid_indices); % Remove reefs with < MIN coral cover 
            else      %if not reefs > min coral cover                      
                       valid_nonpriority = nonpriority;  %keep whole list
            end
         
            % Get larval output at the time step before, tho here 2 before as no data. Check with YM about timings of steps.
            np_larval_output = sum(RESULT.COTS_larval_output(np_valid_indices, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            np_normalized_coral_cover = zscore(np_coral_cover(np_valid_indices));
            np_normalized_larval_output = zscore(np_larval_output);
            
            % Calculate the combined score
            np_score = np_normalized_coral_cover + np_normalized_larval_output;
            
            % Prioritize reefs with META.GreenZone == 0 in case of ties
            blue_zones = META.GreenZone(valid_nonpriority);
            blue_zone_indices = find(blue_zones == 0);

            % Check if there are valid blue zone indices
            if ~isempty(blue_zone_indices)
                % Increment score for reefs in blue zone
                np_score(blue_zone_indices) = np_score(blue_zone_indices) + 0.01;
            end
            
            % Now sort the valid reefs based on the combined score
            [sorted_score, sorted_indices] = sort(np_score, 'descend');
            np_sorted_reefs = valid_nonpriority(sorted_indices);

            %Now need to add in nonpriority reefs removed for not having enough CC.

            % Find the remaining nonpriority reefs
            remaining_nonpriority = setdiff(nonpriority, valid_nonpriority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_nonpriority) == 0
                % Calculate the score for the remaining reefs only 
                non_larval_output_remaining = sum(RESULT.COTS_larval_output(remaining_nonpriority, (t-2):(t-1)),2);
                non_normalized_larval_output_remaining = zscore(non_larval_output_remaining);
                non_score_remaining = non_normalized_larval_output_remaining;
                
                % Prioritize reefs with META.GreenZone == 0 in case of ties
                blue_zones = META.GreenZone(remaining_nonpriority);
                blue_zone_indices = find(blue_zones == 0);
                non_score_remaining(blue_zone_indices) = non_score_remaining(blue_zone_indices) + 0.01;
            
                % Sort the remaining reefs based on the larval output
                [sorted_score_remaining, sorted_indices_remaining] = sort(non_score_remaining, 'descend');
                non_sorted_reefs_remaining = remaining_nonpriority(sorted_indices_remaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_nonpriority_reefs = [np_sorted_reefs; non_sorted_reefs_remaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_nonpriority_reefs = np_sorted_reefs;
            end

            nonpriority = sorted_nonpriority_reefs; % Keep same naming for continuity
            RESULT.COTS_nonpriority_list = nonpriority; % Update the priority list with the new list

        end  

        priority=RESULT.COTS_priority_list;
        nonpriority=RESULT.COTS_nonpriority_list;
        %otherreefs=RESULT.COTS_otherreefs_list(1).list;
        %RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        %NOTE: create a strategy without nonprioty reefs if we exclusively want to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else% here removed randommness and keep on with list which might now contain different reefs if some removed because of high COTS
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority,nonpriority);
                %sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
                %sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
               
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;

%end of case 17 here     


 case 18  % CoTS connec and coral cover, same as case 16 and 17 above, but no preferential weighting to blue or green zones.         
         if t >= (META.COTS_control_start + 1); % 
            gbrmpalist=load('GBRMPAprioritylist.mat'); %Tina - this has been updated for new GBRMPA 2023 PR list
            priority=gbrmpalist.GBRMPAlist; %Tina - changed from csirolist to gbrmpalist
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            %RESULT.COTS_priority_list(1).list=priority;
                
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%           nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
            
            %Going to have to do this whole process with priority reefs first, then NP. 
            %First, exclude any priority with < META.MIN Control cover. 
            p_coral_cover = sum(squeeze(RESULT.coral_pct2D(priority, t, :)), 2); %Get CC at PR at that timestep
            p_valid_indices = find(p_coral_cover > META.min_control_cover); %Find indices with coral cover over minimum threshold
           
            %Remove reefs with < minimum coral cover 
            if isempty(p_valid_indices)== 0
                       valid_priority = priority(p_valid_indices); 
            else      %if no reefs > min coral cover                      
                       valid_priority = priority;  %keep whole list
            end
         
            % Get larval output at the time step before. Larvae only in summer, so combine summer and winter values
            p_larval_output = sum(RESULT.COTS_larval_output(valid_priority, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            if ~isempty(p_valid_indices)
                p_normalized_coral_cover = zscore(p_coral_cover(p_valid_indices));
            else
                % Set p_normalized_coral_cover to a matrix of zeros with the same size as p_coral_cover
                p_normalized_coral_cover = zeros(size(p_coral_cover));
            end

            p_normalized_larval_output = zscore(p_larval_output);
       
            % Calculate the combined score
            p_score = p_normalized_coral_cover + p_normalized_larval_output;
                      
            % Now sort the valid reefs based on the combined score
            [sorted_score, sorted_indices] = sort(p_score, 'descend');
            p_sorted_reefs = valid_priority(sorted_indices);

            %Now need to add in priority reefs removed for not having enough CC. Q: do we prefer to control those next, or non PR
            %with > minimum coral? For now, lets say all PR first.

            % Find the remaining priority reefs
            remaining_priority = setdiff(priority, valid_priority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_priority) == 0
                % Calculate the score for the remaining reefs only on COTS larvae, as we know not enough CC
                larval_output_remaining = sum(RESULT.COTS_larval_output(remaining_priority, (t-2):(t-1)),2);
                normalized_larval_output_remaining = zscore(larval_output_remaining);
                score_remaining = normalized_larval_output_remaining;
                           
                % Sort the remaining reefs based on the larval output
                [sorted_score_remaining, sorted_indices_remaining] = sort(score_remaining, 'descend');
                sorted_reefs_remaining = remaining_priority(sorted_indices_remaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_priority_reefs = [p_sorted_reefs; sorted_reefs_remaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_priority_reefs = p_sorted_reefs;
            end

            priority = sorted_priority_reefs; % Keep same naming for continuity
            RESULT.COTS_priority_list = priority; % Update the priority list with the new list
           

            %Now NonPR reefs, same process for all.
            np_coral_cover = sum(squeeze(RESULT.coral_pct2D(nonpriority, t, :)), 2); %Get CC at NPR at that timestep
            np_valid_indices = find(np_coral_cover > META.min_control_cover); %Find indices with coral cover over minimum threshold
           
            %Remove reefs with < minimum coral cover 
            if isempty(np_valid_indices)== 0
                       valid_nonpriority = nonpriority(np_valid_indices); 
            else      %if no reefs > min coral cover                      
                       valid_nonpriority = nonpriority;  %keep whole list
            end
         
            % Get larval output at the time step before. Larvae only in summer, so combine summer and winter values
            np_larval_output = sum(RESULT.COTS_larval_output(valid_nonpriority, (t-2):(t-1)),2);
            
            % Standardize both coral cover and larval output
            if ~isempty(np_valid_indices)
                np_normalized_coral_cover = zscore(np_coral_cover(np_valid_indices));
            else
                % Set p_normalized_coral_cover to a matrix of zeros with the same size as p_coral_cover
                np_normalized_coral_cover = zeros(size(np_coral_cover));
            end

            np_normalized_larval_output = zscore(np_larval_output);
       
            % Calculate the combined score
            np_score = np_normalized_coral_cover + np_normalized_larval_output;
                      
            % Now sort the valid reefs based on the combined score
            [sorted_score_non, sorted_indices_non] = sort(np_score, 'descend');
            np_sorted_reefs = valid_nonpriority(sorted_indices_non);

            %Now need to add in nonpriority reefs removed for not having enough CC.

            % Find the remaining priority reefs
            remaining_nonpriority = setdiff(nonpriority, valid_nonpriority);
            
            % Check if there are remaining priority reefs
            if isempty(remaining_nonpriority) == 0
                % Calculate the score for the remaining reefs only on COTS larvae, as we know not enough CC
                larval_output_nonremaining = sum(RESULT.COTS_larval_output(remaining_nonpriority, (t-2):(t-1)),2);
                normalized_larval_output_nonremaining = zscore(larval_output_nonremaining);
                score_nonremaining = normalized_larval_output_nonremaining;
                           
                % Sort the remaining reefs based on the larval output
                [sorted_score_nonremaining, sorted_indices_nonremaining] = sort(score_nonremaining, 'descend');
                sorted_reefs_nonremaining = remaining_nonpriority(sorted_indices_nonremaining);
            
                % Combine the sorted valid priority reefs and remaining priority reefs
                sorted_nonpriority_reefs = [np_sorted_reefs; sorted_reefs_nonremaining];
            else
                % If there are no remaining priority reefs, set sorted_priority_reefs to sorted valid reefs
                sorted_nonpriority_reefs = np_sorted_reefs;
            end

            nonpriority = sorted_nonpriority_reefs; % Keep same naming for continuity
            RESULT.COTS_nonpriority_list = nonpriority; % Update the priority list with the new list

       end  
  
        priority=RESULT.COTS_priority_list;
        nonpriority=RESULT.COTS_nonpriority_list;
        %otherreefs=RESULT.COTS_otherreefs_list(1).list;
        %RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        %NOTE: create a strategy without nonprioty reefs if we exclusively want to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef visited on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else% here removed randommness and keep on with list which might now contain different reefs if some removed because of high COTS
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority,nonpriority);
                %sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
                %sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
               
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
        end
        
        %criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        criteria1=zeros(size(sorted_indices_all,1),5);%change for less storage 
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        %criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        %criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        %criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;


          
%end of case 18 here  
             

end