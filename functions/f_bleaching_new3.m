% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
%
% INHERITS FROM f_bleaching_new+f_bleaching_new2
% Modified (08/2013) to integrate species-specific bleaching mortalities
% Modified (08/2018) to integrate genetic adaptation to warming
% Modified (06/2023) to integrate inter-individual heat tolerance
% Modified (10/2023) after adding heritability in HT
% Modified (03/2026) for tracking realised mortality per species
% -------------------------------------------------------------------------

function [coral, genes, algal, coral_cover_loss, coral_mortality, total_mortality] = f_bleaching_new3(coral, genes, algal, bleaching_whole_mortality,...
    CORAL, doing_3D, nb_coral_types, doing_clades, doing_genetics, bleaching_whole_offset, bleaching_partial_offset,Topt_baseline, Topt2index)

% Extract data from the structures (need to be filled again when leaving)
algal_cm2 = [algal.cover_cm2] ;
[coral_cm2, coral_HT, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy (coral);

id0 = zeros(size(coral_cm2)) ;

sensitivity_bleaching = id0 ;
extent_bleaching = id0 ;
proba_switching = id0 ;

id1 = id0+1;
id1(coral_cm2 <= 0) = 0 ; % Remove the already dead ones

col_start = 1;
col_stop = 0;

%% Assign sensitivity to bleaching to each coral colony
if doing_genetics == 0  
        
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s)*exp(-0.347*coral_HT(:,col_start:col_stop));    
        extent_bleaching(:,col_start:col_stop)= CORAL.bleaching_partial_extent(s);
        proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
        
        col_start = col_start + species_ID(s) ;     
    end 
    
else % code with genetics needs to be re-visited in light of new HT
    
    MORT = id0 ;
    
    for s = 1:nb_coral_types
        
        if species_ID(s)>0
            
            col_stop = col_stop + species_ID(s) ;
            
            id1_tmp = id1(:,col_start:col_stop);
            s_mort = zeros(size(id1_tmp)) ;
            colony_ID_tmp = colony_ID(:,col_start:col_stop);
            
            % First update the list of QTLs to only keep the survivors from previous mortality
            list_old = genes(s).list_coral_ID;
            list_new = colony_ID_tmp(id1_tmp==1);
            
            check = ismember(list_old,list_new) ;
            
            genes(s).QTLs(check==0,:,:)=[];
            genes(s).list_coral_ID(check==0,:)=[];
            genes(s).phenotypes(check==0,:)=[];
            
            % Phenotype is Topt so need to find out the corresponding DHW
            % Then convert into whole colony mortality using Hughes relationship
            % Output must be a matrix of proba of whole colony-mortality
            All_Topts = round(10*(Topt_baseline + genes(s).phenotypes))/10 ;
            
            if isempty(All_Topts)==0     
               s_mort(id1_tmp==1) = bleaching_whole_mortality(10*All_Topts - Topt2index) ;              
            end
            
            % New (06/2023): heat tolerance effect size for each colony affects their sensitivity coefficient to bleaching
            % Note this variability is independent of coral genotypes
            % sensitivity_bleaching(:,col_start:col_stop) = CORAL.sensitivity_bleaching(s);
            sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s)*exp(-0.347*coral_HT(:,col_start:col_stop));
            extent_bleaching(:,col_start:col_stop) = CORAL.bleaching_partial_extent(s);
            proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
            MORT(:,col_start:col_stop) = s_mort ;
            
            col_start = col_start + species_ID(s) ;
            
        end
    end    
end

% if doing_clades == 1
    % Clade-induced tolerance to thermal stress
    % THIS NEEDS TO BE REVISITED DUE TO NEW IMPLEMENTATION OF CORAL TOLERANCE PER INDIVIDUAL
%     sensitivity_bleaching(clade==2) = sensitivity_bleaching(clade==2) * CORAL.bleaching_tolerance_clade ;
% end

%________________________________
%
% Whole-colony mortality
%________________________________

% id1(coral_cm2 < CORAL.adol_size) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)
id1(coral_cm2 < CORAL.size_threshold_wcm) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)

% Proba of initial mortality is obtained by adjusting the expected mortality for an average reef
% (bleaching_whole_mortality, determined by DHW) to the sensitivity to bleaching of each coral
prob_initial_mortality = sensitivity_bleaching*CORAL.bleaching_depth*bleaching_whole_mortality; 

prob_initial_mortality(prob_initial_mortality>1)=1; % cap to 1 before extrapolating to 6 month
prob_whole_mortality = 1-(1-prob_initial_mortality).^bleaching_whole_offset;

id_dead = id0 ;
id_dead(id1==1) = binornd(1, prob_whole_mortality(id1==1)); %independent drawing for each colony, winners=0, losers=1
coral_loss = coral_cm2.*id_dead ;

algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_loss,2) ;
% coral_cm2(id_dead==1) = - coral_cm2(id_dead==1);  % now dead (negatives)
coral_cm2(id_dead==1) = 0; % we do not track dead colonies (negatives) any more
total_mortality = sum(sum(id_dead))/sum(sum(id1)); % total mortality across groups

%________________________________
%
% Partial mortality
%________________________________

id2 = id1 - id_dead ; % update the survivors
prob_partial_mortality = id2.*(1-(1-prob_initial_mortality).^bleaching_partial_offset); % May 2023 -> use initial mortality which is already capped to 1

id_part = id0 ;
id_part(id2==1) = binornd(1, prob_partial_mortality(id2==1)); %independent drawing for each colony, winners=0, losers=1
bleach_extent = floor(extent_bleaching.* coral_cm2 .* id_part) ; 
algal_cm2(:,1) = algal_cm2(:,1) + sum(bleach_extent, 2) ;
coral_cm2 = coral_cm2 - bleach_extent ;

%________________________________
%
% Clade switching
%________________________________
if doing_clades == 1
    rand_switch = rand(size(id2)) ; % Generates random probability of switching to the thermally-tolerant clade (clade 2)
    clade(rand_switch < proba_switching & coral_cm2>0) = 2 ; % NOTE THIS INDEPENDENT OF BLEACHING MORTALITY
end


%%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, coral_HT, surface_cm2, volume_cm3, colony_ID, clade, species_ID, nb_coral_types, doing_clades, doing_3D);

for a=1:size(algal_cm2, 2) 
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end

%%%%%%%% Keep track of losses per species
count = 1;
coral_cover_loss = zeros(1,nb_coral_types); % loss of coral cover per species (as aresult of whole-colony + partial mortality
coral_mortality = coral_cover_loss; % realised mortality per species (proportion of colonies lost due to bleaching)

for s = 1:nb_coral_types
    select_sp = count:(count + species_ID(s)-1);
    coral_cover_loss(s) = sum(sum(coral_loss(:,select_sp)))+sum(sum(bleach_extent(:,select_sp)));
    coral_mortality(s) = sum(sum(id_dead(:,select_sp)))/sum(sum(id1(:,select_sp)));
    count = count + species_ID(s);
end    
