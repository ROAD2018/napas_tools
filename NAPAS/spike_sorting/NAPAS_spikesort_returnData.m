% 2019-09-15 Dylan Royston
%
% From 2018-04-10 Dylan Royston
%
% Function that initializes and returns PCA-based spike sorting
% Copied to clean up and optimize for Somatomapping
%
% === INPUTS ===
% DataTime:     duration of data (in minutes); this is used for something?
% numofsess:    number of sessions (days) to be sorted at a time (this will be 1 for SM data, could be removed later to optimize)
%
%
%
% === UPDATES ===
%
%
%
%%

function [exp_snip, exp_time, Data, idata] = FUNC_SpikeSort_returnData(subject, DataTime, numofsess, Data, idata, channel_choice, color_flag)

switch subject
    case 'subj1'
        max_channels = 256;
end

% assign default values if not provided
if nargin < 2 || isempty(numofsess)
    numofsess = 1;
end

if nargin < 1 || isempty(DataTime)
    DataTime = 3;
end

if isempty(channel_choice)
    channel_choice = 1:max_channels;
end

InitParams.mu = [-100 100 0 50 -50;100 -100 0 50 -50];
InitParams.Sigma = [100,0;0,100];

DAT = cell(numofsess,1);
for iSess=1:numofsess
    
    % if data files aren't provided, prompt for manual load; otherwise, load data
    if ~exist('Data')
        [Data, idata] = prepData();
    else
        DAT{iSess} = {idata};
    end% IF, ~exist('Data')
end

exp_snip = cell(numofsess,max_channels,5);
exp_time = cell(numofsess,max_channels,5);

total_channels = length(channel_choice);

% progress bar graphic
H = waitbar(1/total_channels,['Channel 1/' num2str(total_channels)]);
set(H, 'Position', [1988 64 270 56]);

for sess = 1:numofsess
    
    % load data for each session
    idata = DAT{sess}{1};
    
    iChannel = 0;
    lastSortedChannel = 0;
    
    %% while loop containing sorting activity

colorspec = jet(48);
% colorspec = jet(5);

    while iChannel < total_channels

        
        iChannel =      iChannel+1;
        channel_num =   channel_choice(iChannel);
        
        % load snips, spiketimes (loading format depends on pedestal)
        % to undo edits replace channel_num with iChannel
        if channel_num<129
            snips = idata.QL.Data.SPIKE_SNIPPET.ss.snippet(:,idata.QL.Data.SPIKE_SNIPPET.ss.source_index==0&...
                idata.QL.Data.SPIKE_SNIPPET.ss.channel==channel_num);
            spike_times = idata.QL.Data.SPIKE_SNIPPET.ss.source_timestamp(idata.QL.Data.SPIKE_SNIPPET.ss.source_index==0&...
                idata.QL.Data.SPIKE_SNIPPET.ss.channel==channel_num);
        else
            snips = idata.QL.Data.SPIKE_SNIPPET.ss.snippet(:,idata.QL.Data.SPIKE_SNIPPET.ss.source_index==1&...
                idata.QL.Data.SPIKE_SNIPPET.ss.channel==channel_num-128);
            spike_times = idata.QL.Data.SPIKE_SNIPPET.ss.source_timestamp(idata.QL.Data.SPIKE_SNIPPET.ss.source_index==1&...
                idata.QL.Data.SPIKE_SNIPPET.ss.channel==channel_num-128);
        end
        
        
        % 2018-08-09 Royston: adding outlier removal to keep plots from getting fucked
        snip_maxes =            max(snips);
        outliers =              find(snip_maxes > 4*median(snip_maxes) );
        snips(:, outliers) =    [];
        spike_times(outliers) = [];
        
        snip_mins =             min(snips);
        outliers =              find(snip_mins < 4*median(snip_mins) );
        snips(:, outliers) =    [];
        spike_times(outliers) = [];
        
        
        if color_flag == 1
            
            [peak_vals, peak_loc] = min(snips);
            peak_vals = (-1 * peak_vals) + max(snips);
            scaled_peaks = floor( ( (peak_vals - min(peak_vals) ) / (max(peak_vals) - min(peak_vals) ) ) *48 );
            scaled_peaks(scaled_peaks ==0) = 1;
            
            % 2018-08-14 Royston: added an exception for channes with a single spike, which happens somehow
            if isnan(scaled_peaks)
                scaled_peaks = 1;
            end
            
            % plot snips
            subset_size = 1000;
%             subset_size = size(snips, 2);
            figure(2);
            if size(snips,2)<subset_size
                randSnip = [1:1:size(snips,2)];
                for j = 1:size(snips,2)
                    plot(1:48,snips(:,j), 'Color', colorspec(scaled_peaks(j), :));hold on;
                end
            else% if there are more than 400 snippets, randomly plot some
                randSnip = randperm(size(snips,2));
                for j = 1:subset_size
                    plot(1:48,snips(:,randSnip(j)), 'Color', colorspec(scaled_peaks(randSnip(j)), :));hold on;
                end
            end
            hold off;
            set(gcf, 'Position', [2798 371 804 625]);
            set(gca, 'Color', 'k');
            
        else
            scaled_peaks = [];
            
            % plot snips
            subset_size = 1000;
            figure(2);
            if size(snips,2)<subset_size
                randSnip = [1:1:size(snips,2)];
                for j = 1:size(snips,2)
                    plot(1:48,snips(:,j));hold on;
                end
            else% if there are more than 400 snippets, randomly plot some
                randSnip = randperm(size(snips,2));
                for j = 1:subset_size
                    plot(1:48,snips(:,randSnip(j)));hold on;
                end
            end
            hold off;
            set(gcf, 'Position', [2798 371 804 625]);
            
        end% IF color_flag
        
        % increment progress bar
        waitbar(iChannel/total_channels,H,['Channel ',num2str(iChannel),'/' num2str(total_channels)]);
        
        
        % this is where DataTime takes place, seems to be a frequency limiter? number of samples?
        if size(snips,2)>DataTime*60/20 && size(snips,2)<DataTime*60*100 % low limit set for 1/20 Hz, high limit for 100 Hz

            % working on a try-catch to stop PCA from crashing and dumping all data
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    % call function to run PCA on current snippets and re-plot sorted snippets
                    [mus,sigmas, ppis, K ,sorted] = NAPAS_spikesort_runPCA(InitParams,snips, scaled_peaks, color_flag); % sort
                catch MException
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
            
            
            if K == 2 || K == 3 || K == 4 || K == 5
                figure(2);
                Sort = sorted(randSnip);
                
                % 2018-04-17 Royston: check and reassign unit-ID in descending snippet counts (most=1)
                for c_idx = 1 : K
                    cluster_points{c_idx} =     find(Sort == c_idx);
                end
                num_points =            cell2mat( cellfun(@length, cluster_points, 'uni', false) );
                [~, cluster_order] =    sort(num_points, 'descend');
                new_sorted =            zeros(1, length(Sort));
                for new_c_idx = 1 : K
                    new_sorted(cluster_points{new_c_idx}) = find(cluster_order == new_c_idx);
                end
                Sort = new_sorted;
                
                if color_flag == 1
                    set(gca, 'Color', 'k');
                end
                
                % plot snippets according to cluster
%                 
                for n = 1:min(size(snips,2),400)
                    if Sort(n) == 1;
                        if color_flag == 1
                            plot(1:48,snips(:,randSnip(n)),'w');hold on;
                            set(gca, 'Color', 'k');
                        else
                            plot(1:48,snips(:,randSnip(n)),'k');hold on;
                        end
                    elseif Sort(n) == 2;
                        if color_flag == 1
                            plot(1:48,snips(:,randSnip(n)),'g');hold on;
                            set(gca, 'Color', 'k');
                        else
                            plot(1:48,snips(:,randSnip(n)),'g');hold on;
                        end
                    elseif Sort(n) == 3;
                        plot(1:48,snips(:,randSnip(n)),'b');hold on;
                    elseif Sort(n) == 4;
                        plot(1:48,snips(:,randSnip(n)),'c');hold on;
                    elseif Sort(n) == 5;
                        plot(1:48,snips(:,randSnip(n)),'m');hold on;
                    end
                end
                
                %pause for visual inspection of clusters
                pause(3);
                hold off;
            end% IF, K
            
            if K>1 && K < 9
                for j = 1:K
                    exp_snip(sess,channel_num,j) = {snips(:,sorted' == j)};
                    exp_time(sess,channel_num,j) = {spike_times(sorted' == j)};
                end
            elseif K == 9
                if lastSortedChannel > 0
                    iChannel = lastSortedChannel - 1;
                    exp_snip(sess,lastSortedChannel,:) = cell(1,1,5);
                    exp_time(sess,lastSortedChannel,:) = cell(1,1,5);
                else
                    iChannel = iChannel - 1;
                end% IF, lastSortedChannel > 0
                continue
            else
                exp_snip(sess,channel_num,1) = {snips};
                exp_time(sess,channel_num,1) = {spike_times};
            end% IF, K>1
            lastSortedChannel = iChannel;
        else
            exp_snip(sess,iChannel,1) = {snips};
            exp_time(sess,iChannel,1) = {spike_times};
        end% IF, size(snips, 2)
        
        
    end% WHILE, iChannel
    
    %%
    
end% FOR, sess

close(H);

end% FUNCTION
