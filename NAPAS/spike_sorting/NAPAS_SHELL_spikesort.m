% 2019-09-10 Dylan Royston
% 
% 2018-04-10 Dylan Royston
% Shell script for using PCA to sort single-unit spikes
%
% For SM data, an entire session (day) should be sorted at once, which will be slow/data-intensive
% Possible alternate for future update: sort 1 session and just apply that sort to subsequent trials
%
%
%
%
%
%% Initialize and load target data

clear; clc;

subject =               	'subj1';

date_to_analyze =           'date';

raw_data_path = 			'raw_path';
save_data_path = 			'save_path';

switch subject
    case 'subj1'
        directory =             fullfile(save_data_path, date_to_analyze];
        raw_file_directory =    raw_data_path;
end% SWITCH subject


switch subject
    case 'subj1'
        max_channels =          256;
end% SWITCH subject


load_data =     'auto_multi';

switch load_data
    case 'manual'
        [raw_data, raw_idata] = prepData();
        
    case 'auto_multi'
        
        all_filenames =     {};
        
        for f_idx = 1 : length(raw_file_folders)
            
            clearvars files_to_load
            
            folder_string =     [raw_file_directory raw_file_folders{f_idx}];
            contents =          dir(folder_string);
            files =             {contents(:).name};
            bin_files =         FUNC_find_string_in_cell(files, '.bin');
            
            for file_idx = 1 : length(bin_files)
                full_file_string =              [folder_string '\' files{bin_files(file_idx)}];
                files_to_load{file_idx} =       full_file_string;
            end% FOR, file_idx
            
            all_filenames =     [all_filenames files_to_load];
            
        end% FOR, f_idx
        
        [raw_data, raw_idata] = prepData('files', all_filenames);
        
end


cd([directory]);

disp('*** DATA LOADED ***');

%% Perform PCA sorting on sessions

numMinutes =            60;
numofsess =             1;



channel_choice =        1:max_channels;

color_flag =            1;
% perform spike sorting
% this function doesn't actually alter data or idata
tic;
[snippets, timing, new_data, new_idata] =   NAPAS_spikesort_returnData(subject, numMinutes, numofsess, raw_data, raw_idata, channel_choice, color_flag);
toc

snippets =                          squeeze(snippets);
timing =                            squeeze(timing);



% insert edited sorts into full sort data

% this should be spun off into a separate "edit spike sort" function, but it lives here for now
% if only editing select channel sorts, regenerate snippets/timing from existing idata so new sort can be inserted
if length(channel_choice)<max_channels
    
    new_snippets =  cell(max_channels, 5);
    new_timing =    cell(max_channels, 5);
    
    existing_snippet_data =         new_idata.QL.Data.SPIKE_SNIPPET.ss;
    
    for chan_idx = 1 : max_channels
        
        % retain the edited channel data, insert original for everything else
        if isempty(find(channel_choice == chan_idx) )
            
            if chan_idx<(max_channels/2)
                src_idx =                   0;
                chan_idx_temp =             chan_idx;
            else
                src_idx =                   1;
                chan_idx_temp =             chan_idx-128;
            end
            
            correct_source_spikes =             find(existing_snippet_data.source_index == src_idx);
            correct_chan_spikes =               find(existing_snippet_data.channel == chan_idx_temp);
            found_spikes =                      intersect(correct_source_spikes, correct_chan_spikes);
            found_spike_units =                 existing_snippet_data.unit(found_spikes);
            num_units =                         unique(found_spike_units);
            
            for unit_idx = 1 : length(num_units)
                unit_spike_idx =                    found_spikes(found_spike_units == (unit_idx-1) );
                unit_snippets =                     existing_snippet_data.snippet(:, unit_spike_idx );
                unit_times =                        existing_snippet_data.source_timestamp(unit_spike_idx);
                new_snippets(chan_idx, unit_idx) =  {unit_snippets};
                new_timing(chan_idx, unit_idx) =    {unit_times};
            end% FOR, unit_idx
            
        else
            new_snippets(chan_idx, :) =         snippets(chan_idx, :);
            new_timing(chan_idx, :) =           timing(chan_idx, :);
        end% IF, isempty
    end% FOR, chan_idx
else
    new_snippets =      snippets;
    new_timing =        timing;
end% IF, channel_choice

disp('*** STRUCTURES MADE ***');

save([date_to_analyze '_full_sort'], 'raw_data', 'raw_idata', 'new_snippets', 'new_timing', '-v7.3');

disp('*** RAW DATA SAVED ***');

%% edit data structures

%Step 2. Run spikes back into native data structures
% (this function is skipped because its primary value is rebinning spikes into the data.spikecount field, which I don't use, and...
% it doesn't do the activeunitmask correctly, so I'm just editing that here. NOTE: this does mean data.spikecount is unsorted)

remask_method = 'manual';
% remask_method = 'function';


switch remask_method
    case 'function'
        
        tic;
        new_data =          NAPAS_spikesort_makeChanMask(raw_data, squeeze(new_timing), subject);
        toc
        
    case 'manual'
        
        active_units =      cellfun(@isempty, new_timing);
        active_units =      ~active_units;
        active_units =      active_units';
        active_unit_mask =  reshape(active_units, 1, max_channels*5);
        % convert to chan/unit for debugging
        [chans, units] =    unitToChan(find(active_unit_mask));
        
        new_data = raw_data;
        
        new_data.SortedActiveChannelMask = active_unit_mask;
        
end% SWITCH, remask_method

% new_data.ActiveChannelMask =            active_unit_mask;
new_idata =                             raw_idata;
new_idata.QL.Data.SPIKE_SNIPPET.ss =    FUNC_insert_sorted_idata(subject, new_snippets, new_timing);

%Step the next: data is now ready to be run through whatever analysis you already wrote
sorted_task_data.data =       new_data;
sorted_task_data.idata =      new_idata;

disp('*** DATA STRUCTURES ORGANIZED ***');


%% extract template waveforms from each unit

unit_waveform_hold = zeros(48, max_channels, 5);
unit_waveforms =    cell(max_channels, 5);

figure; hold on;

for chan_idx = 1 : max_channels
    
    num_units =         length( find( cellfun(@isempty, new_snippets(chan_idx, :) ) == 0 ) );
    
    for unit_idx = 1 : num_units
        
        unit_snippets =         new_snippets{chan_idx, unit_idx};
        
        unit_template(:, 1) =   mean(unit_snippets, 2);
        unit_template(:, 2) =   std(unit_snippets');
        
        unit_waveforms{chan_idx, unit_idx} = unit_template;
        
        unit_waveform_hold(:, chan_idx, unit_idx) = unit_template(:, 1);
        
        plot(unit_template(:, 1));
        
    end% FOR, unit_idx
    
    
end% FOR, chan_idx


disp('*** WAVEFORMS EXTRACTED ***');


%% filter out non-spike outliers before John's dot-product isolation

scrubbed_unit_waveforms =   unit_waveforms;
scrubbed_snippets =         new_snippets;
scrubbed_isolated_mask =    zeros(max_channels, 5);

figure;
hold on;

fraction_illegal = NaN(size(scrubbed_isolated_mask));

for chan_idx = 1 : max_channels
    
    for unit_idx = 1 : 5
        
        if ~isempty(unit_waveforms{chan_idx, unit_idx})
            
            unit_template_wave = unit_waveform_hold(:, chan_idx, unit_idx);
            
            plot(unit_template_wave, 'k');
            
            % LVL1: spikes should dip negative at ~t=12 (primary peak)
            if (unit_template_wave(12) > -15)
                
                scrubbed_isolated_mask(chan_idx, unit_idx) = -2;
                plot(unit_template_wave, 'Color', 'r');
            else
            end
            
            % LVL2: spikes should start around 0
            if (unit_template_wave(1) > 10)
                scrubbed_isolated_mask(chan_idx, unit_idx) = -2;
                plot(unit_template_wave, 'Color', 'r');
            else
            end
            
            % LVL3: spikes should rise above 0 at t=24 (rebound)
            if (unit_template_wave(24) < 0)
                scrubbed_isolated_mask(chan_idx, unit_idx) = -2;
                plot(unit_template_wave, 'Color', 'r');
            else
            end
            
            % LVL4: spikes should peak (subj max, obj negative) later than t18 (implemented for weird units in Jan's data
            [wave_max, max_loc] = max(unit_template_wave);
            if max_loc == 18
                scrubbed_isolated_mask(chan_idx, unit_idx) = -2;
                plot(unit_template_wave, 'Color', 'r');
            end
            
            %             % 2018-08-21 Royston: LVL5, ISI calc
            current_unit_timing =   new_timing{chan_idx, unit_idx};
            ISI =                   diff(current_unit_timing)*1000;
            ISIs{chan_idx, unit_idx} = ISI;
            fraction_illegal(chan_idx, unit_idx) =      (length(find(ISI<1))/length(ISI) )*100;
            
            if fraction_illegal(chan_idx, unit_idx) > 10
                scrubbed_isolated_mask(chan_idx, unit_idx) = -3;
                plot(unit_template_wave, 'Color', 'r');
            end
            
        end% IF isempty
        
        
    end% FOR, unit_idx
    
end% FOR, chan_idx

disp('*** OUTLIERS MARKED ***');

%% John's unit isolation analysis

exp_snip = squeeze(new_snippets);
exp_time = squeeze(new_timing);

numOfChans = size(exp_snip,1);

meanDotProd = nan(size(exp_snip));

for chan = 1:numOfChans
    for unit = 1:5
        if ~isempty(exp_snip{chan,unit})
            dotProd = nan(1,size(exp_snip{chan,unit},2)-1);
            for i = 1:size(exp_snip{chan,unit},2)-1
                dotProd(i) = dot( double(exp_snip{chan,unit}(:,i)) ./ norm(double(exp_snip{chan,unit}(:,i))),...
                    double(exp_snip{chan,unit}(:,i+1)) ./ norm(double(exp_snip{chan,unit}(:,i+1))));
            end% FOR, snippet
            meanDotProd(chan,unit) = mean(dotProd);
        end% IF, ~isempty
    end% FOR, unit
end% FOR, chan

figure();
h = histogram(meanDotProd(~isnan(meanDotProd)),50,'Normalization','pdf');
x = h.BinEdges(1:50)+.5*h.BinWidth;
y = h.Values;


%%


% 2018-08-15 Royston: originally John used a gauss2, that hasn't been working so I changed to a gauss3 and pick whichever 2 work best
clearvars dist_coeffs gaussFit

% gaussFit = fit(x',y','gauss1');
% dist_coeffs = [gaussFit.a1 gaussFit.b1 gaussFit.c1];
%
% gaussFit = fit(x',y','gauss2');
% dist_coeffs = [gaussFit.a1 gaussFit.b1 gaussFit.c1; gaussFit.a2 gaussFit.b2 gaussFit.c2];
%
gaussFit = fit(x',y','gauss3');
dist_coeffs = [gaussFit.a1 gaussFit.b1 gaussFit.c1; gaussFit.a2 gaussFit.b2 gaussFit.c2; gaussFit.a3 gaussFit.b3 gaussFit.c3];
%
hold on;
dists(1) = plot(x,gaussFit.a1*exp(-((x-gaussFit.b1)/gaussFit.c1).^2), 'Color', 'b');
dists(2) = plot(x,gaussFit.a2*exp(-((x-gaussFit.b2)/gaussFit.c2).^2), 'Color', 'r');
dists(3) = plot(x,gaussFit.a3*exp(-((x-gaussFit.b3)/gaussFit.c3).^2), 'Color', 'g');

% dist_values =   [dists(1).YData];
% dist_values =   [dists(1).YData; dists(2).YData];
dist_values =   [dists(1).YData; dists(2).YData; dists(3).YData];

threshold_def =     'std';
% threshold_def =     'dist';


switch threshold_def
    case 'dist'
        % switched threshold generation to AOC
        areas =     trapz(dist_values');
        
        [~, dist_order] = sort(areas, 'descend');
        
        primary =   dists(dist_order(1));
        second =    dists(dist_order(3));
        
        high_num =  3;
        low_num =   1;
        
        % find the x-val of the point where the higher dist (well-isolated) surpasses the lower dist (poorly isolated)
        thresh =    x( find( dist_values(high_num, :)>dist_values(low_num, :), 1, 'last') )
        
    case 'std'
        thresh =    dist_coeffs(2) + dist_coeffs(3)
end

thresh
line([thresh thresh], [0 10], 'Color', 'k', 'LineWidth', 2);

disp('*** UNIT THRESHOLD FOUND ***');

figure;

for chan_idx = 1 : max_channels
    num_units = length( find( ~isnan(meanDotProd(chan_idx, :) ) ) );
    
    for unit_idx = 1 : num_units
        
        % if not already rejected
        if scrubbed_isolated_mask(chan_idx, unit_idx) == 0
            
            unit_dotprod = meanDotProd(chan_idx, unit_idx);
            
            % isolated unit
            if unit_dotprod >= thresh
                scrubbed_isolated_mask(chan_idx, unit_idx) = 2;
                plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'b');
            else
                scrubbed_isolated_mask(chan_idx, unit_idx) = 1;
                plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'm');
            end
           
            
        end% IF
        
    end% FOR, unit_idx
    
end% FOR, chan_idx

disp('*** UNIS LABELED ***');


%% pull out snippets from different unit classes

figure; hold on;

for chan_idx = 1 : max_channels
    
    num_units = length( find( ~isnan(meanDotProd(chan_idx, :) ) ) );% this is a stupid way to do this here
    
    for unit_idx = 1 : num_units
        % isolated
        if scrubbed_isolated_mask(chan_idx, unit_idx) == 2
            plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'b');
        end
        % rejected by inconsistent waveforms
        if scrubbed_isolated_mask(chan_idx, unit_idx) == 1
            plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'm');
        end
        % rejected by average waveform shape
        if scrubbed_isolated_mask(chan_idx, unit_idx) == -2
            plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'r');
        end
        % rejected by ISI
        if scrubbed_isolated_mask(chan_idx, unit_idx) == -3
            plot(unit_waveforms{chan_idx, unit_idx}(:, 1), 'Color', 'k');
        end
        
    end% FOR, unit_idx
    
end% FOR, chan_idx


%%

sort_save_dir = save_data_path;

switch subject
    case 'subj1'
        save([sort_save_dir '\' date_to_analyze '\' date_to_analyze '_sort_isolation'], 'scrubbed_isolated_mask', 'scrubbed_unit_waveforms', '-v7.3');
end

scrubbed_active_mask = reshape(scrubbed_isolated_mask', max_channels*5, 1);

sorted_task_data.data.SortedActiveChannelMask = scrubbed_active_mask;

disp('*** MASK SAVED ***');


%% process sorted data for analysis

chain_snips = cell(max_channels*5, 1);
chain_times = cell(max_channels*5, 1);
counter = 1;


for chan_idx = 1 : max_channels
    
    for unit_idx = 1 : 5
        
        if ~isempty(new_snippets(chan_idx, unit_idx) )
            chain_snips(counter) = new_snippets(chan_idx, unit_idx);
            chain_times(counter) = new_timing(chan_idx, unit_idx);
        end
        counter = counter + 1;
        
    end
    
end


sorted_task_data.snippets = chain_snips;
sorted_task_data.snip_times = chain_times;

is_craniux = 1;

time_window = [2 2 0.001];


if is_craniux == 1
        task_names = {'task_list'};
else
    task_names = {'Motor_All'};
end

NAPAS_spikesort_preprocData(subject, date_to_analyze, directory, sorted_task_data, time_window, task_names);

disp('*** DONE ***');



