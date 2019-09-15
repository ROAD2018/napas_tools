% 2019-09-15 Dylan Royston
% From 2018-07-06 Dylan Royston
%
% Adapted from Process_open_loop_spikes.m to work with whole-session data spike-sorted with SCRIPT_sort_spikes_PCA
% Said function uses FUNC_spikesort_returndata()
%
%
% INPUTS:
%   Date: date of data recording (year-month-day)
%   Task: specific movement to be processed
%   Directory: Analysis folder, parent directory containing folder where variables are to be saved
%   Time_Window: 3x1 vector of [seconds before move; seconds after move; bin width in seconds]
%
% OUTPUTS:
%   processed_data.binned_spikes: binned spikes organized around stimulus time
%   processed_data.raw_spikes: unbinned, raw spike times organized around stimulus time
%   processed_data.time: time vector of bins
%   processed_data.units: list of active units for indexing spikes
%   processed_data.smoothed_FR: instantaneous firing rate from binned spike counts, smoothed with spike_filter from online decoder
%
%
%%

function data_to_return = FUNC_Preprocess_sorted_spikes(subject, date_to_analyze, analysis_directory, sorted_task_data, time_window, task_names)

%% 1.Load data and set flags
% clear;clc;
display('*** INITIALIZING VARIABLES ***')

clearvars processed_data


% Extracts raw data files from Analysis folder structure (required to exist)

time_before =           time_window(1);
time_after =            time_window(2);
bin_size =              time_window(3);

%% blocked out because data is passed directly, can adapt to load from folder if PCA isn't performed online

data =  sorted_task_data.data;
idata = sorted_task_data.idata;

%% 2.Get indexing information
% pure logistics, dealing with Blackrock Neuroport data structures and timing differences between pedestals
display('*** RETRIEVING INDEX INFORMATION ***')

% get electrode references
channel_index =     double(idata.QL.Data.SPIKE_SNIPPET.ss.channel);
unit_index =        double(idata.QL.Data.SPIKE_SNIPPET.ss.unit);
pedestal_index =    double(idata.QL.Data.SPIKE_SNIPPET.ss.source_index);

switch subject
    case 'subj1'
        max_units =         1280;% number of possible units on arrays
        max_chans =         256;
end


%% pull out snippet data for storage

display('*** EXTRACTING SNIPPETS ***');

% 2018-07-06 Royston: pulling out snippets for plotting later
unit_snippet_holder = sorted_task_data.snippets;

%% THIS RUNS ON WHOLE-SESSION DATA NOW

switch subject
    case 'subj1'
        
        pre_time_S =        time_before;% time before stimulus code
        post_time_S =       10 + time_after;% time after stimulus code
        
        % event codes per Neuroport
        snippet_events =        idata.QL.Data.RAW_DIGITAL_EVENT;
        
        if ~isempty(find( unique(snippet_events.data(1, :)) == 2))
            event_code =            2;% video-onset Craniux code
        else
            event_code =            max(snippet_events.data(1, :) );
        end
        
        
        % 2018-12-05 Royston: patch to deal with sensory tasks having different (offset by 240?) event codes
        if event_code == 2
            temp =              snippet_events.data(1, :);
            temp(temp>240) =    temp(temp>240) - 240;
            event_idx_both =    find(temp == event_code)';
        else
            event_idx_both =        find(snippet_events.data(1, :) == event_code)';
        end
        
 
        
        event_source_idx =      snippet_events.source_index(event_idx_both);
        event_idx_boxes{1} =    event_idx_both(event_source_idx == 0);
        event_idx_boxes{2} =    event_idx_both(event_source_idx == 1);
        
        event_timeS_boxes{1} = snippet_events.source_timestamp(event_idx_boxes{1});
        event_timeS_boxes{2} = snippet_events.source_timestamp(event_idx_boxes{2});
        
        
        
        
        % get spike times
        spike_struct =          idata.QL.Data.SPIKE_SNIPPET.ss;
        raw_spike_times =       spike_struct.source_timestamp;
        
        spike_timeS_boxes{1} =  raw_spike_times(spike_struct.source_index == 0);
        spike_timeS_boxes{2} =  raw_spike_times(spike_struct.source_index == 1);
        
        % 2018-12-12 Royston: added event time resort to deal with (I think) auto-file loading across sessions?
        event_timeS_boxes{1} = sort(event_timeS_boxes{1}, 'ascend');
        event_timeS_boxes{2} = sort(event_timeS_boxes{2}, 'ascend');

        % removes erroneous timing due to Neuroport clock resets
        [~, event_base(1)] =    min(event_timeS_boxes{1});
        [~, event_base(2)] =    min(event_timeS_boxes{2});
        [~, time_base(1)] =     min(spike_timeS_boxes{1});
        [~, time_base(2)] =     min(spike_timeS_boxes{2});
        
        event_timeS_boxes{1} = event_timeS_boxes{1}(event_base(1):end);
        event_timeS_boxes{2} = event_timeS_boxes{2}(event_base(2):end);
        spike_timeS_boxes{1} = spike_timeS_boxes{1}(time_base(1):end);
        spike_timeS_boxes{2} = spike_timeS_boxes{2}(time_base(2):end);
        
        %% 3.Organize spikes around stimulus by unit
        display('*** ORGANIZING SPIKES AROUND STIMULUS ***')
        
        spikes_by_unit = cell(1, 2);
        
        for ibox = 1:2
            
            n_events(ibox) =                size(event_timeS_boxes{ibox},2);
            
            unit_spikes = cell(1, max_units);
            spikes_by_unit{ibox} =          cell(n_events(ibox), max_units);% holds all spike times by box, unit and event
            channel_idx_boxes{ibox} =       channel_index(  pedestal_index == (ibox-1)  );
            unit_idx_boxes{ibox} =          unit_index(  pedestal_index == (ibox-1)  );% creates indexing information for box-specific spikes
            pedestal_idx_boxes{ibox} =      pedestal_index(  pedestal_index == (ibox-1)  ) ;
            
            
            for ievent = 1:n_events(ibox)
                
                trial_spike_idx = find( spike_timeS_boxes{ibox} > (event_timeS_boxes{ibox}(ievent)-pre_time_S)...
                    & spike_timeS_boxes{ibox} <= (event_timeS_boxes{ibox}(ievent)+post_time_S) );% finds spikes in event window
                
                trial_channel_idx =         channel_idx_boxes{ibox}(trial_spike_idx);
                trial_unit_idx =            unit_idx_boxes{ibox}(trial_spike_idx);
                
                % converts ped/chan/unit referencing to 1280-unit index (equation is true, the logic below may be suspect)
                % block 1 compensates for pedestal-channel offset
                % block 2 converts into units/channel
                % block 3 compensates for unit indexing changing from 0 to 1
                trial_real_unit_index =    ( (ibox-1)*128*5) + (trial_channel_idx-1)*5 + (trial_unit_idx + 1);
                
                
                trial_spike_timeS =         spike_timeS_boxes{ibox}(trial_spike_idx);% finds corresponding spike times
                trial_spike_timeS =         trial_spike_timeS - (event_timeS_boxes{ibox}(ievent));% normalizes spike times to event
                
                
                unique_units = unique(trial_real_unit_index);
                
                for iunit = 1:length(unique_units)
                    
                    unit_reference =        unique_units(iunit);
                    unit_reference_idx =    find(trial_real_unit_index == unit_reference);
                    unit_spikes =           trial_spike_timeS(unit_reference_idx);% assigns spikes to each corresponding unit
                    
                    spikes_by_unit{ibox}{ievent,unit_reference} = unit_spikes; % cell array containing spikes per unit per event
                end
            end
            
        end
        
        
        spikes_by_unit_all = cell(max(n_events),max_units);
        
        for ibox=1:2
            for ievent=1:n_events(ibox)
                for iunit = 1:max_units
                    % compiles all spikes into box-agnostic array
                    spikes_by_unit_all{ievent,iunit} = [spikes_by_unit_all{ievent, iunit} spikes_by_unit{ibox}{ievent, iunit}];
                end
            end
        end
                
        
        if event_code == 2% if craniux
            
            
            
        else% if PTB
            
            
            if isfield(idata.QL.Data, 'STIM_PRESENT')

                clearvars stim_names
                stimulus_filenames =    cast(idata.QL.Data.STIM_PRESENT.stim_filename', 'char');

                for stim_idx = 1 : size(stimulus_filenames, 1)
                    curr_name =                 stimulus_filenames(stim_idx, :);
                    characters =                unique(curr_name);
                    blank_spaces =              strfind(curr_name, characters(1));
                    curr_name(blank_spaces) =   [];
                    stim_names{stim_idx} =      curr_name;
                end

                % restrict events to video onset
                video_idx =     FUNC_find_string_in_cell(stim_names, '.wmv');
                stim_names =    stim_names(video_idx);

                spikes_by_unit_all = spikes_by_unit_all(video_idx, :);

                n_events =              length(stim_names);
            
            end% IF, isfield
            
        end% IF, event_code == 2
        
        
        % 2019-03-08 Royston: debug figure, patching to verify PsychoPy events work correctly
%         figure; 
%         line([event_timeS_boxes{1}' event_timeS_boxes{1}'], [0 1]);
        time_intervals = round( diff(event_timeS_boxes{1}) );
        
        betweens =      find(time_intervals > 30);
        between_start = event_timeS_boxes{1}(betweens) + 10;
        
        moves =   		find(time_intervals == 14);
        move_start = 	event_timeS_boxes{1}(moves);
        
        
        
        all_move_onsets = sort(horzcat(moves, betweens, length(event_timeS_boxes{1}) ), 'ascend');
        n_events =      length(all_move_onsets);
        
        
        %% 4.Bin spike data
        display('*** BINNING SPIKES ***')
        
        
        bin_width =     bin_size;% default = 0.02 (20ms, same as pre-binned data)
        bin_time =      -pre_time_S:bin_width:post_time_S;%change to percent of time
%         pc_length =     10/bin_width;% length of video in bin-indexes (default 10s, 500 indexes)
        
        active_units =  find(data.SortedActiveChannelMask > 0);
        
        spikes_binned =     zeros(length(bin_time), max_units, max(n_events));
        spikes_binned =     zeros(length(bin_time), length(active_units), max(n_events));
        zero_marker =       find(bin_time == 0);% index of event
        
        for bindex = 1:(length(bin_time)-1)
            
            disp([num2str(bindex) '\' num2str(length(bin_time)) ] );
            
            bin_top =       bin_time(bindex+1);
            bin_bottom =    bin_time(bindex);
            
%             for iunit = 1:max_units
            for iunit = 1: length(active_units)
                
                current_unit = active_units(iunit);%active_units(iunit, 1);% bins only active units
                
                for ievent = 1:max(n_events)
                    
                    % counts number of spikes per bin per unit for each event
                    n_spikes = sum( spikes_by_unit_all{ievent,current_unit} > bin_bottom &...
                        spikes_by_unit_all{ievent,current_unit} <= bin_top );
                    
                    spikes_binned(bindex,iunit, ievent) = n_spikes;
                end
            end
        end
        
        
        
        %% 2018-07-06 Royston: organize data into tasks
        
        % 2018-12-05 Royston: patch to deal with n_events being more than one value
        if length(n_events)>1
            n_events = n_events(1);
        end
        
        if event_code == 2
            
            task_labels = task_names;
            
            unique_tasks = unique(task_labels);
            
            task_indexing =     reshape( [1:n_events], [8 length(task_names)] );
                        
            for cond_idx = 1 : length(task_names)
                                
                structured_task_data{cond_idx} = spikes_binned(:, :, task_indexing(:, cond_idx)' );
                structured_raw_spikes{cond_idx} = spikes_by_unit_all(task_indexing(:, cond_idx)', :);
                
            end% FOR, cond_idx
            
        else% PTB
            
            clearvars structured_task_data structured_raw_spikes task_labels
            counter = 0;
            
            for task_idx = 1 : length(task_names)
                
                current_task =  task_names{task_idx};
                
                name_parts =    strsplit(current_task, '_');
                
                % if task is not "All", it will have 8 repetitions
                if ~strcmp(name_parts{2}, 'All')
                    
                    counter = counter+8;
                    
                    task_max_ind =      counter;
                    task_min_ind =      counter-7;
                    
                 if ~exist('stim_names')
                        test_names =    name_parts(2);
                    else
                        test_names =    stim_names(task_min_ind:task_max_ind);
                 end
                    
                    if length(uniqueStrCell(test_names))==1
                        
                        structured_task_data{task_idx} =    spikes_binned(:, :, task_min_ind:task_max_ind);
                        structured_raw_spikes{task_idx} =   spikes_by_unit_all(task_min_ind:task_max_ind, :);
                        task_labels{task_idx} =             current_task;
                    end
                    
                % if all, has 2-4 reps of 9 movements
                else
                    
                    counter = counter + 18;
                    
                    task_max_ind =      counter;
                    task_min_ind =      counter-17;
                    
                    test_names =    stim_names(task_min_ind:task_max_ind);
                    
                    name_counts =   FUNC_count_string_vals(test_names);
                    name_counts =   [name_counts{:}];
                    
                    if length(unique(name_counts))==1
                        structured_task_data{task_idx} =    spikes_binned(:, :, task_min_ind:task_max_ind);
                        structured_raw_spikes{task_idx} =   spikes_by_unit_all(task_min_ind:task_max_ind, :);
                        task_labels{task_idx} =             test_names;
                    end
                    
                end
                
                
            end% FOR, task_idx
            
            
            
        end% IF, event_code
        
        
        organized_labels = task_labels;
        organized_raw_data = structured_raw_spikes;

end% SWITCH, subject



% %% FIGURE: Plot average raster of all units
%         display('*** PLOTTING AVERAGE RASTER ***')
%
%         all_trials =    mean(spikes_binned(:, :, 4), 3);
%         bin_time =      bin_time(1:end);
%
%         Plot_Rasters_CRS02(all_trials,bin_time,active_units,0);
%         caxis([0 3])
%%

% smoothing = 'Gaussian';

organized_labels = task_names;

for task_idx = 1 : length(organized_labels)
    
    task_name =     organized_labels{task_idx};
    
    if iscell(task_name)
        task_title = task_names{task_idx};
        sublabels = task_name;
    else
        task_title = task_name;
        sublabels = [];
    end
    


    bin_time =      -pre_time_S:bin_size:post_time_S;%change to percent of time
    bin_time(end) = [];
    
    if ~exist('active_units')
        active_units = find(data.SortedActiveChannelMask>0);
    end

        
    data_to_return(task_idx).snippets =       sorted_task_data.snippets;
        
    data_to_return(task_idx).binned_spikes =  structured_task_data{task_idx};
    if exist('organized_raw_data')
        data_to_return(task_idx).raw_spikes =     organized_raw_data{task_idx};
    end
    data_to_return(task_idx).time =           bin_time;
    data_to_return(task_idx).units =          active_units;
    data_to_return(task_idx).unit_label =     data.SortedActiveChannelMask;
    data_to_return(task_idx).subtasks =       sublabels;
    
    processed_data = data_to_return(task_idx);
    
    saved_data_string = [analysis_directory '\task_files\' date_to_analyze '_' task_title '_PREPROCESSED.mat'];

    
    disp([num2str(task_idx) '/' num2str(length(organized_labels))]);
    
    save(saved_data_string, 'processed_data', '-v7.3');
    
    clearvars processed_data
    
end% FOR, task_idx

sort_info_string = [analysis_directory '\' date_to_analyze '_SNIP_SORTS_PREPROCESSED.mat'];
save(sort_info_string, 'unit_snippet_holder', '-v7.3');


display('*** COMPLETED PREPROCESSING ***');

