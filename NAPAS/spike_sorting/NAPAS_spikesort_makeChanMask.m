% 2018-04-10 Dylan Royston
%
% Copied from hst2/Analysis Code/Stability Code/SortedSpikeCounts_SortedActiveChanMask
%
% Function to take spike-sorted data and reassign information in Data fields according to new mask
%
%
%
% === UPDATES ===
%%

function [data] = FUNC_SortedSpikeCounts_SortedActiveChanMask(data,exp_time, subject)

switch subject
    case 'CRS02b'
        max_units =     1280;
        binwidth =      0.02;
    case 'BMI01'
        max_units =     960;
        binwidth =      0.03;
%         binwidth =      0.001;
end

% extract and initialize reassigned spike timing information
SpikeTimes =    exp_time;
SpikeTimes =    SpikeTimes';
SpikeTimes =    reshape(SpikeTimes,max_units,1);
%SpikeTimes(cellfun('length',SpikeTimes)<4)={[]};
%arraySwitch = 640-sum(double(cellfun(@isempty,SpikeTimes(1:640))));
%tempSpikeTimes = SpikeTimes(~cellfun('isempty',SpikeTimes));



min1 =          min( cellfun(@min, SpikeTimes(~cellfun('isempty', SpikeTimes(1:(max_units/2)) ) ) ) );
min2 =          min( cellfun(@min, SpikeTimes( [false((max_units/2),1); ~cellfun('isempty', SpikeTimes((max_units/2 + 1):max_units) ) ] ) ) );

max1 =          max( cellfun(@max, SpikeTimes(~cellfun('isempty', SpikeTimes(1:(max_units/2)) ) ) ) );
max2 =          max(cellfun( @max,SpikeTimes( [ false((max_units/2),1); ~cellfun('isempty', SpikeTimes((max_units/2 + 1):max_units) ) ] ) ) );


time =          max( [max1 - min1 max2 - min2] );

% time =          max( [max(cellfun(@max, SpikeTimes(~cellfun('isempty',SpikeTimes(1:(max_units/2)) ) ) ) ) - min1 ...
%                     max(cellfun(@max,SpikeTimes( [false((max_units/2),1); ~cellfun('isempty', SpikeTimes((max_units/2 + 1):max_units) ) ] ) ) ) - min2] );

numofunits =    length(SpikeTimes);
disp(['Max: ' num2str(numofunits)]);

% bins =          size(data.SpikeCount,1);
% SpikeCount =    zeros(bins,numofunits);
% binwidth =      time/bins;
% 

% binwidth =      0.02;
% binwidth =      0.001;
bins =          time/binwidth;
% SpikeCount =    NaN(round(bins), numofunits);
SpikeCount =    zeros(round(bins), numofunits);
% 
% num_orig_bins = size(data.SpikeCount, 1);
% 
% fake_time =     (1:num_orig_bins)*0.03;
% 
% % 2019-02-26 Royston: testing better way to organize spikes around onset
% prebin_events =     data.TaskStateMasks.state_num;
% prebin_size =       0.03;
% 
% 
% counter = 1;
% counter2 = 1;
% move_set = [];
% for bin_idx = 2 : length(prebin_events)
%     
%     if (prebin_events(bin_idx) == 5) && (prebin_events(bin_idx-1) == 4)
%         move_set(1, counter) = bin_idx;
%         
%         counter = counter + 1;
%     end
%     
%     if (prebin_events(bin_idx) == 4) && (prebin_events(bin_idx-1) == 5)
%         move_set(2, counter2) = bin_idx;
%         
%         counter2 = counter2 + 1;
%         
%     end% IF
%     
% end% FOR, bin_idx
% 
% move_times = (move_set./prebin_size);
%  
% new_time = 0.001:0.001:time;
% 
% [val, idx] = FUNC_find_closest_value(new_time, 27);





if binwidth > .022 || binwidth < .018
    str = sprintf('Binwidth is unexpected: %.4d',binwidth);
    disp(str);
end

% normalize spike times to neuroport zeros
for i = 1:(max_units/2)
    SpikeTimes(i) = {SpikeTimes{i}-min1};
end

for i = (max_units/2 + 1):max_units
    SpikeTimes(i) = {SpikeTimes{i}-min2};
end

tic;
% assign spike times to new unit IDs
for i = 1:numofunits
    
    rounded_spikes =    round(SpikeTimes{i}, 3);
    
    if ~isempty(rounded_spikes)
        
        for j= 1:bins
            
            rounded_bins =      [round((j-1)*binwidth, 3) round(j*binwidth, 3)];
            
            SpikeCount(j,i) = sum( double( rounded_spikes >= rounded_bins(1) & rounded_spikes < rounded_bins(2) ) );
            %             test_vector(j) = sum( double(SpikeTimes{i} >= (j-1)*binwidth & SpikeTimes{i} < j*binwidth) );
        end
    end
    
    disp(i);
end
toc

 


 blankcols = find(all(SpikeCount'==0));
 neighbors = diff(blankcols);
 
 N =    60;% number of consecutive values to find
 x =    neighbors == 1;% indices of target values
 f =    find([false,x]~=[x,false]);% indices of blocks (>=1) of target value
 
 block_size = f(2:2:end)-f(1:2:end-1);
 
 [maxes, locs] = sort(block_size, 'descend');
 
 target_blocks = find(maxes>N);
 
 block_start = f(2*locs([target_blocks])-1);
  to_blank = [];
 for block_idx = 1 : length(target_blocks)
     offset(block_idx) = block_start(block_idx) + maxes(block_idx);
     this_blank = block_start(block_idx) : offset(block_idx);
     to_blank = horzcat(to_blank, this_blank);
 end
 
% to_blank = [ block_start(2) : ( block_start(2) + maxes(2) ), block_start(1) : block_start(1) + maxes(1) ];
 

figure; hold on;
plot(neighbors);
for block_idx = 1 : length(target_blocks)
    line([block_start(block_idx) block_start(block_idx)], [0 18], 'Color', 'r');
    line([offset(block_idx) offset(block_idx)], [0 18], 'Color', 'r');
end


 
%  to_blank =  blankcols(find(neighbors == 1) + 1);

%  figure; line([blankcols' blankcols'], [0 1], 'Color', 'k');
% figure; plot(blankcols); xlim([0 length(SpikeCount)]);

%  to_blank =  blankcols(find(neighbors < 120) + 1);
  to_blank =  blankcols(find(neighbors < 140) + 1);
 
 clean_spikes = SpikeCount;
 clean_spikes(to_blank, :) = [];
 
 new_active_mask = logical(nansum(clean_spikes, 1));
 
%  figure; 
%  imagesc(clean_spikes(:, new_active_mask)');
%  colormap(flipud(gray));
%  caxis([0 3]); 
%  
%  
%  % event marks
%  snippet_events =        data.TaskStateMasks.state_num;
%  
%  onset_marker = [];
%  offset_marker = [];
%  counter1 = 0;
%  counter2 = 0;
%  for idx = 1 : length(snippet_events)-1
%      
%      if (snippet_events(idx) == 4) && (snippet_events(idx+1) == 5)
%          counter1 = counter1 + 1;
%          onset_marker(1, counter1) = idx;
%          onset_marker(2, counter1) = data.trial_num(idx);
%      elseif (snippet_events(idx) == 5) && (snippet_events(idx+1) == 4)
%          counter2 = counter2 + 1;
%          offset_marker(1, counter2) = idx;
%          offset_marker(2, counter2) = data.trial_num(idx);
%      end
%      
%  end
%  
%  Plot_VerticalMarkers(onset_marker(1, :), 'Color', 'g');
%  Plot_VerticalMarkers(offset_marker(1, :), 'Color', 'r');

% create new Active mask
% data.SpikeCount = SpikeCount;
data.new_SpikeCount = clean_spikes;
% data.SortedActiveChannelMask = logical(sum(data.SpikeCount, 1));
data.SortedActiveChannelMask = new_active_mask;


end% FUNCTION

