% 2018-07-06 Dylan Royston
%
% Copied from jeff's spikeSortToIDataSS
%
% Takes in (spikes, times) from FUNC_SpikeSort_returnData and replaces relevant fields in idata for later processing
%
%
%
%
%
%%

function ss = FUNC_insert_sorted_idata(subject, spikes,times)

% get spike sort output back into idata format

switch subject
    case 'BMI01'
        max_channels = 192;
    case 'CRS02b'
        max_channels = 256;
end

% squeeze spikes and times arrays
spikes =                squeeze(spikes);
times =                 squeeze(times);

spksPerUnit =           cellfun(@numel,times);
nSpikes =               sum(sum(spksPerUnit)); % total number of spikes

% initialize ss struct
ss.source_index =       zeros(1,nSpikes,'int32');
ss.channel =            zeros(1,nSpikes,'int16');
ss.unit =               zeros(1,nSpikes,'uint8');
ss.source_timestamp =   zeros(1,nSpikes);
ss.snippet =            zeros(48,nSpikes,'single');

spkIdx =                0;

pedestal_size =         max_channels/2;

% loop and populate ss structure
for iSrc = 0:1
    for iChan = 1:pedestal_size
        chanIdx = iChan + pedestal_size*iSrc;
        for iUnit = 1:5
            if spksPerUnit(chanIdx,iUnit)
                
               nSpksThisUnit =                      spksPerUnit(chanIdx,iUnit);
               idx1 =                               spkIdx+1;
               spkIdx =                             spkIdx + nSpksThisUnit;
               
               ss.source_index(idx1:spkIdx) =       int32(iSrc);
               ss.channel(idx1:spkIdx) =            int16(chanIdx - pedestal_size*iSrc);
               ss.unit(idx1:spkIdx) =               uint8(iUnit-1);
               ss.source_timestamp(idx1:spkIdx) =   times{chanIdx,iUnit};
               ss.snippet(:,idx1:spkIdx) =          spikes{chanIdx,iUnit};
               
            end% IF, spksPerUnit
        end% FOR, iUnit
    end% FOR, iChan
end% FOR, iSrc

% sort by time (this will also sort by source, unfortunately)
[ss.source_timestamp,sortIdx] =     sort(ss.source_timestamp);
ss.source_index =                   ss.source_index(sortIdx);
ss.channel =                        ss.channel(sortIdx);
ss.snippet =                        ss.snippet(:,sortIdx);
ss.unit =                           ss.unit(sortIdx);

end% FUNCTION
