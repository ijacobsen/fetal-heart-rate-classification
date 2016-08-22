% [med, medDev, stv, sti, msti, poinstd1, poinstd2, poinccm higuch, psd1, psd2, psd3, psd4, psd5, ltv, delta, sampEn, fuzzyEn, mFHR, sdFHR, LTI, STV, II]
% [ 1     2      3    4     5      6          7         8      9     10    11    12    13    14    15    16     17      18      19     20    21,  22,  23]  
function fhrFeatures = getFHRfeatures(fhr, featureList, lastMinutes, timeSrsLengthMin, overlap, fs, numMinsDelete)

% get last n minutes of data
nMinutesBack = lastMinutes*fs*60;                     % 10 minutes = 4Hz*60sec*10min = 2400
deleteLastMins = numMinsDelete*4*60;
for iter=1:length(fhr(:, 1))
    entireData(iter, :) = fhr(iter, end-nMinutesBack+1:end-deleteLastMins);
end

numTimeSeries = lastMinutes-overlap-numMinsDelete;        % *** WARNING: not error checking! must be integer (how many time series)
trainingSize = length(entireData(:, 1));                   % number of training examples
timeSrsLengthSamples = timeSrsLengthMin*60*fs;
overlapSamples = overlap*60*fs;


% initial (no overlap)
data = entireData(:, 1:timeSrsLengthSamples);
tempFeat = zeros(trainingSize, 1);                  % this column will be deleted at the end (place holder)
for iter=1:length(featureList)
    fptr = str2func(featureList{iter, 1});
    args = num2cell(featureList{iter, 2});
    fet = fptr(data, args{:});
    tempFeat = [tempFeat, fet];
end
tempFeat = tempFeat(:, 2:end);                    % remove 0's column
fhrFeatures(:, :, 1) = tempFeat;

% remaining (with overlap)
for loopNum = 2:numTimeSeries
    data = entireData(:, ((loopNum-1)*overlapSamples + 1):(loopNum+1)*overlapSamples);
    tempFeat = zeros(trainingSize, 1);                  % this column will be deleted at the end (place holder)
    for iter=1:length(featureList)
        fptr = str2func(featureList{iter, 1});
        args = num2cell(featureList{iter, 2});
        fet = fptr(data, args{:});
        tempFeat = [tempFeat, fet];
    end
    tempFeat = tempFeat(:, 2:end);                    % remove 0's column
    fhrFeatures(:, :, loopNum) = tempFeat;
end
