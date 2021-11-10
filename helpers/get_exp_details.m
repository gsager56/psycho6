function [exp_info, raw_movie] = get_exp_details(dataPath, param)
% this function takes in a cell array of paths and outputs all the
% experimental details of the recording in exp

    % let's check some stuff in the param structure to try and help the
    % user know if they might be making a mistake
    if isfield( param, 'epochs_for_selectivity' ) & isfield( param, 'corr_idxs' )
        
        error(['ERROR: Please either leave corr_idxs empty or epochs_for_selectivity empty. ',...
            'This helps ensure you do not make any unintended mistakes'])
    end
    
    
    % load in the param file information - this may not be necessary
    %tmp = load(fullfile(dataPath, 'stimulusData', 'chosenparams.mat'));
    %exp.params = tmp.params;
    
    % load in the run deatils - this may not be necessary
    %tmp = load(fullfile(dataPath, 'stimulusData', 'runDetails.mat'));
    %exp.runDetails = tmp;
    
    % load in image description
    tmp = LoadImageDescription(dataPath);
    exp_info.state = tmp;%%%
    linescan = ~exp_info.state.acq.scanAngleMultiplierSlow;
    
    % load in photodiode info
    [exp_info.photodiode.data, exp_info.photodiode.highResLinesPerFrame] = ReadInPhotodiode(exp_info.state, dataPath); %%%
    
    % get epoch list
    [epochBegin, epochEnd, endFlash, flashBeginInd] = GetStimulusBounds( exp_info.photodiode.data, ...
        exp_info.photodiode.highResLinesPerFrame, exp_info.state.acq.frameRate, linescan );
    
    if length(endFlash)>1
        error(['ERROR: more than one photodiode flash was 20 projector frames long.\n',...
            ' Check what psycho5 says could be the error in ReadImagingData.m\n',...
            ' when length(endFlash)>1'])
    end
    
    % store when epoch begins and ends
    exp_info.epochBegin = floor(epochBegin) + 1;
    exp_info.epochEnd = floor(epochEnd) + 1;
    
    % load in the stimulus data
    tmp = load(fullfile(dataPath, 'stimulusData', 'stimdata.mat'));
    stimData_time = tmp.stimData(:,1);
    stimData_epoch = tmp.stimData(:,3);
    
    stimData_time( stimData_epoch==0 ) = [];
    stimData_epoch( stimData_epoch==0 ) = [];
    
    % get rid of stuff past the last photodiode flash    
    stimData_epoch( stimData_time > (epochEnd-epochBegin)/exp_info.state.acq.frameRate ) = [];
    stimData_time( stimData_time > (epochEnd-epochBegin)/exp_info.state.acq.frameRate ) = [];
    
    % interpolate new list of epoch
    exp_info.time = linspace(stimData_time(1),stimData_time(end),floor(epochEnd)-floor(epochBegin)+1);
    exp_info.epochVal = round(interp1(stimData_time,stimData_epoch,exp_info.time,'nearest'))';
    
    % find out which channels where being acquired
    acqCell = struct2cell(exp_info.state.acq);
    acqFldNames = fieldnames(exp_info.state.acq);
    chanAcqState = [acqCell{~cellfun('isempty', regexp(acqFldNames, 'saving.*(\d+)'))}];
    acquiredChannels = find(chanAcqState);

    [wasAcq,desiredChanIdxs] = ismember(param.movieChan,acquiredChannels);
    
    if isempty(acquiredChannels)
        warning('The acquired channels was not stored, so I will assume channels 1 and 3 were acquired')
        acquiredChannels = [1, 3];
        desiredChanIdxs = 1;
    elseif ~all(wasAcq)
        error('At least one of the requested channels from variable channelsDesired were not actually acquired');
    end
    numChannelsAcquired = length(acquiredChannels);

    originalDataPath = DirRec(dataPath,'.tif');
    [~,thisFile] = fileparts(originalDataPath{1});
    moviePath = fullfile(dataPath,[thisFile '.tif']);

    for movieChan = 1:length(param.movieChan)
        raw_movie(:, :, :, movieChan) = LoadTiffStack(moviePath,[desiredChanIdxs(movieChan) numChannelsAcquired]);
    end
    
    % only take the part of the movie when stimulus was presented
    raw_movie = raw_movie(:,:, exp_info.epochBegin : exp_info.epochEnd, :);
    
    
end