function [exp_info, raw_movie, param] = get_exp_details(dataPath, param)
% this function takes in a cell array of paths and outputs all the
% experimental details of the recording in exp
    
    
    % load in the param file information - this may not be necessary
    tmp = load(fullfile(dataPath, 'stimulusData', 'chosenparams.mat'));
    exp.params = tmp.params;
    
    % load param_file
    exp.param_file = readmatrix(fullfile(dataPath, 'stimulusData', param.stim), 'OutputType', 'string');
    
    if ischar(param.probe_epochs)
        assert( strcmp(param.probe_epochs,'do it for me'), 'Could not understand what you mean in param.probe_epochs' )
        row = -1;
        for i = 1 : size(exp.param_file,1)
            if strcmp(exp.param_file{i,1}, 'Stimulus.epochName')
                row = i;
            end
        end
        assert(row > 0, 'Could not find the column Stimulus.epochName in the param file. Try hard coding probe epocs to avoid this error')

        for i = 2 : size(exp.param_file,2)
            epoch_name = exp.param_file{row,i}; % name of epoch in param file
            found_epoch = false;
            param_epoch_names = {exp.params.epochName};
            for j = 1 : size(exp.params,2) % loop through all epoch names in the probe and param file
                if strcmp(param_epoch_names{j},epoch_name)
                    found_epoch= true;
                end
            end
            assert( found_epoch, ['Could not find the epoch ', epoch_name, ' in the list of presented epochs. \n', ...
                'Either figure out the issue or avoid this error by hard coding the probe epochs'] )
        end

        num_param_epochs = size(exp.param_file,2)-1;

        num_probe_epochs = size(exp.params,2) - num_param_epochs;
        param.probe_epochs = 1 : num_probe_epochs;
    end
    
    
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
    
    if length(param.movieChan) > 1
        % we'll trust the user knows whatthey're doing
        warning('The acquired channels was not stored, so I will assume channels 1, 2, and 3 were acquired')
        acquiredChannels = [1, 2, 3];
        desiredChanIdxs = param.movieChan;
    elseif isempty(acquiredChannels)
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