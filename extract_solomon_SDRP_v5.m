function [curr_dat, trial_per_trial, rat_ID, time_wall_act, time_freezing_act, time_rearing_act] =  extract_solomon_SDRP_v5(varargin)

%% Nomenclature

% EVENTS
% 9 = Session Start
% 10 = Session Stops

% BEHAVIOR
% 1  = Freezing Actor
% 2  = NosePoke 
% 3  = High Value Lever Press
% 4  = Low Value Lever Press
% 5  = High Value Reward Poke
% 16 = Low Value Reward Poke
% 6  = Grooming
% 7  = Proximity to Wall
% 8  = Rearing
% 17 = Nose Poke Wrong

% PARTNER BEHAVIOR

% 11 = Freezing Partner
% 12 = Proximity To Wall
% 13 = Rearing Partner
% 14 = Grooming

a = figure; % Latency Lever
b = figure; % Latency Rwd Poke
c = figure; % Time Spent Close to Wall
d = figure; % Freezing

%% This section gets LogFiles of interest
if isequal(nargin,1)
    if isfolder(varargin{1})
        files = dir(fullfile(varargin{1}, '*.csv'));
%         files = files(strncmp('Actor',{files(:).name},1));
    elseif exist(varargin{1}, 'file')
        [datadir,fn,~] = fileparts(varargin{1});
        files = dir(datadir);
        files = files(strcmp(fn,{files(:).name}));
    else
        error('input is not a file nor a directory')
    end
elseif nargin > 1
    extract_varargin % read configurable input flags that could overwrite defaults 
    %error('too many input arguments')
else
    [fn, datadir, outc] = uigetfile('', 'Select a Solomon Datafile');
    if isequal(outc, 0)
        error('no file selected')
    else
        files = dir(datadir);
        files = files(strcmp(fn,{files(:).name}));       
    end
end


%% This section extracts the data from logfile for each rat

    rats_per_cage = 4; disp(['There are ' num2str(rats_per_cage) ' rats per cage. '])
    cage_number = 3;   disp(['There are ' num2str(cage_number) ' cages. '])
    % Create an index of rat ID
    rat_ID = linspace(11,10*rats_per_cage+1,rats_per_cage);
    rat_ID = repmat(rat_ID ,cage_number,1);
    rat_ID(2:end,:) = rat_ID(2:end,:) + repmat([1:(cage_number-1)]',1,rats_per_cage);
    rat_ID = rat_ID'; rat_ID = rat_ID(:);
    % Total rat number
    nRats         = length(rat_ID); 
    % Expected number of trials
    nTrials_expect = 24;
    % Expected number of sessions
    nSess_expect = 4;
    % PreAllocation
    curr_dat        = cell(nRats, 4);
    idx_plot        = NaN(nRats,nSess_expect);
    trial_per_trial = struct('Rat_ID', NaN(1,1),...
                             'Choice', NaN(nTrials_expect,1),...
                             'Trial_Dur', NaN(nTrials_expect,1),...
                             'Latency_Lever', NaN(nTrials_expect,1),...
                             'Rearing', NaN(nTrials_expect,1),...
                             'TimeSpentCloseWall', NaN(nTrials_expect,1),...
                             'TimeSpentFreezing', NaN(nTrials_expect,1));
                     
for iRat = 1 : nRats
    curr_files = files(strncmp(['Rat_' num2str(rat_ID(iRat))],{files(:).name},6));
    
for iFile = 1 : length(curr_files)
    
    % Extract filename
    filename  = curr_files(iFile).name;
    % Extract current session data  
    curr_sess = str2double(filename(strfind(filename,'S')+1));
    % Add current session to the number of baseline sessions
    if contains(filename,'PPP')
    % Extract session # from shock filename and add + 1 to plot after baseline sessions 
    curr_sess = curr_sess + 1;
    else 
    curr_sess = 1;
    end
    % Idx used for ploting data later on
    idx_plot(iRat,curr_sess) = curr_sess;
    % Store rat ID in output array
    trial_per_trial(iRat,curr_sess).Rat_ID = rat_ID(iRat);
    % Extract data
    curr_dat{iRat,curr_sess} = xlsread(filename); 
    % Eliminate data after end session
    end_array = find(curr_dat{iRat,curr_sess}(:,3) == 11); 
    % Eliminate data before start session  
    start_array = find(curr_dat{iRat,curr_sess}(:,3) == 10);  
    % Check that last data point corresponds with end session
    if ~isequal(curr_dat{iRat,curr_sess}(end_array,3), 11)
       error('Last data point does not correspond with end session');
    end  
    % Select data of interest
    curr_dat{iRat,curr_sess} = curr_dat{iRat,curr_sess}(start_array:end_array,1:size(curr_dat{iRat,curr_sess},2));  
    % Bring timing data back to 0
    curr_dat{iRat,curr_sess}(:,1) = curr_dat{iRat,curr_sess}(:,1) - curr_dat{iRat,curr_sess}(1,1);

    if ~isequal((idx_plot(iRat,curr_sess)), false)     
    % Use nosepokes to timestamp the trials   
    if isequal(mod(iRat,2),true) % If the rat is an actor, use its own nosepoke to stamp
    nosepoke_idx = [find(curr_dat{iRat,idx_plot(iRat,(idx_plot(iRat,curr_sess)))}(:,2) == 2);...
                    find(curr_dat{iRat,idx_plot(iRat,(idx_plot(iRat,curr_sess)))}(:,2) == 11)];
    else % Otherwise keep the nosepoke timestamp from previous actor
    nosepoke_idx = [find(curr_dat{iRat-1,idx_plot(iRat-1,(idx_plot(iRat-1,curr_sess)))}(:,2) == 2);...
                    find(curr_dat{iRat-1,idx_plot(iRat-1,(idx_plot(iRat-1,curr_sess)))}(:,2) == 11)];
    end 
                        
    % Extract number of trials performed            
    nTrials = length(nosepoke_idx);
    % Detect unperformed trials 
%     if isequal(mod(iRat,2),false)
    if ~isequal(nTrials, nTrials_expect)
        disp(['Uncomplete Session For Rat ' num2str(iRat) ' in Session ' num2str(curr_sess)...
              '. Rat performed ' num2str(nTrials) ' Trials'])
    end
%     end
    % For each trial, extract the variables of interest
    for iTrial = 1 : nTrials-1
    temp = curr_dat{iRat,idx_plot(iRat,(idx_plot(iRat,curr_sess)))}(nosepoke_idx(iTrial):nosepoke_idx(iTrial+1)-1,:);    
    % In case not last columns behavior were not filled, add NAN columns to allow extraction of no behav
    if size(temp,2) ~= 11 % If all behaviors were detected, excel file will have 11 columns 
    temp(:,size(temp,2):11) = NaN(size(temp,1),length(size(temp,2):11));        
    end
         
    trial_per_trial(iRat,curr_sess).Trial_Dur(iTrial,:) = max(temp(:,1)) - temp(1,1);
    trial_per_trial(iRat,curr_sess).Trial_Dur(end+1:nTrials_expect,1) = NaN;  

        % High  Value Press
        curr_time = temp(temp(:,2) == 3,1) - temp(1,1);
        if ~isempty(curr_time)
        trial_per_trial(iRat,curr_sess).Latency_Lever(iTrial,:) = max(curr_time) - curr_time(1,1);
        trial_per_trial(iRat,curr_sess).Choice(iTrial,:) = 2; % High Value Choice
        end

        % Low  Value Press
        curr_time = temp(temp(:,2) == 4,1) - temp(1,1);
        if ~isempty(curr_time)
        trial_per_trial(iRat,curr_sess).Latency_Lever(iTrial,1) = max(curr_time) - curr_time(1,1);
        trial_per_trial(iRat,curr_sess).Choice(iTrial,:) = 1; % Low Value Choice
        end

        % Rearing
        if ~isempty(temp(temp(:,7) == 9,1)) % Is the coding for actor
        curr_time = temp(temp(:,7) == 9,1) - temp(1,1);
        else                                % Is the coding for partner
        curr_time = temp(temp(:,8) == 14,1) - temp(1,1);    
        end
        if ~isempty(curr_time)
        trial_per_trial(iRat,curr_sess).Rearing(iTrial,1) = max(curr_time) - curr_time(1,1);
        end


        % Time Spent Close to Wall
        if ~isempty(temp(temp(:,2) == 8,1)) % Is the coding for actor
        curr_time = temp(temp(:,2) == 8,1) - temp(1,1);
        else                                % Is the coding for partner
        curr_time = temp(temp(:,4) == 13,1) - temp(1,1);    
        end
        if ~isempty(curr_time)
        trial_per_trial(iRat,curr_sess).TimeSpentCloseWall(iTrial,1) = max(curr_time) - curr_time(1,1);
        end

        % Time Spent Freezing
        if ~isempty(temp(temp(:,6) == 1,1))
        curr_time = temp(temp(:,6) == 1,1) - temp(1,1);
        else
        curr_time = temp(temp(:,5) == 12,1) - temp(1,1);    
        end
        if ~isempty(curr_time)
        trial_per_trial(iRat,curr_sess).TimeSpentFreezing(iTrial,1) = max(curr_time) - curr_time(1,1); 
        end     
    end
    end
    
if isequal(length(curr_files) , false)
disp(['Rat ' num2str(iRat) ' cannot be extracted.'])    
else
disp(['Rat ' num2str(iRat) ' extracted correctly with nTrials = ' num2str(nTrials) ' for nSess = ' num2str(iFile)])
end
end   

end


%% This loop plots the variable
RAT = true;
% after_shock = true if trials after / before shock should be analyzed
after_shock = false;
    if isequal(after_shock,true)
    disp('Data is showed for trials after shock & no shock trials')
    else
    disp('Data is showed for shock & no shock trials')
    end
for iRat = 1 : nRats;
while RAT

for iSess = 1 : size(idx_plot,2);
    % Do not process that session if not performed or data absent
    if isnan(idx_plot(iRat,iSess))
        break
    end       
    % Create index of shock trials
    if isequal(mod(iRat,2),true) % If the rat is an actor, use its own nosepoke to stamp
       idx = trial_per_trial(iRat,iSess).Choice == 2;
    else % Otherwise keep the nosepoke timestamp from previous actor
       idx = trial_per_trial(iRat-1,iSess).Choice == 2;
    end 
    
    % idx_shift will select trials after shock, not shock trials 
    if isequal(after_shock,true)
    idx_shift = zeros(size(idx)); 
    n = 1; % Shift units
    idx_shift(n+1:end) = idx(1:end-n);
    idx = logical(idx_shift);
    disp('Data showed for trials after shock, not shock trials')
    end
     
    % Latency Lever
    figure(a)
    trial_per_trial(iRat,iSess).Latency_Lever(end+1:nTrials_expect,1) = NaN;
    scatter(iSess+(iSess-1), mean(trial_per_trial(iRat,iSess).Latency_Lever(idx)),  'r', 'filled');hold on
                    errorbar(iSess+(iSess-1), mean(trial_per_trial(iRat,iSess).Latency_Lever(idx)), ...
                    std(trial_per_trial(iRat,iSess).Latency_Lever(idx)) / power(numel(trial_per_trial(iRat,iSess).Latency_Lever(idx)),0.5), 'r')
                    scatter(iSess+(iSess-1)+1, mean(trial_per_trial(iRat,iSess).Latency_Lever(~idx)),  'g', 'filled')
                    errorbar(iSess+(iSess-1)+1, mean(trial_per_trial(iRat,iSess).Latency_Lever(~idx)), ...
                    std(trial_per_trial(iRat,iSess).Latency_Lever(~idx)) / power(numel(trial_per_trial(iRat,iSess).Latency_Lever(~idx)),0.5), 'g') 
    % Store Data
    latency_lever{1,1}(iRat,iSess) = mean(trial_per_trial(iRat,iSess).Latency_Lever(idx));
    latency_lever{1,2}(iRat,iSess) = mean(trial_per_trial(iRat,iSess).Latency_Lever(~idx));                 
    
    % Time Spent Close to Wall
    figure(c)
    if isequal(mod(iRat,2),true) % Plot UP for Partners
       subplot(2,1,1);
    else % Plot down for Partners
       subplot(2,1,2);
    end 
    trial_per_trial(iRat,iSess).TimeSpentCloseWall(end+1:nTrials_expect,1) = zeros;   
    trial_per_trial(iRat,iSess).TimeSpentCloseWall(isnan(trial_per_trial(iRat,iSess).TimeSpentCloseWall)) = 0;
    % Compute Percent of Time 
    to_plot = rdivide((trial_per_trial(iRat,iSess).TimeSpentCloseWall * 100), trial_per_trial(iRat,iSess).Trial_Dur);
    scatter(iSess+(iSess-1), mean(to_plot(idx)),  'r', 'filled');hold on
                    errorbar(iSess+(iSess-1), mean(to_plot(idx)), ...
                    std(to_plot(idx)) / power(numel(to_plot(idx)),0.5), 'r')
                    scatter(iSess+(iSess-1)+1, mean(to_plot(~idx)),  'g', 'filled')
                    errorbar(iSess+(iSess-1)+1, mean(to_plot(~idx)), ...
                    std(to_plot(~idx)) / power(numel(to_plot(~idx)),0.5), 'g') 
    % Store Data
    time_wall_act{1,1}(iRat,iSess) = mean(to_plot(idx));
    time_wall_act{1,2}(iRat,iSess) = mean(to_plot(~idx));                 
         
    % Time Spent Freezing
    trial_per_trial(iRat,iSess).TimeSpentFreezing(end+1:nTrials_expect,1) = NaN;      
    trial_per_trial(iRat,iSess).TimeSpentFreezing(isnan(trial_per_trial(iRat,iSess).TimeSpentFreezing)) = 0;      
    % Compute Percent of Time 
    figure(d)
    if isequal(mod(iRat,2),true) % Plot UP for Actors
       subplot(2,1,1);
    else % Plot down for Partners
       subplot(2,1,2);
    end 
    to_plot = rdivide((trial_per_trial(iRat,iSess).TimeSpentFreezing * 100), trial_per_trial(iRat,iSess).Trial_Dur);        
    scatter(iSess+(iSess-1), mean(to_plot(idx)),  'r', 'filled');hold on
                    errorbar(iSess+(iSess-1), mean(to_plot(idx)), ...
                    std(to_plot(idx)) / power(numel(to_plot(idx)),0.5), 'r')
                    scatter(iSess+(iSess-1)+1, mean(to_plot(~idx)),  'g', 'filled')
                    errorbar(iSess+(iSess-1)+1, mean(to_plot(~idx)), ...
                    std(to_plot(~idx)) / power(numel(to_plot(~idx)),0.5), 'g')                 
    % Store Data
    time_freezing_act{1,1}(iRat,iSess) = mean(to_plot(idx));
    time_freezing_act{1,2}(iRat,iSess) = mean(to_plot(~idx));  
    
    % Time Spent Rearing
    trial_per_trial(iRat,iSess).Rearing(end+1:nTrials_expect,1) = NaN;      
    trial_per_trial(iRat,iSess).Rearing(isnan(trial_per_trial(iRat,iSess).Rearing)) = 0;      
    % Compute Percent of Time 
    figure(b)
    if isequal(mod(iRat,2),true) % Plot UP for Actors
       subplot(2,1,1);
    else % Plot down for Partners
       subplot(2,1,2);
    end 
    to_plot = rdivide((trial_per_trial(iRat,iSess).Rearing * 100), trial_per_trial(iRat,iSess).Trial_Dur);        
    scatter(iSess+(iSess-1), mean(to_plot(idx)),  'r', 'filled');hold on
                    errorbar(iSess+(iSess-1), mean(to_plot(idx)), ...
                    std(to_plot(idx)) / power(numel(to_plot(idx)),0.5), 'r')
                    scatter(iSess+(iSess-1)+1, mean(to_plot(~idx)),  'g', 'filled')
                    errorbar(iSess+(iSess-1)+1, mean(to_plot(~idx)), ...
                    std(to_plot(~idx)) / power(numel(to_plot(~idx)),0.5), 'g')                 
    % Store Data
    time_rearing_act{1,1}(iRat,iSess) = mean(to_plot(idx));
    time_rearing_act{1,2}(iRat,iSess) = mean(to_plot(~idx)); 
    
end
break
end
end
    figure(a);
    title('Latency Lever')   
    xlabel('# Sessions'); ylabel('Latency Lever (s)');     
    
    figure(b); title('Latency Reward Poke')
    xlabel('# Sessions'); ylabel('Rwd Poke Duration (s)');    
    
    figure(c); subplot(2,1,1); title('Time Spent Close to Wall Actor')
    xlabel('# Sessions'); ylabel('Time Spent Close to Wall (s)');   
               subplot(2,1,1); title('Time Spent Close to Wall Partner ')
    xlabel('# Sessions'); ylabel('Time Spent Close to Wall (s)');
    
    figure(d); subplot(2,1,1); title('Time Spent Freezing Actor')
    xlabel('# Sessions'); ylabel('% Freezing');  
               subplot(2,1,1); title('Time Spent Freezing Partner')
    xlabel('# Sessions'); ylabel('% Freezing');  

end

