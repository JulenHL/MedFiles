%# MedFiles
function [dat] = getMedFiles(file_idx,varargin)


%%             Düsseldorf University - 05.02.2013
%                   Julen Hernandez Lallement
% This function will extract the data from the Med PC logfile(s) into the
% structure dat with size(nSessions, nBoxes).
% Modifyable options are the size of the MedPC arrays
% (default Nx5) could be differently set by the user. In that case, change
% the nColumns parameter.

%% Netherlands Institute for Neuroscience
%  Julen Hernandez-Lallement - 19-10-2016

%% initialize defaults

% nColumns % number of columns entries in the diskvars; see MedPC script.

% Input is:
% file_idx = file name (excluding session number) that should be extracted
% varargin = folder of data of interest ("cd" in case folder should be scanned. DO
%            not enter anything if you want one single file to be extracted
% if you need to pick one file, enter the file_idx input only
nColumns = 12;

%% getting datafiles

% The Input must be a MedAssociates Operant Boxes Output Logfile or
% directory with logfiles

if isequal(nargin,2)
    if isdir(varargin{1})
        files = dir(varargin{1});
      % Update the last input to be more specific regarding screened files
%         files = files(strncmp('Baseline_S',{files(:).name},5)); % take
        files = files(strncmp(file_idx,{files(:).name},5)); % take only filenames specified in input
        datadir = varargin{1};
    elseif exist(varargin{1}, 'file')
        [datadir,fn,~] = fileparts(varargin{1});
        files = dir(datadir);
        files = files(strcmp(fn,{files(:).name}));
    else
        error('input is not a file nor a directory')
    end
elseif nargin > 2
    error('too many input arguments')
else
    [fn, datadir, outc] = uigetfile('', 'Select a MedAssociates Datafile');
    if isequal(outc, 0)
        error('no file selected')
    else
        files = dir(datadir);
        files = files(strcmp(fn,{files(:).name}));
        
    end
end

%% Read in of datafiles

number_readstring = repmat('%f', 1,nColumns);
dat = struct('filename', '',...
    'startdate', '',...
    'enddate', '',...
    'subject', '',...
    'experiment', '',...
    'group', '',...
    'boxnumber', '',...
    'starttime', '',...
    'endtime', '',...
    'scriptname', '');

nSessions = length(files);
% tmp_dat_empty = dat(1,1);

for iSess=1:nSessions
    
    iDat = 0;
    file_done = false;
    fn = {files(iSess).name};
    fn_str = fn{1}; 
    real_sess_numb = str2double(fn_str(strfind(fn_str, 'S') + 1: end)); 
    
    fid = fopen(fullfile(datadir,fn{1}));
    fseek(fid,0,'eof');
    end_pos = ftell(fid);
    frewind(fid);
    filename = textscan(fid, 'File: %s',1, 'delimiter', '\t');
    while ~file_done
        
        iDat = iDat + 1;
                
        dat(real_sess_numb,iDat).filename    = filename{1}{1};
        startdate                   = textscan(fid, 'Start Date: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).startdate   = startdate{1}{1};
        enddate                     = textscan(fid, 'End Date: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).enddate     = enddate{1}{1};
        subject                     = textscan(fid, 'Subject: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).subject     = subject{1};
        experiment                  = textscan(fid, 'Experiment: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).experiment  = experiment{1};
        group                       = textscan(fid, 'Group: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).group       = group{1};
        boxnumber                   = textscan(fid, 'Box: %f',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).boxnumber   = boxnumber{1};
        starttime                   = textscan(fid, 'Start Time: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).starttime   = starttime{1}{1};
        endtime                     = textscan(fid, 'End Time: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).endtime     = endtime{1}{1};
        scriptname                  = textscan(fid, 'MSN: %s',1, 'delimiter', '\t');
        dat(real_sess_numb,iDat).scriptname  = scriptname{1}{1};     
        
        reading_diskvars = true;
        disp('Start reading DiskVars')
        while reading_diskvars
            varname = textscan(fid, '%s',1, 'delimiter', '\t'); % ignore the var name
            if isempty(varname{1})
                reading_diskvars = false;
            else
                varname = varname{1}{1}(1);
                ok = true;
                out = [];
                linesread = 0;                               
            if ~isequal(isnan(str2double(varname)),true)
                error('Problem with Array Name Extraction. Check Colum Number')
            end
                
                while ok
                    curr_pos = ftell(fid);
                    % ignore the line indicator
                    tmp = textscan(fid, ['%*s' number_readstring],1,...
                        'MultipleDelimsAsOne', 1, 'CollectOutput', 1, 'ReturnOnError', true); 
                    linesread = linesread + 1;
                   % disp(['lines read for ' varname ': ' num2str(linesread)])
                    if all(isnan(tmp{1}))
                        ok = false;
                        new_pos = ftell(fid);
                        fseek(fid, curr_pos - new_pos, 'cof');
                    else
                        if any(isnan(tmp{1}))
                            ok = false;
                        end
                        out = [out,tmp{1}];
                    end
                    curr_pos = ftell(fid);
                    test = textscan(fid, '%s',1, 'MultipleDelimsAsOne', 1);
                    
                    if strncmp('Start',test{1},5)
                        ok = false;
                        reading_diskvars = false;
                    end
                    new_pos = ftell(fid);
                    fseek(fid, curr_pos - new_pos, 'cof');
                end
                dat(real_sess_numb, iDat).(varname) = out;
                curr_pos = ftell(fid);
            end
        end
        disp(['Finished reading DiskVars for boxdata: ' num2str(iDat)])
        if abs(curr_pos - end_pos)<10
            file_done = 1;
        end
    end       
end

% [extracted_data, rat_ID] = extract_data_MedFiles_v2(dat);

end# MedFiles
