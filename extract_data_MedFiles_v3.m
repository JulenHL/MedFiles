function [extracted_data, rat_ID, rwd_pos] = extract_data_MedFiles_v3(varargin)

%% This function will extract the data from Trial per Trial, Block, Event, Time Event array.

% Input:
% This should be the output array from getMedFiles_v2

% Output
% extracted_data is the reorganized trial per trial & block per block data
% rat_ID is the order of rats trained (rows are sessions)
% rwd_pos is the position (left or right) of the high value reward for each
% rat (columns are rats, reordered)


% The script creates a LogFile ('output_log.txt') containing all information about missed
% trials / unperformed sessions / inconsistencies, per rat, per session

% \  A(0)     = Side of the High Effort Reward (1 = Right, 2 = Left)
% \  A()      = Trial by Trial Data
% \  A(P)     = Trial Number
% \  A(P+1)   = Trial Number in this Block
% \  A(P+2)   = Trial Type (1 = Free, 2 = Forced Left, 3 = Forced Right)
% \  A(P+3)   = Chosen Side (1 = Left, 2 = Right)
% \  A(P+4)   = Chosen Value (1 = Low Value, 2 = High Value)
% \  A(P+5)   = Omission Type (-1, -2, -3)
% \  A(P+6)   = Trial Completed Successfully (1=Yes, 0=No)
% \  A(P+7)   = Latency Nose Poke
% \  A(P+8)   = Latency Lever Press
% \  A(P+9)   = Latency Reward Poke
% \  A(P+10)  = Duration Reward Poke
% \  A(P+11)  = Short Nose Poke duration (1 = Yes, 0 = No)
% 
% \  B()      = Block by Block Data
% \  B(W)     = Block Number
% \  B(W+1)   = Number of High Value Reward Chosen in this Block
% \  B(W+2)   = Number of Low  Value Reward Chosen in this Block
% \  B(W+3)   = Number of Successful Free Trials this Block
% \  B(W+4)   = Number of Nose Poke   Omission     (Type -1)
% \  B(W+5)   = Number of Lever Press Omission     (Type -2)
% \  B(W+6)   = Number of Failed Forced Trials     (Type -3) 


%% Script

global output_log

nSess            = size(varargin{1,1},1);
nRats            = size(varargin{1,1},2);
nTrials          = 24;
nTrialsPerBlock  = 24;
nBlocks          = nTrials / nTrialsPerBlock;
nForced          = 4;
size_trial_array = 12; % Columns number in A Array
size_block_array = 7;  % Columns number in B Array


extracted_data = struct; 

% Pre- allocation
A = cell(nSess, nRats);
B = cell(nSess, nRats);



reward_position = NaN(nSess,nRats); % Whether High Value was Left or Right
rat_ID = NaN(nSess,nRats);

for iSess = 1 : nSess;
    
    sum_rat = 0;
%     rat_ID = [];

plot_choice     = NaN(nTrials,nRats);                            % Produces cumsum dat for Choice Data
plot_side       = NaN(nTrials,nRats);                            % Produces cumsum dat for Side   Data
absolut_choice  = NaN(nTrials,nRats);                            % Produces Absolut Choice Data
absolut_side    = NaN(nTrials,nRats);                            % Produces Absolut Side   Data
trial_per_trial_all = NaN(size_trial_array,nTrials,nRats);
block_per_block_all = NaN(size_block_array,nBlocks,nRats);
    
   for iBox = 1 : nRats;
        
A{iSess,iBox} = varargin{1,1}(iSess,iBox).A;
A{iSess,iBox} = A{iSess,iBox}(2:(size_trial_array*nTrials)+1);

B{iSess,iBox} = varargin{1,1}(iSess,iBox).B;
B{iSess,iBox} = B{iSess,iBox}(2:(size_block_array*nBlocks)+1);

% size_trial_array = findpattern(A{iSess,iBox}, [2 2]) - 1; % Columns number in A Array. Caution: This is specific to my data storage.
% size_block_array = 7;  % Columns number in B Array. Caution: This is specific to my data storage.

% Extract Position of High Value Reward
reward_position(iSess,iBox) = varargin{1,1}(iSess,iBox).A(1,1);  
iRat = str2double(varargin{1,1}(iSess,iBox).subject{1});
rat_ID(iSess, iBox) = iRat;
    
% Check if user has correctly entered the subjects number  
    if isequal(iBox,nRats);       
        for idx_rat = 1 : nRats;
            sum_rat = sum_rat + (nRats - (nRats-idx_rat));
        end
        
        if sum(rat_ID(iSess, :)) ~= sum_rat
         cprintf('-err', ['Subject(s) number not consistent - Look Sess #' num2str(iSess)]);disp(' ')         
%        disp(['Subject(s) number not consistent - Look Sess #' num2str(iSess)])
%         error (['Subjects are: ' num2str(rat_ID(iSess,:))])
        end       
    end
    
    trial_per_trial_all(:,:,iRat) = reshape(A{iSess,iBox},size_trial_array,nTrials);
    block_per_block_all(:,:,iRat) = reshape(B{iSess,iBox},size_block_array,nBlocks);
    
%     if ~isequal(cumsum(trial_per_trial_all(1,:,iRat)),cumsum(1:size(trial_per_trial_all(1,:,iRat) > 0,2)))
%         error('trial per trial array is not correctly organized. Check Column number')
%     end
                 
    % Double Check Data by Cross Referencing Arrays A & B
    if ~isequal(size(find(trial_per_trial_all(5,nForced + 1 : end,iRat) == 2),2),sum(block_per_block_all(2,:,iRat))) % Compares High Value Choice
       disp(['Discrepancy in data between Arrays A & B in High Value Choices for Rat #' num2str(iRat) ' at Session #' num2str(iSess)])
    end
    
    if ~isequal(size(find(trial_per_trial_all(5,nForced + 1 : end,iRat) == 1),2),sum(block_per_block_all(3,:,iRat))) % Compares Low Value Choice
       disp(['Discrepancy in data between Arrays A & B in Low Value Choices for Rat #' num2str(iRat) ' at Session #' num2str(iSess)])
    end
    
    extracted_data(iSess, 1).trial = trial_per_trial_all;
    extracted_data(iSess, 1).block = block_per_block_all;

   u = 0;
   v = 0;
   temp = trial_per_trial_all(:,:,iRat);
   
   for i = 1 : nTrials;
       
      if temp(5,i) == 1;                          % Low Value Reward Chosen
          absolut_choice(i,iRat) = 1;
             plot_choice(i,iRat) = u + 1;         % +1 For Low Value Reward
             u = u + 1;
         if reward_position(iSess,iBox) == 1;     % Right Lever is High Value Reward  
             absolut_side(i,iRat) = 2;            % 2 For Left Side
             plot_side(i,iRat) = v - 1;           % -1 For Right Side
             v = v - 1;   
         elseif reward_position(iSess,iBox) == 2; % Left Lever is High Value Reward 
             absolut_side(i,iRat) = 1;            % 1 For Right Side             
             plot_side(i,iRat) = v + 1;           % +1 For Left Side
             v = v + 1;
             
         elseif reward_position(iSess,iBox) == 0; % No Defined High Value Rwd Position. Should not happen ;)                 
                 output_log = (['Inconsistency in LogFile at session # ' num2str(iSess) 'Wrong position of Rwd']);
                 disp (output_log)
                 write_output  
                 
                 try
         plot_choice(i,iRat) = plot_choice(i-1,iRat);
         plot_side(i,iRat)   = plot_side(i-1,iRat);
         u = u -1;
         absolut_choice(i,iRat)  = 0;
                 catch
         plot_choice(i,iRat) = 0;
         plot_side(i,iRat) = 0;    
         u = u -1;
         absolut_choice(i,iRat)  = 0;
                end
                 
         else output_log = (['Unexpected problem: check LogFile Trial # ' num2str(i) ' at session # ' num2str(iSess)]);
              disp (output_log)
              write_output  
         end
          
 
       elseif temp(5,i) == 2;                   % High Value Reward Chosen        
            plot_choice(i,iRat) = u - 1;          % -1 For High Value Reward
            u = u - 1;
            absolut_choice(i,iRat) = 2;
         if reward_position(iSess,iBox) == 2 ;    % Left Lever is High Value Reward 
             absolut_side(i,iRat) = 2;            % 2 For Left Side                   
             plot_side(i,iRat) = v - 1;           % -1 For Left Side
             v = v - 1;
             
         elseif reward_position(iSess,iBox) == 1; % Right Lever is High Value Reward
             absolut_side(i,iRat) = 1;            % 1 For Right Side            
             plot_side(i,iRat) = v + 1;           % +1 for left side
             v = v + 1;
             
        elseif temp(5,i) == 0;
                 output_log = (['Inconsistency in LogFile at session # ' num2str(iSess) 'Wrong position of Rwd']);

                 write_output  
                 
              try
       plot_choice(i,iRat) = plot_choice(i-1,iRat);
       plot_side(i,iRat) = plot_side(i-1,iRat);
       u = u -1;
       absolut_choice(i,iRat)  = 0;
              catch
       plot_choice(i,iRat) = 0;
       plot_side(i,iRat) = 0;   
       u = u -1;
       absolut_choice(i,iRat)  = 0;
    
              end
        
         else    output_log = (['Unexpected problem: check LogFile Trial # ' num2str(i) ' at session # ' num2str(iSess)]);
                 write_output      
         end
   
      elseif temp(6,i) < 0;   % Omission recorded
                 output_log = (['Omission for Rat # ' num2str(iRat) ' at Trial # ' num2str(i) ' at session # ' num2str(iSess)]);
                 write_output 
                 
          try
   plot_choice(i,iRat) = plot_choice(i-1,iRat);
      plot_side(i,iRat) = plot_side(i-1,iRat);
      absolut_choice(i,iRat)  = 0;
          catch
   plot_choice(i,iRat) = 0;
       plot_side(i,iRat) = 0;       
       absolut_choice(i,iRat)  = 0;
          end
          
      elseif temp(1,i) == 0;  % Not performed trial         
             output_log = (['Rat #' num2str(iRat) ' stopped at trial # ' num2str(i) ' at session # ' num2str(iSess)]);
             write_output  
             
          try
       plot_choice(i,iRat) = plot_choice(i-1,iRat);
       plot_side(i,iRat) = plot_side(i-1,iRat);
       absolut_choice(i,iRat)  = 0;
          catch
       plot_choice(i,iRat) = 0;
       plot_side(i,iRat) = 0;       
       absolut_choice(i,iRat)  = 0;
          end
          
      else
          disp('Choice data not correctly collected. Maybe session stopped by experimentator or wrong trial per trial arraz size')
          try
       plot_choice(i,iRat) = plot_choice(i-1,iRat);
       plot_side(i,iRat) = plot_side(i-1,iRat);
       absolut_choice(i,iRat)  = 0;
          catch
       plot_choice(i,iRat) = 0;
       plot_side(i,iRat) = 0;       
       absolut_choice(i,iRat)  = 0;
          end
   
      end
                 
   end
 
         if isequal(sum(isnan(plot_choice(:,iBox))),nTrials) % Session was not performed
             
             plot_choice(length(plot_choice)-sum(isnan(plot_choice(:,iBox))) + 1 : ...
              length(plot_choice),iBox) = 0;
          
            plot_side(length(plot_side)-sum(isnan(plot_side(:,iBox))) + 1 : ...
              length(plot_side),iBox) =  0;
             
             
         elseif sum(isnan(plot_choice(:,iBox))) > 0 % some trials were not performed
              
             plot_choice(length(plot_choice)-sum(isnan(plot_choice(:,iBox))) + 1 : ...
              length(plot_choice),iBox) =  plot_choice(length(plot_choice)-...
              sum(isnan(plot_choice(:,iBox))), iBox);
          
            plot_side(length(plot_side)-sum(isnan(plot_side(:,iBox))) + 1 : ...
              length(plot_side),iBox) =  plot_side(length(plot_side)-...
              sum(isnan(plot_side(:,iBox))), iBox);
          end
    end
    
    extracted_data (iSess, 1).cumsumchoice     = plot_choice;
    extracted_data (iSess, 1).cumsumside       = plot_side;
    extracted_data (iSess, 1).absolut_choice   = absolut_choice; 
    extracted_data (iSess, 1).absolut_side     = absolut_side; 
end

% Check Look for side of HV reward
[~,idx_rwd_pos] = sort(rat_ID');
idx_rwd_pos = idx_rwd_pos';

for iSess = 1 : nSess;
curr_dat = reward_position(iSess,:);
rwd_pos(iSess,:) = curr_dat(idx_rwd_pos(iSess,:));
end

rwd_pos = [sort(rat_ID(1,:));rwd_pos];

for iRat = 1 : nRats;
    if size(unique(rwd_pos(2:end,iRat)),1) > 1
       cprintf('-err', ['Problem with HV side allocation for subject #' num2str(iRat)]);disp(' ')
    end
end

    
    function write_output
        
                 disp(output_log)
                 fid = fopen('output_log.txt','a+');
                 fprintf(fid, ' %s\r\n', output_log);
                 fclose(fid);        
    end

end