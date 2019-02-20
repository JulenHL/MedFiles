function [dat_signi, dat_raw, dat_ref_upper, dat_ref_down, indiv_scores, percent_choice_real] = permutation_group(varargin)

% This function will generate a reference distribution of nReps with
% randomly chosen data points from input data
% Individual scores are computed as percent change from baseline using the
% randomly selected data, i.e., virtual rats 
% then, each real data point (rat) is compared to the distribution and
% categorized as switchers or non switchers.

% Input is trial per trial choice data
% 2 = High Value Choice (i..e, + shock in shock sessions)
% 1 = Low Value Choice
% 0 = Non performed trial

% Rows should be rats
% Columns should be trials
% Concatenate sessions in one input array

%%

dat = varargin{1};

nRats = size(dat,1);
nTrials = size(dat,2);
nSess = 4;
nBlocks = 1;
nRep = 10000;
nSessUsed = [2:4];

disp (['Number of animals:  ' num2str(nRats)])
disp (['Number of trials:  ' num2str(nTrials)])
disp (['Number of sessions:  ' num2str(nSess)])
indiv_scores = NaN(nRep, 1,nRats);
percent_choice = NaN(nRep,4,nRats);
percent_choice_real = NaN(4,nRats);
dat_ref_upper = NaN(1,nRats);
dat_ref_down  = NaN(1,nRats); 

dat_raw = NaN(1,nRats);
dat_signi = NaN(1,nRats);
    
for iRat = 1 : nRats
    curr_dat = dat(iRat,:);

    for iRep = 1 : nRep

             progressBar(iRep/nRep, 50, ['Permuting...Rat # ' num2str(iRat)])     

        % Generate mask to select values to be mixed
            mask = randperm(nTrials);
        % Shuffle values
            curr_dat = curr_dat(mask);
        % Compute Percent preferences    
         for iSess = 1 : nSess;
             percent_choice(iRep,iSess, iRat) = numel(find(curr_dat(1,(iSess*(nTrials/nSess*nBlocks))-((nTrials/nSess*nBlocks)-1) : iSess*(nTrials/nSess*nBlocks)) == 2)) / ...
         numel(curr_dat(1,(iSess*(nTrials/nSess*nBlocks))-((nTrials/nSess*nBlocks)-1) : iSess*(nTrials/nSess*nBlocks))) * 100;
         end    
        % Compute invidivual scores as percent change scores
%             indiv_scores(iRep,1,iRat) = ((percent_choice(iRep,1, iRat) - mean(percent_choice(iRep,nSessUsed, iRat)))) / ...
%                                    ((percent_choice(iRep,1, iRat) + mean(percent_choice(iRep,nSessUsed, iRat))));
%             indiv_scores(iRep,1,iRat) = ((mean(percent_choice(iRep,nSessUsed, iRat)) - percent_choice(iRep,1, iRat))) / ...
%                                    ((percent_choice(iRep,1, iRat) + mean(percent_choice(iRep,nSessUsed, iRat))));
            indiv_scores(iRep,1,iRat) = ((mean(percent_choice(iRep,nSessUsed, iRat)) - percent_choice(iRep,1, iRat))) / ...
                                   ((100 - percent_choice(iRep,1, iRat)));
    %         indiv_scores(iRep,1) = ((curr_dat(1,iRep) - nanmean(curr_dat(2:PPP_sess,iRep)))) / ...
    %                                ((nanmean(curr_dat(2:PPP_sess,iRep)))) * 100;
    end
    
curr_indiv_scores = sort(indiv_scores(:,1,iRat));
         % Reference data
     dat_ref_upper(1,iRat) = curr_indiv_scores(ceil(size(curr_indiv_scores,1)*0.975),1);
     dat_ref_down(1,iRat) = curr_indiv_scores(ceil(size(curr_indiv_scores,1)*0.025),1);
% dat_ref_upper(1,iRat) = curr_indiv_scores(ceil(size(curr_indiv_scores,1)*0.95),1);
% dat_ref_down(1,iRat)  = curr_indiv_scores(ceil(size(curr_indiv_scores,1)*0.05),1);

     % Compute real scores
    
      for iSess = 1 : nSess;
             percent_choice_real(iSess, iRat) =  numel(find(dat(iRat,(iSess*(nTrials/nSess*nBlocks))-((nTrials/nSess*nBlocks)-1) : iSess*(nTrials/nSess*nBlocks)) == 2)) / ...
         numel(dat(iRat,(iSess*(nTrials/nSess*nBlocks))-((nTrials/nSess*nBlocks)-1) : iSess*(nTrials/nSess*nBlocks))) * 100;                 
      end    
      
%      dat_raw(1,iRat) = ((percent_choice_real(1, iRat) - mean(percent_choice_real(nSessUsed, iRat)))) / ...
%                        ((percent_choice_real(1, iRat) + mean(percent_choice_real(nSessUsed, iRat))));
%      dat_raw(1,iRat) = ((mean(percent_choice_real(nSessUsed, iRat)) - percent_choice_real(1, iRat))) / ...
%                        ((percent_choice_real(1, iRat) + mean(percent_choice_real(nSessUsed, iRat))));
     dat_raw(1,iRat) = ((mean(percent_choice_real(nSessUsed, iRat)) - percent_choice_real(1, iRat))) / ...
                       ((100- percent_choice_real(1, iRat)));
%      dat_raw(1,iRat) = ((dat(iRat,1) - nanmean(dat(iRat,2:PPP_sess)))) / ...
%                        ((nanmean(dat(iRat,2:PPP_sess)))) * 100;                   
     dat_signi(:,iRat) = dat_raw(1,iRat) >  dat_ref_upper(1,iRat); % >= ?
     
end
     
% end

% figure;hist(indiv_scores(:))
% line([0 0],ylim, 'Color', 'r')


figure;
for iRat = 1 : nRats;
    subplot(4,7,iRat); scatter(1, dat_raw(1,iRat), 'filled', 'k'); hold on
line([1 1],[dat_ref_down(1,iRat) dat_ref_upper(1,iRat)], 'Color', 'r')
end
figure; steps = 100; h = max(hist(dat_raw, steps));
hist(dat_raw, steps)
hold on
% line([dat_ref_down dat_ref_upper], [h h], 'Color', 'r')

iqr0 = prctile(dat_raw,[25 75]);
    
plot(nanmean(dat_raw), h+1,'ko','markerfacecolor','c', ...
                           'markersize',10)
plot([nanmean(dat_raw)-nanstd(dat_raw),...
                           nanmean(dat_raw)+nanstd(dat_raw)],...
                           [h+1,h+1] ,'-c', 'linewidth', 4)
plot(nanmedian(dat_raw),h+2,'ko','markerfacecolor','r',...
                           'markersize',10)
plot(iqr0,[h+2, h+2], '-r', 'linewidth',4);

limy = ylim;
ylim([limy(1) limy(2) + 2])


figure
scatter(ones(nRats,1)+ rand(nRats,1), dat_raw, 'filled', 'k')
xlim([0 5])
hold on

plot([4.5,4.5],[nanmean(dat_raw) - nanstd(dat_raw), nanmean(dat_raw) + nanstd(dat_raw)], '-c', 'linewidth', 4)
plot(4.5, nanmean(dat_raw),'ko','markerfacecolor','c', 'markersize',10)

iqr0 = prctile(dat_raw,[25 75]);
plot(4.2, nanmedian(dat_raw),'ko','markerfacecolor','r','markersize',10)
plot([4.2,4.2], iqr0, '-r', 'linewidth',4);

% plot([2.5,2.5],[nanmean(dat_raw_SBL041) - nanstd(dat_raw_SBL041), nanmean(dat_raw_SBL041) + nanstd(dat_raw_SBL041)], '-c', 'linewidth', 4)
% plot(2.5, nanmean(dat_raw_SBL041),'ko','markerfacecolor','c', 'markersize',10)
% 
% iqr0 = prctile(dat_raw_SBL041,[25 75]);
% plot(2.2, nanmedian(dat_raw_SBL041),'ko','markerfacecolor','r','markersize',10)
% plot([2.2,2.2], iqr0, '-r', 'linewidth',4);

