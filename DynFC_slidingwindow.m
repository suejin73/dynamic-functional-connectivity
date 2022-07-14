function [NV, ICdyn, WinFC] = DynFC_slidingwindow(data,WL)

%%% The function calculates dynamic functional connectivity features
%%% using a sliding window approach with pairwise correlation.
%%% Although it accounts for different length of fMRI timecourses (if filled with 0),
%%% it is suggested to use the same length of fMRI timeserious.
%%% If there are any, a warning message will show up.
%%% For IC dynamics, the number of ROIs has to be even and the ROI arrangement sould be 
%%% [Left ROIs, followed by Right ROIs] in the same order or vice versa [R,L]. 
%%%
%%% input:
%%%     data is an M*N*T matrix of fMRI signals (in double),
%%%     M is the number of ROIs
%%%     N is the number of subjects
%%%     T is the fMRI volumes
%%%     data can also be an 1*N cell,
%%%     each cell element contains T (volume) * M (ROI) signals (in double)
%%%     WL is the length of windown (default = 20)
%%%
%%% output:
%%%     NV is network variation: the average FC differences across windows of global matrix
%%%     ICdyn is inter-hemispheric dynamics: the average FC differences across windows of IC (i.e. L/R ROIs) pairs
%%%     WinFC is a cellarray of the windowed correlation matrices
%%%
%%% ref (part of the features used in): 
%%%     Education, and the balance between dynamic and stationary functional connectivity jointly support executive functions 
%%%     in relapsing-remitting multiple sclerosis. Hum Brain Mapp. 2018 Dec;39(12):5039-5049.

%%% Sue-Jin Lin @ MNI, 20200226
%%% 20200430 update to take cell as well
%%% 20220714 update to calculate IC dynamics

%%% decide data format
checkdata = iscell(data);
if checkdata ==1
    % cell to matrix
    for ii = 1:size(data,2)
        dataMAT(:,:,ii) = data{ii}';
    end
    clear data
    data = permute(dataMAT, [1 3 2]);
end


%%% set up parameters for sliding window
if ~exist('WL','var')
    WL = 20;
end


%%% calculate windowed correlation matrices
for i=1:size(data,2)  % per subject
    signals = squeeze(data(:,i,:))';
    for k=1:size(signals,1)-WL  % per window
        [r,p] = corrcoef(signals([k:WL+k-1],:)) ; % move 1 timepoint
        CorrMat = r - diag(diag(r));
        WinFC{i}(:,:,k) = CorrMat; % WinFC in a cell arraw, each cell has windowed FC matrices per sub
        clear r p CorrMat
    end
end

%%% calculate dynamic features
if mod(size(data,1),2) ~= 0
    disp('warning! # of ROI is not even, skip calculating IC dynamics')
    ICdyn = [];
    
    % only calculate NV
    for m = 1:length(WinFC) % per sub's windowed FC
        WinMat = WinFC{m};
        WinMat = atanh(WinMat); % z transform
        for l= 1:size(WinMat,3)-1 % per window
            temp1 = sqrt((triu(WinMat(:,:,l)) - triu(WinMat(:,:,l+1))).^2);
            NV_all(1,l) = nansum(nansum(temp1));
        end
        if sum(isinf(NV_all))>0
            sprintf('the %d th subject has different length of fMRI timepoints',m)
            EndTPind = find(isinf(NV_all));
            NV_all(isinf(NV_all)) = 0; % if there are 0 in the timepoints-> generate inf, make it 0
            % adjust NV measure for shorter time course
            NV_ave = nansum(NV_all)/EndTPind; % average across windows per subject
            NV(m,1) = double(NV_ave); % contains 1 value of all subejcts
        else
            NV_ave = nansum(NV_all)/size(WinMat,3)-1; % average across windows per subject
            NV(m,1) = double(NV_ave); % contains 1 value of all subejcts
        end
    end
    
    
else
    
    % calculate all features
    for m = 1:length(WinFC) % per sub's windowed FC
        WinMat = WinFC{m};
        WinMat = atanh(WinMat); % z transform
        for l= 1:size(WinMat,3)-1 % per window
            temp1 = sqrt((triu(WinMat(:,:,l)) - triu(WinMat(:,:,l+1))).^2);
            NV_all(1,l) = nansum(nansum(temp1));
            temp2 = sqrt((diag(WinMat([1:size(data,1)/2], [size(data,1)/2+1:size(data,1)], l)) - diag(WinMat([1:size(data,1)/2], [size(data,1)/2+1:size(data,1)], l+1))).^2);
            ICdyn_all(1,l) = nansum(nansum(temp2));
        end
        if sum(isinf(NV_all))>0
            sprintf('the %d th subject has different length of fMRI timepoints',m)
            EndTPind = find(isinf(NV_all));
            NV_all(isinf(NV_all)) = 0; % if there are 0 in the timepoints-> generate inf, make it 0
            ICdyn_all(isinf(ICdyn_all)) = 0;
            % adjust NV measure for shorter time course
            NV_ave = nansum(NV_all)/EndTPind; % average across windows per subject
            NV(m,1) = double(NV_ave); % contains 1 value of all subejcts
            ICdyn_ave = nansum(ICdyn_all)/EndTPind;
            ICdyn(m,1) = double(ICdyn_ave);
        else
            NV_ave = nansum(NV_all)/size(WinMat,3)-1; % average across windows per subject
            NV(m,1) = double(NV_ave); % contains 1 value of all subejcts
            ICdyn_ave = nansum(ICdyn_all)/size(WinMat,3)-1;
            ICdyn(m,1) = double(ICdyn_ave);
        end
    end
end

disp('DONE!! ^o^')

