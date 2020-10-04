%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     2D MS RF Shimming Code
%     Zhipeng Cao, Ph.D
%     Vanderbilt University Medical Center, 2020.09.27
% 
%     Note:
%     1) It interfaces the MRCodeTool "External Waveform Definition".
%     2) It outputs the shim ".dat" file for patch to read (separately
%     provided by Zhipeng) on Philips 7T.
%     3) It needs Fessler's Image Recon Toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear all; clc; close all;

    load('b1data.mat');
    [Nx,Ny,Nslice,Nc] = size(B1p);
      
    %% Algorithm Settings
    betaMin = 0.1; % or 0
    betaMax = 5;
    betaIncr = 0.1;
    
    %% Init with Quad Phase
    phsinit = angle(sum(B1p,3));

    %% Initialize
    betaCtr0 = uint8((betaMax-betaMin)/betaIncr + 1);
    
    % Prepare error report
    errOut_shim1 = zeros(betaCtr0,1);
    errOut_shim_RMSE1 = zeros(betaCtr0,1);
    errOut_shim_SAR1 = zeros(betaCtr0,1);
    
    errOut_shim2 = zeros(betaCtr0,1);
    errOut_shim_RMSE2 = zeros(betaCtr0,1);
    errOut_shim_SAR2 = zeros(betaCtr0,1);
    
    %  Build system matrix for each 2D slice
    tmp = permute(B1p,[3 1 2]);
    tmp = tmp(:,:).';
    tmp = tmp(logical(mask),:);
    A = bsxfun(@times,exp(-1i*angle(tmp(:,1))),tmp);
    
    % Init RF matrix
    rf1 = zeros(8,betaCtr0);
    
    % This gives slightly lower SAR !!
    beta = betaMax : -betaIncr : betaMin; % reg high to low
    
    %% Run RF Shim
    for betaCtr = 1:betaCtr0        
        %  Init Ainv with current beta
        if betaCtr == 1 % Important!
            tmp = phsinit;
            phs = tmp(mask);
        end
        
        [tmpRF,err,phs,errRMSE,errSAR] = shimMSfun1(A,phs,beta,betaCtr);
        
        rf1(:,betaCtr) = tmpRF;
        
        errOut_shim_RMSE1(betaCtr,:) = errRMSE;
        errOut_shim_SAR1(betaCtr,:) = errSAR;
        errOut_shim1(betaCtr,:) = err;
    end
        
    figure(1);
    hold on;
    plot(sum(errOut_shim_SAR1,2), sum(errOut_shim_RMSE1,2),'-*');
    hold off;
    title('L-curves');
    xlabel('RF Power');
    ylabel('FA Accuracy');
    
    %% Evaluate Shim Result and Redo Slices with FD Regularization
    
    % Change this to select different solutions on the L-curve, the idx = 1 start with beta = 0, thus highest power and best homogeneity.
    idxShim = betaCtr0 * 0.8; % reg high to low

    rf_out = rf1(:,idxShim);

    % Adjust the overall flip angle of each slice for the output solution
    rf_out = rf_out / mean(abs(A*rf_out));
    
    % Show Conventional Shim
    tmpRF = abs(rf_out) .* exp(1i*round(angle(rf_out)/pi*180)/180*pi);
        
    tmp = A*tmpRF;
    mcompp2 = embed(tmp,squeeze(mask));
    figure(2); im(fliplr(mcompp2)); caxis([0,1.2]); colorbar; title('w/o FD');
        
    %% Null detection
    b1Threshold = 0.15;
    
    % Evaluate Null
    fdFlag = zeros(1,2); % second index is similar to spokes

    if ~isempty(A)
        tmpM = abs(A*rf_out);
        fdFlag = (sum((tmpM<b1Threshold) & (tmpM>0))>1);
    end
    
    % Init RF matrix
    rf2 = zeros(8,betaCtr0);
    
    % Run RF Shim Again, with fdFlag
    for betaCtr = 1:betaCtr0        
        %  Init Ainv with current beta
        if betaCtr == 1 % Option
            tmp = phsinit;
            phs = tmp(mask);
        end

        [tmpRF,err,phs,errRMSE,errSAR] = shimMSfun2(A,phs,beta,betaCtr,fdFlag,mask,B1p);
        
        rf2(:,betaCtr) = tmpRF;
        
        errOut_shim_RMSE2(betaCtr,:) = errRMSE;
        errOut_shim_SAR2(betaCtr,:) = errSAR;
        errOut_shim2(betaCtr,:) = err;
    end
    
    figure(1);
    hold on;
    plot(sum(errOut_shim_SAR2,2), sum(errOut_shim_RMSE2,2),'-*');
    hold off;
    title('RF Shimming');
    xlabel('RF Power');
    ylabel('FA Accuracy');
    legend({'RF Shim MLS','RF Shim FD-MLS'});
    
    %% Quad Mode Comparison - scaled fair comparison
    m = A*(ones(8,1)/mean(abs(A*ones(8,1))));
    mcomp = embed(m,squeeze(mask));
    figure(4); im(fliplr(mcomp)); caxis([0,1.2]); colorbar; title('Quad');
    
    %% Select Shim Result for export and compare
    rf_out = rf2(:,idxShim);
    
    % Adjust the overall flip angle of each slice for the output solution
    rf_out = rf_out / mean(abs(A*rf_out));
    
    %% Shim Mode Comparison    
    tmpRF = abs(rf_out) .* exp(1i*round(angle(rf_out)/pi*180)/180*pi);
        
    tmp = A*tmpRF;
    mcompp2 = embed(tmp,squeeze(mask));

    figure(3); im(fliplr(mcompp2)); caxis([0,1.2]); colorbar; title('w/ FD');