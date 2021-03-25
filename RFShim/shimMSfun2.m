function [rfOut,errTot,phs,errRMSE,errSAR] = shimMSfun2(A,phs,beta,betaCtr,fdFlag,faThresh,maxItr,mask,B1p)

        if maxItr == 1
            fdThresh = faThresh;
            maxItr = 20;
        else
            minThresh = 0;
            faIncr = (faThresh-minThresh) / (maxItr-1); 
            fdThresh = faThresh : -faIncr : minThresh;
        end
        
        if isempty(A)
            rfOut = zeros(8,1); % if empty slice, dont calculate in cost.
            errTot(idx) = 0;
            errRMSE(idx) = 0;
            errSAR(idx) = 0;

        else
            
            errTmpTot = [];
            errTmpRMSE = [];
            errTmpSAR = [];

            itr = 1;
            errTmpTot(itr) = inf;
            errTmpRMSE(itr) = inf;
            errTmpSAR(itr) = inf;

            if fdFlag                   
                % Use TV to avoid nulls
                TV0 = FD;
                ctrFlag = 1;

                %%%%%%%%%%%%%%%%%%

                TV0 = TV0 * (B1p .* mask);

                %%%%%%%%%%%%%%%%%%

                TV0 = permute(TV0,[3,1,2,4]);
                TV0 = permute(TV0(:,:),[2,1]);
            end

            flag = 1;
            while flag % for target phase update
                itr = itr + 1;

                % Use TV to avoid nulls
                TV = 0;
                if fdFlag && itr > 2                      
                    tmp = abs(mtmp);

                    % Important
                    if length(fdThresh) == 1
                        % optimized
                        tmp = (tmp < fdThresh) & (tmp>0);
                        TV = sqrt( 1.5 * sqrt(sum(tmp(:))) ) * TV0;
                    else
                        % general recommended
                        tmp = (tmp < fdThresh(ctrFlag)) & (tmp>0);
                        TV = sqrt( 10 * sqrt(sum(tmp(:))) ) * TV0;
                    end    
                end

                tmpRF = (A'*A+beta(betaCtr)*speye(8)+TV'*TV)\(A'*exp(1i*phs));

                tmpRF = tmpRF / mean(abs(A*tmpRF));
                tmpRF = tmpRF * exp(-1i*angle(tmpRF(1)));
            
                % calculate the excitation pattern
                mtmp = A*tmpRF;

                % update the target phase pattern
                phs = angle(mtmp);

                errTmpRMSE(itr) = norm(abs(mtmp)-1)^2;
                errTmpSAR(itr) = real(tmpRF'*tmpRF);
                errTmpTot(itr) = errTmpRMSE(itr) + beta(betaCtr) * errTmpSAR(itr);

                if (errTmpTot(itr) < 0.99999*errTmpTot(itr-1))
                    tmpRFOld = tmpRF;
                else
                    % if TV triggered, then restart the rf update
                    if (sum(TV(:)) ~= 0) && ctrFlag <= maxItr
                        flag = 1; ctrFlag = ctrFlag + 1;

                        itr = 1;
                        errTmpTot = [];
                        errTmpRMSE = [];
                        errTmpSAR = [];

                        errTmpTot(itr) = inf;
                        errTmpRMSE(itr) = inf;
                        errTmpSAR(itr) = inf;
                    else
                        flag = 0;
                        if (errTmpTot(itr) <= errTmpTot(itr-1))
                            rfOut = tmpRF;
                            errTot = errTmpTot(itr);
                            errRMSE = errTmpRMSE(itr);
                            errSAR = errTmpSAR(itr);
                        else
                            rfOut = tmpRFOld;
                            errTot = errTmpTot(itr-1);
                            errRMSE = errTmpRMSE(itr-1);
                            errSAR = errTmpSAR(itr-1);
                        end
                    end                     
                end
            end
        end
end