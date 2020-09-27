function [rfOut,errTot,phs,errRMSE,errSAR] = shimMSfun1(A,phs,beta,betaCtr)

        % MS Shim Optimization
        rfOut = zeros(8,1);
                        
        if isempty(A)
            rfOut = zeros(8,1); % if empty slice, dont calculate in cost.
            errTot = 0;
            errRMSE = 0;
            errSAR = 0;

        else

            % new
            errTmpTot = [];
            errTmpRMSE = [];
            errTmpSAR = [];

            itr = 1;
            errTmpTot(itr) = inf;
            errTmpRMSE(itr) = inf;
            errTmpSAR(itr) = inf;

            flag = 1;
            while flag % for target phase update
                itr = itr + 1;

                % Turns out the initial guess is not important since
                % there is effect from target phase update, thus pinverse is good enough. (?)

                tmpRF = (A'*A+beta(betaCtr)*speye(8))\(A'*exp(1i*phs));

                % calculate the excitation pattern
                m = A*tmpRF;

                % update the target phase pattern
                phs = angle(m);

                errTmpRMSE(itr) = 1/2*norm(m-exp(1i*phs))^2;
                errTmpSAR(itr) = 1/2*real(tmpRF'*tmpRF);
                errTmpTot(itr) = errTmpRMSE(itr) + beta(betaCtr) * errTmpSAR(itr);

                if (errTmpTot(itr) < 0.99999*errTmpTot(itr-1))
                    tmpRFOld = tmpRF;
                else
                    flag = 0;
                    if (errTmpTot(itr) <= errTmpTot(itr-1))
                        rfOut(:,1) = tmpRF;
                        errTot = errTmpTot(itr);
                        errRMSE = errTmpRMSE(itr);
                        errSAR = errTmpSAR(itr);
                    else
                        rfOut(:,1) = tmpRFOld;
                        errTot = errTmpTot(itr-1);
                        errRMSE = errTmpRMSE(itr-1);
                        errSAR = errTmpSAR(itr-1);
                    end
                end
            end
        end
end