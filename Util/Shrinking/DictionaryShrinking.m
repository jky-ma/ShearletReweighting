function [s,info] = DictionaryShrinking(x,y,d,alpha,truC,adaptive,lambda,epsilon)
    %% Shrinkage-Operator for optimization subproblem in ADMM
    % Shrinks analysis coefficients adaptively for each subband
    %
    %% Input    x       -   analysis coefficients
    %           y       -   2nd part of signal
    %           d       -   dimensions of coefficient space
    %           alpha   -   shrinking parameter
    %           truC    -   analysis coefficients of the true signal
    %           adaptive-   string to select shrinking method
    %           lambda  -   regularization parameter

    %
    %% Output    s       -   shrinked coefficients
    %            info    -   vector conating information about e.g. threshold
    %
    % written by Jackie Ma & Maximilian Maerz, April 2016
    
    %% initialize and define functions
    vec     = @(x) x(:);
    shrink  = @(x,alpha) (sign(x)).*max(abs(x)-alpha, 0);%(x./abs(x)).*max(abs(x)-alpha, 0);
    % shrink  = @(x,alpha) x.*max(abs(x)-alpha, 0)./(max(abs(x)-alpha, 0) + alpha);
    x       = reshape(x,sum(d,1));
    y       = reshape(y,sum(d,1));
    truC    = reshape(truC,d);
    s       = x + y;

    
    
    if strcmp(adaptive,'TrueIRL1')
        for i=2:d(3)
            a = (epsilon + abs(truC(:,:,i))).^(-1);
            s(:,:,i) = shrink(x(:,:,i) + y(:,:,i),a*lambda*alpha);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(a*lambda*alpha));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
    elseif strcmp(adaptive,'IRL1')
        factor = 1e-2.*Factor;
        for i=1:d(3)
            a = (epsilon + abs(x(:,:,i)+y(:,:,i))).^(-1);
            s(:,:,i) = shrink(x(:,:,i) + y(:,:,i),a*factor*alpha);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(a*factor*alpha));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
    elseif strcmp(adaptive,'NewTrueIRL1')
        for i=1:d(3)
            lam(i) = quantile(vec(abs(truC(:,:,i))),0.98);
            a = lam(i)*(epsilon + abs(truC(:,:,i))).^(-1);
            s(:,:,i) = shrink(x(:,:,i) + y(:,:,i),a*lambda*alpha);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(a*lambda*alpha));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
    elseif strcmp(adaptive,'NewIRL1')
        for i=1:d(3) 
            if i == 1
                lam(i) = 0;
            else
                lam(i) = max(vec(abs(x(:,:,i)+y(:,:,i))));
            end
            a = 1e-1*lam(i)*(epsilon + abs(x(:,:,i)+y(:,:,i))).^(-1);

            s(:,:,i)  = shrink(x(:,:,i) + y(:,:,i),a*alpha*lambda);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(a*alpha*lambda));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
        
     elseif strcmp(adaptive,'NewIRL1CT')
        for i=1:d(3) 
            if i == 1
                lam(i) = max(vec(abs(x(:,:,i)+y(:,:,i))))*1e-5;
            else
                lam(i) = quantile(vec(abs(x(:,:,i))),0.98);%+y(:,:,i)))) deleted for CT;
            end
            a = lam(i)*(epsilon + abs(x(:,:,i))).^(-1); %+y(:,:,i)))) deleted for CT;

            s(:,:,i)  = shrink(x(:,:,i) + y(:,:,i),a*alpha*lambda);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(a*alpha*lambda));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end

    elseif strcmp(adaptive,'L1')
        for i=1:d(3)
            if i ==1
                lambdas(i) = 1e0*lambda;
            else
                lambdas(i) = 1*lambda;
            end
            s(:,:,i)  = shrink(x(:,:,i) + y(:,:,i),lambdas(i)*alpha);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(lambdas(i)*alpha));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
    elseif strcmp(adaptive,'NewL1')
        for i=2:d(3)
            lambdas(i) = quantile(vec(abs(x(:,:,i))),0.98);
            s(:,:,i)  = shrink(x(:,:,i) + y(:,:,i),lambdas(i)*alpha*lambda);
            info(i,1) = mean(vec(abs(x(:,:,i) + y(:,:,i))));
            info(i,2) = mean(vec(lambdas(i)*alpha*lambda));
            info(i,3) = mean(vec(abs(s(:,:,i))));
        end
    elseif strcmp(adaptive,'TV')  
        s       = s*0;
        info    = 0;
    end

    %% output
    s = s(:);
end