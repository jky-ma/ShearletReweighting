function out = TVOnlysolver(y,m,A,D,alpha,beta,mu1,mu2,varargin)
    %% SB/ADMM 
    %  code to solve the following CS-Problem in Analysis Formulation:
    % 
    %   (1)         min_u J(u)  s.t. ||A*u - y||_2^2<=sigma,
    % 
    % where J(u) = lambda_iter*|D*u|_1 +  alpha*||grad*u||_1,
    % for a possibly redundant Frame D and A the forward measurement
    % operator. lambda_iter can be updated iteratively.
    %
    %
    % To solve (1), we transform it into an unconstrained formulation
    % 
    %   (2)        min_u J(u)  + beta/2*||A*u - y||_2^2,
    %
    % Problem (2) is solved by using ADMM iterations. 
    % That is we are "decpoupling the l1 and l2 portions of the problem" by substituting
    %
    %        w = |D*u|_1, dx = |grad_x*u|_1, dy = |grad_y*u|_1
    % 
    % and considering the problem
    %
    %   (3)     (u_(k+1),w_(k+1),dx_(k+1),dy_(k+1)) =  argmin_{u,w,dx,dy}  lambda*|w|_1 + alpha*|(dx,dy)|_1 + ...
    %                   beta/2*|y - A*u + y|_2^2 + mu1/2*|dx - grad_x*u - bx_k|_2^2 + mu1/2*|dy - grad_y*u - by_k|_2^2 ...
    %                       + mu2/2*|w - D*u - bw_k|_2^2
    %              
    % with updates
    %
    %           bw_(k+1) = bw_k + (D*u_(k+1) - w_(k+1))
    %           bx_(k+1) = bx_k + (grad_x*u_(k+1) - dx_(k+1))
    %           by_(k+1) = by_k + (grad_y*u_(k+1) - dy_(k+1))
    %
    % Equation (3) is thereby solved separately for each component. For w,dx, and dy this is
    % done by an (possibly adaptive) shrinkage formula, and for u we solve a linear system,
    % which is (sometimes) diagonalized by fft2.
    %
    %% Input:
    % y         -   measured signal
    % m         -   dimensions, m = (m(1),m(2))
    % A         -   measurements contain A.adj(), A.times() and A.mask,
    %                   set A.mask = NaN if mask unknown
    % D         -   used dictionary, contains D.adj() -> analysis operator
    %                                   and D.times() -> synthesis operator
    %                                   D.d dimension of analysis coefficients, D.spectra
    % alpha     -   reg. parameter (TV)
    % beta      -   reg. parameter (data-fit)
    % mu1       -   reg. parameter (ADMM - TV-term)
    % mu2       -   reg. parameter (ADMM - analysis-term)
    %
    %% Optional input:
    % cMap          -   colormap 
    % doPlot        -   flag for plotting
    % doReport      -   flag for reporting
    % correct       -   tool to encourge real valued reconstruction in
    %                       shrinkage of analysis operator (for Fourier
    %                       measurements)
    % f             -   true signal (if known)
    % zo            -   zoom-part for plotting
    % maxIter  -   number of iterations
    %
    % written by Maximilian März & Jackie Ma, February 2017
    

    %% define technical stuff 
    vec            = @(x)   x(:);    
    array          = @(x,m) reshape(x,m);
    FourierOp      = opDFT2(m(1),m(2),true);
    F.times        = @(x) FourierOp*x;  % fft2, unitary
    F.adj          = @(x) FourierOp'*x; % ifft2
    d              = D.d;
    
    %% set optional parameters
    cMap            =   gray;
    doReport        =   true;
    doTrack         =   false;
    doPlot          =   true;
    correct         =   @(x) (x);
    adaptive        =   'L1';
    f               =   NaN;
    zo{1}           =   m(1)/2:m(1)/2+40;
    zo{2}           =   m(2)/2:m(2)/2+40;
    maxIter         =   25;
    lambda          =   1;
    dynRange        =   1;
    epsilon         =   1e-4;
    Normalize       =   false;
    innerIter       =   2;
    % overwrite them if possible
    for k=1:2:length(varargin),     % overwrites default parameter
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
    % tv only
    lambda          = 0;
    mu2             = 0;
     
    % fetch forward gradient
    [gradX,gradY,Dx,Dy] = getPeriodicGradForw(m);
    % get corresponding Laplacian
    [L,Lx,Ly,DL,DLx,DLy] = getPeriodicLaplacian(m);
    DL = fftshift(DL);
    
    %% check if A is given via mask and normalize measured data
    if ~isnan(A.mask)
        iterative    = false;
        % normalize measured data
        normFac      = max(abs(F.adj(vec(A.mask.*array(y,m)))));
        y            = y./normFac;
        % define measurement operator 
        A.times      = @(x) vec(A.mask.*array(F.times(x),m));
        A.adj        = @(x) F.adj(vec(A.mask.*array(x,m)));
        
        % diagonal matrix for u-subproblem
        K = beta*A.mask + mu1*DL; %  + mu2*D.spectra;        
    else
        iterative   = true;
        % normalize 
        normFac     = max(abs(A.adj(y)));
        y           = y./normFac;
        
        % matrix for u-subproblem
        afun        = @(x) beta*(A.adj(A.times(x))) + mu1*(L*x); % + mu2*vec(D.times(D.adj(x)));       
    end
    
    %% check if f exist, take A.adj(y) otherwise, define truC
    if isnan(f);        
        f    = array(A.adj(y),m)*normFac;
        truC = D.adj(f);
    else
        truC = D.adj(f/normFac);
    end
    
    
    %% initialize for ADMM    
    u           = A.adj(y);
    yk          = y;
    w           = zeros(prod(sum(d,1)),1);
    dx          = zeros(prod(m),1);
    dy          = zeros(prod(m),1);
    bw          = zeros(prod(sum(d,1)),1);
    bx          = zeros(prod(m),1);
    by          = zeros(prod(m),1);
    rnsp        = 0;
    rorre       = 0;
    rorreler    = 0;
    miss        = 0;
    Iter        = 0;
    uOld        = Inf;
    
    if doReport
        fprintf('------------------ CS with SB/ADMM --------------------------\n\n')
        fprintf('Solve min_u lambda_ITER*|D*u|_1 + alpha*|grad u|_1, s.t. ||A*u - y||_2^2 < s^2\n\n')
        fprintf('using SB/ADMM iterations. Regularization parameters:\n\n')
        fprintf('beta for data-fit: beta=%d\n',beta)
        fprintf('alpha for TV-norm: alpha=%d\n',alpha)
        fprintf('lambda for Besov: lambda=%d\n',lambda)
    	fprintf('mu1, mu2 for SB/ADMM: mu1=%d, mu2=%d\n\n',mu1,mu2)
        fprintf('Maximal number of iterations: %d\n',maxIter)
        
        DataMiss    = norm(A.times(u) - y);
        TV          = norm(gradX*u,1) + norm(gradY*u,1); % anisotropic(!)
        Besov       = 0;% norm(D.adj(u),1);
        
        fprintf('Begin of iterations:\n\n')
        fprintf('It No.    |  DataMissFit  | Besov \t   |   TV   \t  | |u-uOld| \t   | Time for Iter\n')
        fprintf('-------------------------------------------------------------------------------------------------\n')
        fprintf('  %2d     | %2.4e\t  | %2.4e\t   | %2.4e   |   %2.4e   | \t %2.2e      \n',0,DataMiss,Besov,TV,0,0);
    end

    while (Iter < maxIter)
           
        tOuter    = tic;
        Iter      = Iter + 1;
     
            for i=1:innerIter           
                %% step 1
                % solve u subproblem
                uOld    = u;
                rhs     = beta*(A.adj(yk)) + mu1*gradX'*(dx-bx) + mu1*gradY'*(dy-by); %... 
                                        % + mu2*D.times(w-bw);   
                % solve linear system with pointwise division or iterative
                % solver
                if iterative
                     [u,~] = pcg(afun,rhs,1e-10,75,[],[],uOld);
                else
                     u = F.adj(vec(array(F.times(rhs),m)./K));
                end
                
                %% step 2
                % updates in w, dx and dy
                % [w,info]   = DictionaryShrinking(D.adj(u),bw,d,1/(mu2),truC,adaptive,lambda,epsilon);
                % w          = correct(w);
                
                % isotropic dx, dy update
                temp        = zeros(m(1),m(2),2);
                temp(:,:,1) = array(gradX*u + bx,m);
                temp(:,:,2) = array(gradY*u + by,m);
                temp2       = IsotropicShrink(temp,alpha/mu1);
                dx          = vec(temp2(:,:,1)); 
                dy          = vec(temp2(:,:,2));         
            end
            
            %% step 3 
            % updates in bw, bx and by
            % bw = bw + (D.adj(u) - w);
            bx = bx + (gradX*u - dx);
            by = by + (gradY*u - dy);

        if doPlot
            %% plot f
            subplot(2,3,1)
            imagesc(abs(f))
            colormap(cMap)
            axis('image')
            c = caxis;
            colorbar()
            title(sprintf('original image f'));

            %% zoomed u
            subplot(2,3,2)
            pl = array(u,m)*normFac;
            if Normalize
                pl = (pl-min(pl(:)))/(max(pl(:))-min(pl(:)));
            end
            zu = abs(pl(zo{1},zo{2}));
            imagesc(zu);
            axis('image');
            colormap(cMap)
            title('zoomed u');

            %% plot gray
            subplot(2,3,3)
            imagesc(abs(pl))
            axis('image') 
            % caxis(c)
            colorbar()
            colormap(cMap)
            title(sprintf('reconstructed image u,\n PSNR=%f,SSIM=%f',psnr(f./max(abs(f(:))),pl./max(abs(pl(:)))),ssim(abs(f),abs(pl),'DynamicRange',dynRange)));
            %% warning: ssim has to be adapted to the dynamic range of signal

            %% plot difference
            subplot(2,3,4)
            diff = f -pl;
            imagesc(abs(diff))
            axis('image')
            colorbar()
            title(sprintf('difference f - u, RE=%f',norm(vec((f))-vec((pl)))/norm(vec(f))));

            %% plot measured data
            subplot(2,3,5)
            g = A.adj(y);
            if ~isnan(A.mask)
                r = norm(double(A.mask(:)),1);
                imagesc(abs(array(g,m))*normFac)
                axis('image')
                colorbar()
                title(sprintf('measured data,\n undersampling rate=%1.2f',r/prod(m)))  
            else
                imagesc(abs(array(g,m))*normFac)
                axis('image')
                colorbar()
                title(sprintf('measured data'))
            end


            %% analyse adaptive truncation
%             subplot(2,3,6)
%             if norm(info) ~= 0
%                 info(:,1) = info(:,1);
%                 plot(info(2:end-1,:));
%                 legend({'mean of truncated signal','mean of threshold','mean of u'});
%                 title(sprintf('Truncation'));  
%             end
            drawnow
        end    
        
        % Bregman update in y
        yk = yk + 1*(y - A.times(u));
            
        
        if doReport
            DataMiss    = norm(A.times(u) - y);
            TV          = norm(gradX*u,1) + norm(gradY*u,1);
            Besov       = 0; % norm(D.adj(u),1);
            tOuter      = toc(tOuter);
            fprintf('  %2d    |  %2.4e\t   | %2.4e   |   %2.4e   |   %2.4e  | \t %2.2e      \n',Iter,DataMiss,Besov,TV,norm(u-uOld)/norm(u),tOuter);
        end
        
        
        if doTrack
            rnsp(Iter)      = psnr(f./max(abs(f(:))), array(u,m)./max(abs(u(:))));
            rorre(Iter) 	= norm(f(:) - normFac*u(:))./norm(f(:));
            rorreler(Iter)  = norm(uOld(:) - u(:))./norm(uOld(:));
            miss(Iter)      = ssim(abs(f),abs(normFac*reshape(u,m)),'DynamicRange',dynRange);
        end
    end

    out.rec         = array(u,m).*normFac;
    % out.info        = info;
    out.psnr        = rnsp;
    out.error       = rorre;
    out.relerror    = rorreler;
    out.ssim        = miss;
    
end
    
