function out =  TGVOnlysolver(y,m,A,D,alpha,beta,mu,varargin)
    %% SB/ADMM 
    %  code to solve the following CS-Problem in Analysis Formulation:
    % 
    %   (1)         min_u J(u)  s.t. ||A*u - y||_2^2<=s^2,
    % 
    % where J(u) = lambda*|D*u|_1 + TGV_alpha^2(u), 
    % for a redundant shearlet Frame D, and A being the measurement operator.
    %
    % To solve (1), we transform it into an unconstrained formulation
    % 
    %   (2)        min_u J(u)  + beta/2*||A*u - y||_2^2,
    %
    % Problem (2) is solved by using ADMM iterations. 
    % For details on the implementation see TVSolver.m and the related paper.
    % 
    %% Input:
    % y         -   measured signal
    % m         -   dimensions, m = (m(1),m(2))
    % A         -   measurements contain A.adj(), A.times() and A.mask,
    %                   set A.mask = NaN if mask unknown
    % D         -   used dictionary, contains D.adj() -> analysis operator
    %                                   and D.times() -> synthesis operator
    %                                   D.d dimension of analysis coefficients, D.spectra
    % lambda    -   reg. parameter Besov
    % alpha     -   (alpha_1,alpha_0), reg. parameter (TGV)
    % beta      -   reg. parameter (data-fit)
    % mu        -   (mu_1,mu_2,mu_3), reg. parameter (ADMM) 
    %
    %% Optional input:
    % cMap      -   colormap 
    % doPlot    -   flag for plotting
    % doReport  -   flag for reporting
    % correct   -   tool to encourge real valued solutions
    % f         -   true signal (if known)
    %
    % written by Maximilian März & Jackie Ma, February 2017

    %% set optional parameters
    cMap            =   gray;
    doReport        =   true;
    doPlot          =   true;
    doTrack         =   false;
    correct         =   @(x) (x);
    adaptive        =   'L1';
    f               =   NaN;
    zo{1}           =   m(1)/2:m(1)/2+40;
    zo{2}           =   m(2)/2:m(2)/2+40;
    maxIter    =   20;
    lambda          =   1;
    dynRange        =   1;
    iterative       = true;
    epsilon         =   1e-4;
    Normalize       = false;
    innerIter       = 2;
    
    % overwrite them if possible
    for k=1:2:length(varargin),     % overwrites default parameter
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
    % TGV only
    lambda  = 0;
    mu(1)   = 0;


    %% define operators for ADMM
    % get periodic backward gradient
    [BgradX,BgradY,BDx,BDy] = getPeriodicGradBack(m);
    % get periodic forward gradient
    [FgradX,FgradY,FDx,FDy] = getPeriodicGradForw(m);
    
    % get corresponding Laplacian
    [L,Lx,Ly,DL,DLx,DLy] = getPeriodicLaplacian(m);
    DL      = fftshift(DL);
    DLx     = fftshift(DLx);
    DLy     = fftshift(DLy);
    FDx     = fftshift(FDx);
    FDy     = fftshift(FDy);
    BDx     = fftshift(BDx);
    BDy     = fftshift(BDy);
    % operator for 2nd derivatives in TGV
    E              = @(p) getE(p,m,BgradX,BgradY);
    
    % Fourier operator
    FourierOp      = opDFT2(m(1),m(2),true);
    F.times        = @(x) FourierOp*x;  % fft2, unitary
    F.adj          = @(x) FourierOp'*x; % ifft2

    % reshaping oprators            
    vec     = @(x) x(:);
    array   = @(x,m) reshape(x,m);
    % arrayzed versions of D and F
    dSize   = D.d;
    Dt      = @(x) array(D.adj(vec(x)),sum(dSize,1));
    Df      = @(x) array(D.times(vec(x)),m);
    Ff      = @(x) array(F.times(vec(x)),m);
    Ft      = @(x) array(F.adj(vec(x)),m);
    
    
    %% check if A is given via mask and normalize measured data
    if ~isnan(A.mask)
        iterative    = false;
        % normalize 
        normFac      = max(abs(F.adj(vec(A.mask.*array(y,m)))));
        y            = y./normFac;
        % define measuremnt operator 
        A.times      = @(x) vec(A.mask.*array(F.times(x),m));
        A.adj        = @(x) F.adj(vec(A.mask.*array(x,m)));
               
    else
        iterative   = true;
        % normalize 
        normFac     = max(abs(A.adj(y)));
        y           = y./normFac;
    end
    
    %% check if f exist, take A.adj(y) otherwise, define truC
    if isnan(f);
        f    = array(A.adj(y),m)*normFac;
        truC = D.adj(f/normFac);
    else
        truC = D.adj(f/normFac);
    end
    

    
    %% initialize
    yk      = y;
    w       = zeros(sum(dSize,1));
    bw      = zeros(sum(dSize,1));
    d       = zeros(m(1),m(2),2);
    bd      = zeros(m(1),m(2),2);
    t       = zeros(m(1),m(2),4);
    bt      = zeros(m(1),m(2),4);
    u       = A.adj(y);
    Du      = zeros(m(1),m(2),2);
    p       = zeros(m(1),m(2),2);
    Ep      = E(p);
    rnsp    = 0;
    rorre   = 0;
    miss    = 0;

    Iter    = 0;
    uOld        = u;
    DataMiss    = norm(vec(u) - vec(f));
    TGV         = norm(vec(Ep(:,:,1)),1) + 2*norm(vec(Ep(:,:,2)),1) + norm(vec(Ep(:,:,4)),1) ...
                            + norm(vec(Du-p),1);
    Besov       = 0; % norm(D.adj(u),1);


    %% construct (diagonal) blocks of linear system
    d2 = mu(2)  + mu(3)*(DLx + 0.5*DLy);
    d3 = mu(2)  + mu(3)*(0.5*DLx + DLy);

    d4 = -1*mu(2)*FDx;
    d5 = -1*mu(2)*FDy;
    d6 =  mu(3)*0.5*conj(BDx).*BDy;
    
    % pre-calculated part of upper-left entry of linear system
    XX   = (-d5.*conj(d5)./d3  - (d4-(d5.*conj(d6)./d3).*(conj(d4) - d6.*conj(d5)./d3))./(d2-(d6.*conj(d6))./d3));
    
    if iterative
        % define upper-left block without using a mask
        d1   = @(x) beta*F.times((A.adj(A.times(F.adj(x))))) + mu(2)*vec(DL.*array(x,m)) + mu(1)*vec(D.spectra.*array(x,m));
        afun = @(x) d1(x) + vec(XX.*reshape(x,m)); % upper-left entry
        X    = {afun, d1, 0;...
                    (d4 - d5.*conj(d6)./d3), (d2 -d6.*conj(d6)./d3), 0};
    else
        d1   = beta*A.mask + mu(2)*DL;
        X    = {(d1-(d5.*conj(d5))./d3  - (d4-(d5.*conj(d6))./d3).*(conj(d4) - d6.*conj(d5)./d3)./(d2-(d6.*conj(d6))./d3)), 0, 0;...
                d4 - d5.*conj(d6)./d3, (d2 - (d6.*conj(d6))./d3), 0};      
    end
    dstruct  = {d1,d2,d3,d4,d5,d6};

    if doReport
        fprintf('------------------ CS with SB/ADMM --------------------------\n\n')
        fprintf('Solve min_u gamma*|D*u|_1 + TGV_alpha^2(u) s.t. ||A*u - y||_2^2 < s^2\n\n')
        fprintf('using SB/ADMM iterations. Regularization parameters:\n\n')
        fprintf('beta for data-fit: beta=%d\n',beta)
        fprintf('lambda for Besov: lambda=%d\n',lambda)
        fprintf('alpha_1 for TGV-norm: alpha1=%d\n',alpha(1))
        fprintf('alpha_0 for TGV-norm: alpha2=%d\n',alpha(2))
        fprintf('mu for SB/ADMM: mu1=%d, mu2=%d, mu3 =%d\n\n',mu(1),mu(2),mu(3))
        fprintf('Maximal number of outer iterations: %d\n',maxIter)

        fprintf('Begin of iterations:\n\n')
        fprintf('   It No. |  DataMissFit  | Besov \t   |    TGV   \t   |   |u-uOld| \t |  Time for  Iter\n')
        fprintf('-------------------------------------------------------------------------------------------------\n')
        fprintf('  %2d      | %2.2e\t  | %2.2e\t    | %2.2e      |   %2.2e      | \t %2.2e      \n',0,DataMiss,Besov,TGV,0,0);
    end
    
    while (Iter < maxIter)  

        tOuter  = tic;
        Iter    = Iter + 1;

            for i=1:innerIter
                %% step 1
                % update in z = (u, px, py)  
                % construct RHS
                B1  = beta*array(A.adj(yk),m) + ...
                            mu(2)*array((FgradX'*vec((d(:,:,1) - bd(:,:,1))) + FgradY'*vec(d(:,:,2) - bd(:,:,2))),m);
                FB1 = Ff(B1);

                B2  = mu(2)*(bd(:,:,1) - d(:,:,1)) + ...
                            mu(3)*array(BgradX'*vec(t(:,:,1) - bt(:,:,1)) +  BgradY'*vec(t(:,:,2)-bt(:,:,2)),m);
                FB2 = Ff(B2);

                B3  = mu(2)*(bd(:,:,2) - d(:,:,2)) +    ...
                            mu(3)*array(BgradX'*vec(t(:,:,2) - bt(:,:,2)) + BgradY'*vec(t(:,:,4) - bt(:,:,4)),m);
                FB3 = Ff(B3);

                uOld = u;
                % solve linear system
                [u,p(:,:,1),p(:,:,2)] = (solveSystem(dstruct,X,FB1,FB2,FB3,iterative,uOld));
                

                u         = (Ft(u));         
                Du(:,:,1) = (array(FgradX*vec(u),m));
                Du(:,:,2) = (array(FgradY*vec(u),m));
                p(:,:,1)  = Ft(p(:,:,1));
                p(:,:,2)  = Ft(p(:,:,2));
                Ep        = E(p);
               
                %% step 2
                % update in w 
                % [w,info]   = DictionaryShrinking(D.adj(u),bw,D.d,1/mu(1),truC,adaptive,lambda,epsilon);
                % w          = correct(array(w,dSize));
                
                
                % update in d 
                d = IsotropicShrink(Du - p + bd,alpha(1)/mu(2));

                % update in t
                t = IsotropicShrink(Ep + bt,alpha(2)/mu(3));
            end

            %% step 3
            % Bregman updates
            % bw = bw + (Dt(u) - w);
            bd = bd + (Du - p - d);
            bt = bt + (Ep - t);


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
            pl = array(u,m);
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
            %caxis(c)
            colorbar()
            colormap(cMap)
            title(sprintf('reconstructed image u,\n PSNR=%f,SSIM=%f',psnr(f./max(abs(f(:))),pl./max(abs(pl(:)))),ssim(abs(f),abs(pl),'DynamicRange',dynRange)));

            %% plot difference
            subplot(2,3,4)
            diff = f -pl;
            imagesc(abs(diff./norm(vec(f))))
            axis('image')
            colorbar()
            title(sprintf('difference f - u, RE=%f',norm(vec((f))-vec((pl)))/norm(vec((f)))));

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
        yk = yk + 1*(y - A.times(vec(u)));

        
        if doReport
            DataMiss    = norm(vec(u) - vec(f));
            Besov       = 0; % norm(D.adj(u),1);
            TGV         = norm(vec(Ep(:,:,1)),1) + 2*norm(vec(Ep(:,:,2)),1) + norm(vec(Ep(:,:,4)),1) ...
                            + norm(vec(Du-p),1);          
            tOuter = toc(tOuter);
            fprintf('  %2d      |  %2.4e\t | %2.4e\t   | %2.4e     |   %2.4e     | \t %2.2e      \n',Iter,DataMiss,Besov,TGV,norm(vec(u)-vec(uOld)),tOuter);
        end
        
        if doTrack
            rnsp(Iter)      = psnr(f./max(abs(f(:))), array(u,m)./max(abs(u(:))));
            rorre(Iter) 	= norm(vec((f))-vec((u*normFac)))/norm(vec((f)));
            miss(Iter)      = ssim(abs(f),abs(array(u,m)*normFac),'DynamicRange',dynRange);
        end
    end
    
    out.rec     = array(u,m).*normFac;
    % out.info    = info;
    out.psnr    = rnsp;
    out.error   = rorre;
    out.ssim    = miss;
end

