function D = getShearletOperator(m,ndir)

    if length(m) == 2
        if length(ndir) == 1
            Sys         = SLgetShearletSystem2D(0,m(1),m(2),(ndir));
        else
            Sys         = SLgetShearletSystem2D(0,m(1),m(2),length(ndir),ndir);
        end

        D.d         = size(Sys.shearlets);
        D.spectra   = (sum(abs(Sys.shearlets).^2,3));

        array       = @(x,m) reshape(x,m);
        vec         = @(x) x(:);

        D.adj       = @(x) vec(SLsheardec2Dswapped(array(x,m),Sys));
        D.times     = @(x) vec(SLshearrec2Dswapped(array(x,D.d),Sys));


        % swapp elemtns in Sys
        t1              = Sys.RMS(1);
        Sys.RMS(1)      = Sys.RMS(end);
        Sys.RMS(end)    = t1;    
        t2              = Sys.shearletIdxs(1,:);
        Sys.shearletIdxs(1,:) = Sys.shearletIdxs(end,:);
        Sys.shearletIdxs(end,:) = t2;
        D.Sys = Sys;
    elseif length(m) == 3
        if length(ndir) == 1
            Sys         = SLgetShearletSystem3D(0,m(1),m(2),m(3),(ndir));
        else
            Sys         = SLgetShearletSystem3D(0,m(1),m(2),m(3),length(ndir),ndir);
        end
        D.d             = size(Sys.shearlets);
        D.spectra       = sum(abs(Sys.shearlets).^2,4);
        
        array           = @(x,m) reshape(x,m);
        vec             = @(x) x(:);
        
        D.adj           = @(x)  vec(SLsheardec3D(array(x,m),Sys));  
        D.times         = @(y)  vec(SLshearrec3D(array(y,D.d),Sys));
    end
end

function y =  SLsheardec2Dswapped(x,Sys)
   y            = SLsheardec2D(x,Sys);
   temp1        = y(:,:,1);
   y(:,:,1)     = y(:,:,end);
   y(:,:,end)   = temp1;
end

function y =  SLshearrec2Dswapped(x,Sys) 
   temp1        = x(:,:,1);
   x(:,:,1)     = x(:,:,end);
   x(:,:,end)   = temp1;
   y            = SLshearrec2D(x,Sys);
end