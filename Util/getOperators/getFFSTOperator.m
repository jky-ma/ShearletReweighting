function D = getFFSTOperator(m,noScales)
     temp   = zeros(m);
     vec    = @(x) x(:);
     array  = @(x,m) reshape(x,m);
     
     [ST, Psi] = shearletTransformSpect(temp,noScales);
     d         = size(ST);
     
     D.adj      = @(x) vec(shearletTransformSpect(array(x,m),Psi));
     D.times    = @(x) vec(inverseShearletTransformSpect(array(x,d),Psi));
     D.d        = d;
     D.spectra  = ones(m);
     D.Psi      = Psi;
end
