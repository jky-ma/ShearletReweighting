function D = getWaveletOperator(m,type,noScales)

    d = [m(1) m(2) 3*noScales+1];   
    W = opWavelet2(m(1),m(2),'Daubechies',type,noScales,1);

    D.adj           = @(x) (W*(x(:)));
    D.times         = @(y) (W'*y);
    
    [a,b,c]         = CoeffDims(m(1), m(2), noScales);
    temp            = size(D.adj(zeros(prod(m),1)),1)/(a*b);
    D.d             = [a,b,temp];
    
    D.spectra       = ones(m);
end

