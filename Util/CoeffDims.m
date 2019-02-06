function [pext, qext, levels] = CoeffDims(p, q, levels)
         
         if p >= 2^levels
            plevels = levels;
            if q >= 2^levels
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               qlevels = floor(log2(q));
               levels = min(plevels,qlevels);
               qext = ceil(q/(2^levels))*2^levels;
            else
               qext = q;
            end
            pext = ceil(p/(2^levels))*2^levels;
         elseif p > 1
            plevels = floor(log2(p));
            if q >= 2^levels
               levels = min(levels,plevels);
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               qlevels = floor(log2(q));
               levels = min(plevels,qlevels);
               qext = ceil(q/(2^levels))*2^levels;
            else
               levels = min(levels,plevels);
               qext = q;
            end
            pext = ceil(p/(2^levels))*2^levels;
         else
            pext = p;
            if q >= 2^levels
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               levels = floor(log2(q));
               qext = ceil(q/(2^levels))*2^levels;
            else
               qext = q;
            end
         end
end