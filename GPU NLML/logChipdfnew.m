function [P,dP] = logChipdfnew(x, A, sigma, L, selGrHess)
%
% Returns the logArithm of the rice probAbility density of the measurement X, 
% with noise free signal amplituded A And a  noise level of sigma. 


% Jelle VerAArt, University of Antwerp, 9-11-2011


if isscalar(sigma)
    sigma = sigma * ones(size(A));
end;

Arg = (A .* x) ./ (sigma.^2);
besL1 = besselmx(double('I'), L-1, Arg, 1);
P = log(x.^L) - log(sigma.^2) + log(A.^(1-L)) - (A.^2 + x.^2)./(2*sigma.^2) + log(besL1)+abs(real(Arg));


if nargout>1

  
    if nargin<5
        selGrHess = true(1,3);
    elseif ~islogical(selGrHess) || numel(selGrHess)~=3
        error('incorrect selGrHess, should be 3 element logical vector.');
    end
    
    besL2 = besselmx(double('I'), L-2, Arg, 1);
    besL = besselmx(double('I'), L, Arg, 1);
    besRat = 0.5*(besL2 + besL)./besL1;
    
    
    if selGrHess(2)
        dPdA = (1-L)./A - A./(sigma.^2) + besRat.*x./(sigma.^2);   
    else
        dPdA = [];
    end
    if selGrHess(1)
        dPdx = L./x - x./(sigma.^2) + besRat.*A./(sigma.^2);
    else
        dPdx = [];
    end;
    if selGrHess(3)
        dpdsigma = -2./sigma + (A.^2 + x.^2)./(sigma.^3) - 2*besRat.*A.*x./(sigma.^3);
    else
        dpdsigma = [];
    end;
    dP = [dPdx(:) dPdA(:) dpdsigma(:)];
end;
