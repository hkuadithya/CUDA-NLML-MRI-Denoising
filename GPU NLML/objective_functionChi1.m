function [f,df] = objective_functionChi1(m,x1)
    x=x1(1);
    sigma=x1(2);
    L=1; % Number of Coils
    [f,df] = logChipdfnew(m,x,sigma,L,[false true true]);
    f = -sum(f);
    df = -sum(df);
end