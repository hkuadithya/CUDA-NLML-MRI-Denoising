function [output, smap]=NLMLbuild(input, D1) 
 
[m, n, z]=size(input);						% x, y and z dimension of the input image.
smap=zeros(m,n,z);
output=zeros(m,n,z);
options = optimset('Display','off', 'GradObj', 'on', 'LargeScale', 'off');
 

for k=1:z  
    for j=1:n 
        for i=1:m 
         if(input(i,j,k)==0) 
             continue; 
         end; 
         
         off=i-1 + (j-1) * m + (k-1) * m * n ;
         D2=D1(25*off+1 : 25*(off+1));
         M1=double(D2(D2>0));
         if(size(M1)==0)
             continue;
         end;
         R = fminunc(@(x)objective_functionChi1(M1(:),x),[mean(M1(:)) std(M1(:))],options);	% log likelihood of observations
         output(i,j,k)=R(1);	% finds a local minimum R(1) of the function objecttive_functionChi1
         smap(i,j,k)=R(2);		% finds the function value fval, i.e, R(2)
       end;
    end;
    fprintf('%d ', k);
end;
end
