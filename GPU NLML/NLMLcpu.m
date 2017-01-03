function [D1]=NLMLcpu(input, t, f) 
%***************************************************************************************************
 %Non local ML method for denoising single and multicoil MRI
 %Ref : J Rajan et. al.," Nonlocal maximum likelihood estimation method for denoising multiple-coil 
 %      magnetic resonance images", Magnetic Resonance Imaging, 2012.  % 
 % Usage : D = NLML2012(input,t,f)
 % eg :    D = NLML2012(input,5,3)   
 % D - Denoised image, input - Noisy image, t - search window, f - local
 % window;
 % Set the numbeer of coils N in the function objective_functionChi1.m
 %**************************************************************************************************
 
 %[m n z]=size(input);						% x, y and z dimension of the input image.
 %smap=zeros(m,n,z); 						% output image= zeroes in the three dimensions.
 F = padarray(input,[f f f],'symmetric');	% symmetric padding across the three dimensions.
 %output=zeros(m,n,z);
 %options = optimset('Display','off', 'GradObj', 'on', 'LargeScale', 'off');
 
D1=diffneigh1_cpu(single(F), [t, t, t], [f, f, f]);

%{
for k=1:z  
    for j=1:n 
        parfor i=1:m 
         if(input(i,j,k)==0) 
             continue; 
         end; 
         
         off=i-1 + (j-1) * m + (k-1) * m * n ;
         D2=D1(25*off+1 : 25*(off+1));
         M1=double(D2(D2>0));
     
         R = fminunc(@(x)objective_functionChi1(M1(:),x),[mean(M1(:)) std(M1(:))],options);	% log likelihood of observations
         output(i,j,k)=R(1);	% finds a local minimum R(1) of the function objecttive_functionChi1
         smap(i,j,k)=R(2);		% finds the function value fval, i.e, R(2)
       end;
       j
    end;
end;
%}
end