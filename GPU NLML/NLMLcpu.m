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
 
F = padarray(input,[f f f],'symmetric');	% symmetric padding across the three dimensions.
D1=diffneigh1_cpu(single(F), [t, t, t], [f, f, f]);

end
