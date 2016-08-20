function [noisy]=addNoise(ima, std) 

	noisy= round(abs(double(ima) + std * (randn(size(ima)) + sqrt(-1) * randn(size(ima)))));

end