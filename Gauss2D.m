function out=Gauss2D(x,sigma);
% Gauss2D	apply a 2 dimensional gauss filter
%
%    out = Gauss2D(x,sigma);
%
%    INPUT: x      image
%           sigma  of gauss filter [sigmaX sigmaY]
%
%    OUTPUT: out   filtered image
%



R = 1+ceil(3*sigma);   % cutoff radius of the gaussian kernel  
for i = -R:R,
   for j = -R:R,
      M(i+R+1,j+R+1) = exp(-(i*i+j*j)/2/sigma/sigma);
   end
end
M = M/sum(M(:));   % normalize the gaussian mask so that the sum is
                   % equal to 1

% Convolute matrices
%out=filter2(M,x);
out=imfilter(x,M,'replicate');

