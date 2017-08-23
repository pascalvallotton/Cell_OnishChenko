function [ imOUt ] = bfimread( path, numb )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imbf=bfopen(path);
imOUt=imbf{1,1}{numb,1};%imshow(imOUt,[])

%input_args=
%output_args=
end

