function [ S ] = FFT_PSD( x,y,T )
%FFT_PSD computes the PSD of a set of sampled data utilizing fft()
%   Detailed explanation goes here

xdft = fft(x);
xdft = xdft(1:length(x)/2);         % "folding" causes mirror image
xdft(2:end-1) = 2*xdft(2:end-1);    % account for magnitude of "de-folding"

ydft = fft(y);
ydft = ydft(1:length(y)/2);
ydft(2:end-1) = 2*ydft(2:end-1);

S = (1/T)*conj(xdft).*ydft;

end

