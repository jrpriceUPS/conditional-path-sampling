function ACF = autocorrelation(data)
%Computes the autocorrelation as a function of lag times for a given matrix
%of data
%
%
%%%%%%%
%INPUT%
%%%%%%%
%
%data      =  a matrix of the data - we compute the autocorrelation of each
%             row in the matrix

%compute the number of samples
[~,n] = size(data);

%save the number of lag times we will compute
nLags = floor(n/2-1);

%normalize the data
data = (data - repmat(mean(data,2),1,n))./(repmat(std(data,0,2),1,n));

%compute the power spectrum
F    =  fft(data,n,2);
F    =  F .* conj(F);

%take the inverse Fourier transform of the power spectrum (this is the
%autocorrelation by the Wiener-Khinchin Theorem)
ACF  =  ifft(F,n,2);

% Retain non-negative lags
ACF  =  ACF(:,1:(nLags + 1));

%Normalize
ACF  =  ACF ./ repmat(ACF(:,1),1,nLags+1);

%Eliminate spurious imaginary parts
ACF  =  real(ACF);