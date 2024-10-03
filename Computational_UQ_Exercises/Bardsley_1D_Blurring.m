%  
%  1d image deblurring with zero boundary conditions.
%
%  Tikhonov regularization is implemented with the regularization parameter 
%  alpha chosen using one of UPRE, GCV, DP, or L-curve. 
%
%  This m-file was used to generate Figures in Chapter 2 of the book.
%
%  written by John Bardsley 2016.
%
%  First, create a Toeplitz matrix A using zero BCs and a Gaussian kernel.
clear all, close all
n      = 256;  % No. of grid points
h      = 1/n;
t      = [h/2 :h: 1-h/2]';
sig    = .03; % kernel width;
%Shifting Gausian Curve to ensure sampled points extracted are w/respect to
%symmetric curve
kernel = (1/sqrt(2*pi)/sig) * exp(-(t-0.5).^2/2/sig^2);


ahat = fft(fftshift(kernel*h));


% Next, set up true solution x_true and data b = A*x_true + error.
x_true  = .75*(.1<t&t<.25) + .25*(.3<t&t<.32) + (.5<t&t<1).*sin(2*pi*t).^4;
x_true  = x_true/norm(x_true);
Ax1      = real(ifft(ahat.*fft(x_true)));
err_lev = 2; % percent error in data.
sigma = err_lev/100 * norm(Ax1)/ sqrt(n);
eta     =  sigma * randn(n,1);
b       = Ax1 + eta;

figure(1), 
  plot(t,x_true,t,b,'o')

alpha = 0.01;
%
% % Now compute the regularized solution and plot it.
xfilt = ifft(diag(conj(fft(fftshift(kernel)*h))./(abs(conj(fft(fftshift(kernel)*h)).^2) + alpha*1)))*fft(b);
rel_error = norm(xfilt-x_true)/norm(x_true);
figure(2)
  plot(t,x_true,'b-',t,xfilt,'k-')


