G = imread('iogray.tif');
Gd = double(G);
X = reshape(Gd,[],1);


[Psmall,center] = psfGauss(32,6);
Pbig = padPSF(Psmall, size(Gd));
center = [17,17];

S = fft2(circshift(Pbig,1-center));

j = fft2(X);

B = real(ifft2(S.*fft2(Gd)));
B = B + 1e-2*randn(size(B));

Xrecon = real(ifft2(fft2(B)./S));

%figure()
subplot(1,3,1)
imagesc(Xtrue),colormap('gray')
subplot(1,3,2)
imagesc(B),colormap('gray')
subplot(1,3,3)
imagesc(Xrecon),colormap('gray')