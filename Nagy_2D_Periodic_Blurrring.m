G = imread('iogray.tif');
Gd = double(G);
X = reshape(Gd,[],1);


[Psmall,center] = psfGauss(32,6);
Pbig = padPSF(Psmall, size(Gd));
center = [17,17];

S = fft2(circshift(Pbig,1-center));

j = fft2(X);

B = real(ifft2(S.*fft2(Gd)));

Xtrue = real(ifft2(fft2(B)./S));

%figure()
imagesc(Xtrue),colormap('gray')