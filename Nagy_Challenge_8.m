F = imread('butterflies.tif');
G = imread('iogray.tif');
Gd = double(G);
%figure();
%imagesc(Gd), colormap("gray");
%figure();
%imshow(G);

%Blurring with zero boundary conditions
zeromat = zeros(512,512);
Gzero = [zeromat zeromat zeromat; zeromat Gd zeromat; zeromat zeromat zeromat];

P = psfGauss(32,1);
Bzero = conv2(Gzero, P, 'same');
Bzeroext = Bzero(513:1024, 513:1024);

figure, subplot(1,2,1)
imagesc(Bzeroext), colormap('gray');


%Blurring with periodic boundary conditions
Gperiod = [Gd Gd Gd; Gd Gd Gd; Gd Gd Gd];
Bperiod = conv2(Gperiod, P, "same");
Bperiodext = Bperiod(513:1024, 513:1024);


subplot(1,2,2)
imagesc(Bperiodext), colormap('gray');







