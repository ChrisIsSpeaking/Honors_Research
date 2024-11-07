function noise = whiteNoise(b, noiseLevel, rngSeed)
%
%      N = whiteNoise(b, level, seed);
%
%  This function generates Gaussian white noise for the 
%  data b. 
%
%  Input:  b - array containing data
%      level - scalar in [0, 1] specifiying level (percentage) of
%              noise.  For example, level = 0.01 implies
%                  norm(N)/norm(b) = 0.01, or 1% noise
%              Default is level = 0.01.
%       seed - Used to set the random number generator.
%              Default is seed = 0.
%
%  Output: N - array same dimension as b, containing pseudo-random
%              values drawn from a normal distribution with mean zero
%              and standard deviation one, and scaled as described above.
%

% Check inputs and set default values.
if nargin < 2 
  noiseLevel = 0.01;
  if nargin > 2
    rng(rngSeed);
  end
end

% Generate noise.
noise = randn(size(b));
noise = noise / norm(noise(:));
noise = noiseLevel*norm(b(:))*noise;
end
