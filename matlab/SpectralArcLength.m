function S = SpectralArcLength( speed, Ts, parameters )
% SPECTRALARCLENGTH computes the smoothness of the give movement speed 
% profile using the spectral arc length method.
% The this function takes three inputs and provides one output.
% Inputs: { speed, Ts*, parameters* } ('*' optional parameters)
%         SPEED: Speed it the speed profile of the movement. This is a 1xN
%         row vecotr. N is the total number of points in the speed profile.
%         The function assumes that the movement speed profile is already
%         filtered and segemented.
%
%         TS*: Sampling time in seconds. (DEFAULT VALUE = 0.01sec) NOTE: IF
%         YOUR DATA WAS NOT SAMPLED AT 100HZ, YOU MUST ENTER THE
%         APPROPRIATE SAMPLING TIME FOR ACCURATE RESULTS.
%
%         PARAMETERS*: This contains the parameters to be used spectral arc
%         lenght computation. This is a 1x2 column vector. This input
%         argument is option
%           - PARAMETER(1): The amplitude threshold to be used to choose
%           the cut-off frequency. The default value is chosen to be 0.05.
%           - PARAMETER(2): Maximum cut-off frequency for the spectral arc 
%           length calcualtion. (DEFAULT VALUE = 10HZ) NOTE: 20Hz IS USED
%           TO REPRESENT THE MAXIMUM FREQUENCY COMPONENT OF A MOVEMENT.
%           THIS WAS CHOSEN TO COVER BOTH NORMAL AND ABNORMAL MOTOR
%           BEAHVIOUR. YOU CAN USE A VALUE LOWER THAN 20Hz IF YOU ARE AWARE
%           OF THE MAXIMUM FREQUENCY COMPONENT IN THE MOVEMENT OF INTEREST.
%           - PARAMETER(3): Zero padding index. This parameter controls the
%           resolution of the movement spectrum calculated from the speed
%           profile. (DEFAULT VALUE = 4). NOTE: IT IS NOT ADVISABLE TO USE
%           VALUES LOWER THAN 4.
%
% Outputs: { S }
%          S: This is smoothness of the given movement.
%
% For any queries about the method or the code, or if you come across any 
% bugs in the code, feel free to contact me at siva82kb@gmail.com
% Sivakumar Balasubramanian. July 02, 2014. 

% Check input arguments.
if nargin == 0
    disp('Error! Input at least the movement speed profile for the smoothness calculation.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
elseif nargin == 1
    % Default sampling time.
    Ts = 1/100; % 10ms.
elseif nargin == 2
    % Default parameters are use for the spectral arc length caclulations.
    parameters = [0.05, 10, 4];    
end;

% Check if the input argument are of the appropriate dimensions.
% Speed profile.
sz = size(speed);
if (sz(1) == 1) && (sz(2) > 1)
    disp('Error! speed must be a row vector.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;
% Sampling time.
if ~isscalar(Ts)
    disp('Error! Ts must be a scalar.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;
% Parameters.
if length(parameters) ~= 3 
    disp('Error! parameter is a vector with two elements.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;

% Calculate the spectrum of the speed profile.
N = length(speed);
Nfft = 2^(ceil(log2(N))+parameters(3));
speedSpectrum = abs(fft( speed, Nfft ));

% Normalize spectrum with respect to the DC component.
freq = 0:(1/Ts)*(1/Nfft):(1/Ts)*((Nfft-1)/Nfft);
speedSpectrum = speedSpectrum'/max(speedSpectrum);

% Choose the spectrum that is always above the amplitude threshold, and
% within the cut-off frequency.
inxFc = find((freq(1:end) <= parameters(2)) & ...
             (speedSpectrum(1:end) >= parameters(1)), 1, 'last');

% Calculate the spectral arc length.
% 1. select the spectrum of interest.
speedSpectrum = speedSpectrum(1:inxFc);
% 2. Calculate the incremental arc lengths.
dArcLengths = sqrt((1/(inxFc-1))^2 + (diff(speedSpectrum)).^2);
% 4. Compute movement smoothness.
S = -sum(dArcLengths);

return;
