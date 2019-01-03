%
% Programmed by Chanwoo Kim for the IEEETran Speech, Audio, and Langauge Processing and ICASSP 2012
%
% (chanwook@cs.cmu.edu)
%
% Important: The input should be mono and be sampled at "16 kHz".
%
%
% * If you want to use logarithmic nonlinearity instead of the power
% nonlinearity, change bPowerLaw to 0 (lilne 28)
%
% PNCC_IEEETran(OutFile, InFile)
%
% iFiltTyep 2 Gamma
% iFiltType 1 HTK MEL
% iFiltType 0 Snaley MEL
% Default
% 0.5, 0.01, 2, 4
%
function [aadDCT] = PNCC(rawdata, fsamp)

	ad_x = rawdata;
 
	%addpath voicebox/;	% With Spectral Subtraction - default parameters
        %ad_x = specsub(rawdata, fsamp);	

	dLamda_L = 0.999;
	dLamda_S = 0.999;

	dSampRate   = fsamp;
	dLowFreq      = 200;% Changed to 40 from 200 as low freq is 40 in gabor as well
	dHighFreq     = dSampRate / 2;
	dPowerCoeff = 1 / 15;

	iFiltType = 1;
	dFactor = 2.0;

	dGammaThreshold = 0.005;

	iM = 0; % Changed from 2 to 0 as number of frames coming out to be different due to queue
	iN = 4;

	iSMType = 0;
    
	dLamda  = 0.999;
	dLamda2 = 0.5;
	dDelta1 = 1;

	dLamda3 = 0.85;
	dDelta2 = 0.2;
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Flags
       %
	bPreem         = 1;	% pre-emphasis flag
	bSSF             = 1;
	bPowerLaw    = 1;
	bDisplay        = 0;
     

	dFrameLen     = 0.025;  % 25.6 ms window length, which is the default setting in CMU Sphinx
	dFramePeriod = 0.010;   % 10 ms frame period
	iPowerFactor  = 1;

	global	iNumFilts;
	iNumFilts = 40;
	
        if iNumFilts<30
           iFFTSize  = 512;
        else
           iFFTSize  = 1024;
        end

	% For derivatives
	deltawindow = 2;	% to calculate 1st derivative
	accwindow   = 2;        % to calculate 2nd derivative

       %	numcoeffs = 13;		% number of cepstral coefficients to be used
	numcoeffs = 13;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Flags
	%
	%
	% Array Queue Ring-buffer
	%
	global Queue_aad_P;
	global Queue_iHead;
	global Queue_iTail;
	global Queue_iWindow;
	global Queue_iNumElem;

	Queue_iWindow  = 2 * iM + 1;
	Queue_aad_P    = zeros(Queue_iWindow, iNumFilts);
	Queue_iHead    = 0;
	Queue_iTail    = 0;
	Queue_iNumElem = 0;
   
	iFL        = floor(dFrameLen    * dSampRate);
	iFP        = floor(dFramePeriod * dSampRate);
	iNumFrames = floor((length(ad_x) - iFL) / iFP) + 1;
	iSpeechLen = length(ad_x);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Pre-emphasis using H(z) = 1 - 0.97 z ^ -1
	%
	if (bPreem == 1)
	    ad_x = filter([1 -0.97], 1, ad_x);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Obtaning the gammatone coefficient. 
	%
	% Based on M. Snelly's auditory toolbox. 
	% In actual C-implementation, we just use a table
	%
	bGamma = 1;
    
	[wts,binfrqs]  = fft2melmx(iFFTSize, dSampRate, iNumFilts, 1, dLowFreq, dHighFreq, iFiltType);
	wts = wts';
	wts(size(wts, 1) / 2 + 1 : size(wts, 1), : ) = [];
	aad_H = wts;
    
	i_FI     = 0;
	i_FI_Out = 0;

	if bSSF == 1
	    adSumPower = zeros(1, iNumFrames - 2 * iM);
	else
	    adSumPower = zeros(1, iNumFrames);
	end
	 
	%dLamda_L   = 0.998;
	aad_P      = zeros(iNumFrames,      iNumFilts);
	aad_P_Out  = zeros(iNumFrames - 2 * iM,      iNumFilts);
	ad_Q       = zeros(1,               iNumFilts);
	ad_Q_Out   = zeros(1,               iNumFilts);
	ad_QMVAvg  = zeros(1,               iNumFilts);
	ad_w       = zeros(1,               iNumFilts);
	ad_w_sm    = zeros(1,               iNumFilts);
	ad_QMVAvg_LA = zeros(1,               iNumFilts);

	MEAN_POWER = 1e10;

	dMean  = 5.8471e+08;
	dPeak = 2.7873e+09 / 15.6250;
	% (1.7839e8, 2.0517e8, 2.4120e8, 2.9715e8, 3.9795e8) 95, 96, 97, 98, 99
	% percentile from WSJ-si84
	            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dPeakVal = 4e+07;% % 4.0638e+07  --> Mean from WSJ0-si84  (Important!!!)
	                %%%%%%%%%%%
	dMean = dPeakVal;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Obtaining the short-time Power P(i, j)
	%
	for m = 0 : iFP : iSpeechLen  - iFL 
		ad_x_st                = ad_x(m + 1 : m + iFL) .* hamming(iFL);
		adSpec                 = fft(ad_x_st, iFFTSize);
		ad_X                   = abs(adSpec(1: iFFTSize / 2));
		aadX(:, i_FI + 1)      = ad_X; 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% Calculating the Power P(i, j)
		%
		for j = 1 : iNumFilts
		        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		        %
		        % Squared integration
		        %
		        
		        if iFiltType == 2
		            aad_P(i_FI + 1, j)  = sum((ad_X .* aad_H(:, j)) .^ 2);
		        else
		            aad_P(i_FI + 1, j)  = sum((ad_X .^ 2 .* aad_H(:, j)));
		        end
		        
		end
        
	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        %
	        % Calculating the Power P(i, j)
	        %
	        
	        dSumPower = sum(aad_P(i_FI + 1, : ));
	             
	        if bSSF == 1
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%
			% Ring buffer (using a Queue)
			%
			if (i_FI >= 2 * iM + 1)
				Queue_poll();
				end
				Queue_offer(aad_P(i_FI + 1, :));

				ad_Q = Queue_avg();

			if (i_FI == 2 * iM)
				ad_QMVAvg     = ad_Q.^ (1 / 15);
				ad_PBias  =  (ad_Q) * 0.9;
			end
          
			if (i_FI >= 2 * iM)  
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%
				% Bias Update
				%
				for i = 1 : iNumFilts,
				    if (ad_Q(i) > ad_PBias(i))
				       ad_PBias(i) = dLamda * ad_PBias(i)  + (1 - dLamda) * ad_Q(i);
				    else
				       ad_PBias(i) = dLamda2 * ad_PBias(i) + (1 - dLamda2) * ad_Q(i);
				    end
				end
			            
			       for i = 1 : iNumFilts,
				ad_Q_Out(i) =   max(ad_Q(i) - ad_PBias(i), 0) ;

				if (i_FI == 2 * iM)
				    ad_QMVAvg2(i)  =  0.9 * ad_Q_Out(i);
				    ad_QMVAvg3(i)  =  ad_Q_Out(i);
				    ad_QMVPeak(i)  =  ad_Q_Out(i);
				end

				if (ad_Q_Out(i) > ad_QMVAvg2(i))
				     ad_QMVAvg2(i) = dLamda * ad_QMVAvg2(i)  + (1 -  dLamda)  *  ad_Q_Out(i);
				else
				     ad_QMVAvg2(i) = dLamda2 * ad_QMVAvg2(i) + (1 -  dLamda2) *  ad_Q_Out(i);
				end

				dOrg =  ad_Q_Out(i);

				ad_QMVAvg3(i) = dLamda3 * ad_QMVAvg3(i);
			          
				if (ad_Q(i) <  dFactor * ad_PBias(i))
				    ad_Q_Out(i) = ad_QMVAvg2(i);
				else
				     if (ad_Q_Out(i) <= dDelta1 *  ad_QMVAvg3(i))
				        ad_Q_Out(i) = dDelta2 * ad_QMVAvg3(i);
				     end
				end
				ad_QMVAvg3(i) = max(ad_QMVAvg3(i),   dOrg);

				ad_Q_Out(i) =  max(ad_Q_Out(i), ad_QMVAvg2(i));
			end
			      ad_w      =   ad_Q_Out ./ max(ad_Q, eps);

			for i = 1 : iNumFilts,
				 if iSMType == 0
			                ad_w_sm(i) = mean(ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts)));
			        elseif iSMType == 1
			                ad_w_sm(i) = exp(mean(log(ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts)))));
			        elseif iSMType == 2
			        	 ad_w_sm(i) = mean((ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts))).^(1/15))^15;
			        elseif iSMType == 3
			 		ad_w_sm(i) = (mean(  (ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts))).^15 )) ^ (1 / 15); 
			        end
			end        
			        
		        aad_P_Out(i_FI_Out + 1, :) = ad_w_sm .* aad_P(i_FI - iM + 1, :);
		        adSumPower(i_FI_Out + 1)   = sum(aad_P_Out(i_FI_Out + 1, :));

		        if  adSumPower(i_FI_Out + 1) > dMean
		             dMean = dLamda_S * dMean + (1 - dLamda_S) * adSumPower(i_FI_Out + 1);
		        else
		             dMean = dLamda_L * dMean + (1 - dLamda_L) * adSumPower(i_FI_Out + 1);
		        end
		        
		        aad_P_Out(i_FI_Out + 1, :) = aad_P_Out(i_FI_Out + 1, :) / (dMean)  * MEAN_POWER;
		        i_FI_Out = i_FI_Out + 1;
		        
		end
           
	else % if not SSF
		adSumPower(i_FI + 1)   = sum(aad_P(i_FI + 1, :));
             
		if  adSumPower(i_FI_Out + 1) > dMean
		     dMean = dLamda_S * dMean + (1 - dLamda_S) * adSumPower(i_FI_Out + 1);
		else
		     dMean = dLamda_L * dMean + (1 - dLamda_L) * adSumPower(i_FI_Out + 1);
		end

		aad_P_Out(i_FI + 1, :) = aad_P(i_FI + 1, :) / (dMean)  * MEAN_POWER;
		end
		i_FI = i_FI + 1;
	end

	
	%adSorted  = sort(adSumPower);
	%dMaxPower = adSorted(round(0.98 * length(adSumPower)));
	%aad_P_Out = aad_P_Out / dMaxPower * 1e10;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Apply the nonlinearity
	%
	%dPowerCoeff
	if bPowerLaw == 1
	    aadSpec = aad_P_Out .^ dPowerCoeff;
	else
	    aadSpec = log(aad_P_Out + eps);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% DCT
	%
	aadDCT                  = dct(aadSpec')';
	
	aadDCT(:, numcoeffs+1:iNumFilts) = [];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% MVN
	%
	for i = 1 : numcoeffs
	       aadDCT( :, i ) = (aadDCT( : , i ) - mean(aadDCT( : , i)))/std(aadDCT(:,i));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Temporal Derivatives
	% calculate 1st derivative (velocity)
	dt1 = deltacc(aadDCT', deltawindow);

	% calculate 2nd derivative (acceleration)
	dt2 = deltacc(dt1, accwindow);
 	% append dt1 and dt2 to mfcco
	aadDCT = [aadDCT'; dt2];
        % aadDCT = [aadDCT'; dt2];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Display
	%
	if bDisplay == 1
	    figure
	    
	    aadSpec = idct(aadDCT', iNumFilts);
	    imagesc(aadSpec); axis xy;
	end

	aadDCT = aadDCT';
	
	%{
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Writing the feature in Sphinx format
	%
	[iM, iN] = size(aadDCT);
	iNumData = iM * iN;
	fid = fopen(szOutFeatFileName, 'wb');
	fwrite(fid, iNumData, 'int32');
	iCount = fwrite(fid, aadDCT(:), 'float32');
	fclose(fid);
	%}
    
end


function dt = deltacc(input, winlen)
% calculates derivatives of a matrix, whose columns are feature vectors

tmp = 0;
for cnt = 1 : winlen
    tmp = tmp + cnt*cnt;
end
nrm = 1 / (2*tmp);

dt   = zeros(size(input));
rows = size(input,1);
cols = size(input,2);
for col = 1 : cols
    for cnt = 1 : winlen
        inx1 = col - cnt; inx2 = col + cnt;
        if inx1 < 1;     inx1 = 1;     end
        if inx2 > cols;  inx2 = cols;  end
        dt(:, col) = dt(:, col) + (input(:, inx2) - input(:, inx1)) * cnt;
    end
end
dt = dt * nrm;
end

function [] = Queue_offer(ad_x)
    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;
    
    Queue_aad_P(Queue_iTail + 1, :) = ad_x;
    
    Queue_iTail    = mod(Queue_iTail + 1, Queue_iWindow);
    Queue_iNumElem = Queue_iNumElem + 1;
    
    if Queue_iNumElem > Queue_iWindow
       error ('Queue overflow'); 
    end
    
  
end


function [ad_x] = Queue_poll()
    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;
    
   
    
    if Queue_iNumElem <= 0
       error ('No elements'); 
    end
    
    
    ad_x =  Queue_aad_P(Queue_iHead + 1, :);
    
    Queue_iHead    = mod(Queue_iHead + 1, Queue_iWindow);
    Queue_iNumElem = Queue_iNumElem - 1;
 
end


function[adMean] = Queue_avg()

    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;
    global iNumFilts;
    
    adMean = zeros(1, iNumFilts);  % Changed from 40 (number of filter banks)

    
    iPos = Queue_iHead;
    
    
    for i = 1 : Queue_iNumElem
        adMean = adMean + Queue_aad_P(iPos + 1 ,: );
        iPos   = mod(iPos + 1, Queue_iWindow);
    end
    
    adMean = adMean / Queue_iNumElem;

end

function [wts,binfrqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
% wts = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
%      Generate a matrix of weights to combine FFT bins into Mel
%      bins.  nfft defines the source FFT size at sampling rate sr.
%      Optional nfilts specifies the number of output bands required 
%      (else one per bark), and width is the constant width of each 
%      band relative to standard Mel (default 1).
%      While wts has nfft columns, the second half are all zero. 
%      Hence, Mel spectrum is fft2melmx(nfft,sr)*abs(fft(xincols,nfft));
%      minfrq is the frequency (in Hz) of the lowest band edge;
%      default is 0, but 133.33 is a common standard (to skip LF).
%      maxfrq is frequency in Hz of upper edge; default sr/2.
%      You can exactly duplicate the mel matrix in Slaney's mfcc.m
%      as fft2melmx(512, 8000, 40, 1, 133.33, 6855.5, 0);
%      htkmel=1 means use HTK's version of the mel curve, not Slaney's.
%      constamp=1 means make integration windows peak at 1, not sum to 1.
% 2004-09-05  dpwe@ee.columbia.edu  based on fft2barkmx

if nargin < 2;     sr = 8000;      end
if nargin < 3;     nfilts = 40;    end
if nargin < 4;     width = 1.0;    end
if nargin < 5;     minfrq = 0;     end  % default bottom edge at 0
if nargin < 6;     maxfrq = sr/2;  end  % default top edge at nyquist
if nargin < 7;     htkmel = 0;     end
if nargin < 8;     constamp = 0;   end


wts = zeros(nfilts, nfft);

% Center freqs of each FFT bin
fftfrqs = [0:nfft-1]/nfft*sr;

% 'Center freqs' of mel bands - uniformly spaced between limits
minmel = hz2mel(minfrq, htkmel);
maxmel = hz2mel(maxfrq, htkmel);
binfrqs = mel2hz(minmel+[0:(nfilts+1)]/(nfilts+1)*(maxmel-minmel), htkmel);

binbin = round(binfrqs/sr*(nfft-1));

for i = 1:nfilts
%  fs = mel2hz(i + [-1 0 1], htkmel);
  fs = binfrqs(i+[0 1 2]);
  % scale by width
  fs = fs(2)+width*(fs - fs(2));
  % lower and upper slopes for all bins
  loslope = (fftfrqs - fs(1))/(fs(2) - fs(1));
  hislope = (fs(3) - fftfrqs)/(fs(3) - fs(2));
  % .. then intersect them with each other and zero
%  wts(i,:) = 2/(fs(3)-fs(1))*max(0,min(loslope, hislope));
  wts(i,:) = max(0,min(loslope, hislope));

  % actual algo and weighting in feacalc (more or less)
%  wts(i,:) = 0;
%  ww = binbin(i+2)-binbin(i);
%  usl = binbin(i+1)-binbin(i);
%  wts(i,1+binbin(i)+[1:usl]) = 2/ww * [1:usl]/usl;
%  dsl = binbin(i+2)-binbin(i+1);
%  wts(i,1+binbin(i+1)+[1:(dsl-1)]) = 2/ww * [(dsl-1):-1:1]/dsl;
% need to disable weighting below if you use this one

end

if (constamp == 0)
  % Slaney-style mel is scaled to be approx constant E per channel
  wts = diag(2./(binfrqs(2+[1:nfilts])-binfrqs(1:nfilts)))*wts;
end

% Make sure 2nd half of FFT is zero
wts(:,(nfft/2+1):nfft) = 0;
% seems like a good idea to avoid aliasing

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = mel2hz(z, htk)
%   f = mel2hz(z, htk)
%   Convert 'mel scale' frequencies into Hz
%   Optional htk = 1 means use the HTK formula
%   else use the formula from Slaney's mfcc.m
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
  htk = 0;
end

if htk == 1
  f = 700*(10.^(z/2595)-1);
else
  
  f_0 = 0; % 133.33333;
  f_sp = 200/3; % 66.66667;
  brkfrq = 1000;
  brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
  logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

  linpts = (z < brkpt);

  f = 0*z;

  % fill in parts separately
  f(linpts) = f_0 + f_sp*z(linpts);
  f(~linpts) = brkfrq*exp(log(logstep)*(z(~linpts)-brkpt));

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = hz2mel(f,htk)
%  z = hz2mel(f,htk)
%  Convert frequencies f (in Hz) to mel 'scale'.
%  Optional htk = 1 uses the mel axis defined in the HTKBook
%  otherwise use Slaney's formula
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
  htk = 0;
end

if htk == 1
  z = 2595 * log10(1+f/700);
else
  % Mel fn to match Slaney's Auditory Toolbox mfcc.m

  f_0 = 0; % 133.33333;
  f_sp = 200/3; % 66.66667;
  brkfrq = 1000;
  brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
  logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

  linpts = (f < brkfrq);

  z = 0*f;

  % fill in parts separately
  z(linpts) = (f(linpts) - f_0)/f_sp;
  z(~linpts) = brkpt+(log(f(~linpts)/brkfrq))./log(logstep);

end
end
