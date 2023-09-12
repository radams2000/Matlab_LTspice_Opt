function [target,errWeights] = setTarget(freqx,matchMode)
%setTarget. User defined target response. The user-defined response must
%be calculated at the same frequencies used in the LTSpice sim
% freqx = freqs returned by initial LTSpice simulation
% matchMode = 1,2 or 3 (ampl only, phase only, or both)

global example; % set in simControl.m


%% example1 circuit 1, match a single op-amp design to a 3rd-order lowpass.

if example==1

    % **************** design a 5-th-order chebychev in S **************
    Wp = 1; % we will scale the axtual cutoff later
    [z,p,k] = cheby1(3,1,Wp,'s'); % 3th-order 1dB chebychev lowpass, w=1rad/sec
    a1 = poly(p); % coefficients of 1/(p(1)*s^2 + p(2)*s + p(3))
    b1 = -1; % use neg sign because it's an inverting filter; only important in case where you are matching phase
    fcutoff = 20e3; % cutoff frequency in Hz
    Hampl = abs(freqs(b1,a1,freqx./(fcutoff)));
    Hphase = unwrap(angle(freqs(b1,a1,freqx./(fcutoff))));
   

    if matchMode==1 % match ampl only
        target_ampl = Hampl/Hampl(1); % DC gain = 1
        target_phase = ones(size(freqx)); % not used
        target = target_ampl; % this is what is returned by function
        errWeights_ampl = ones(size(freqx)); 
        errWeights_phase = ones(size(freqx)); % not used
        errWeights = errWeights_ampl; % errWeights is returned by function
    end
    if matchMode==2 % match phase only
        target_ampl = ones(size(freqx)); % not used
        target_phase = Hphase; % used
        target = target_phase; % returned by function
        errWeights_ampl = ones(size(freqx)); % not used
        errWeights_phase = ones(size(freqx));
        errWeights = errWeights_phase;% errWeights is returned by function
    end
    if matchMode==3 % match both, concatenate ampl and phase for double-length vector
        target_ampl = Hampl;
        target_phase = Hphase;
        target = [target_ampl target_phase];% concanenate ampl + phase, returned by function
        errWeights_ampl = ones(size(freqx));
        errWeights_phase = ones(size(freqx));
        errWeights=[errWeights_ampl errWeights_phase];% concatenate ampl+phase, errWeights is returned by function
    end
     
    % modify errWeights if you want better convergence in some spectral
    % areas at the expense of others
    % figure;
    % if matchMode==1
    %     semilogx(freqx,20*log10(target),'g');
    %     title('chebychev target ampl dB');
    % end
    % if matchMode==2
    %     semilogx(freqx,target,'g');
    %      title('chebychev target phase ');
    % end
    % 
    % if matchMode==3
    %     subplot(2,1,1),semilogx(freqx,20*log10(target(1:end/2)),'g');
    %     title('chebychev target ampl dB');
    %     subplot(2,1,2),semilogx(freqx,target(end/2+1:end),'g');
    %     title('chebychev target phase radians ');
    % end


end



%% example 2. Complicated single-op-amp filter with an LC stopband notch, 200KHz passband.
%% example of classic passband/stopband/ transition-band specification style
if example==2
    target = ones(size(freqx));
    errWeights = ones(size(freqx));
    f1i = find(freqx > 210e3,1); % edge of passband
    f1ai = find(freqx > 150e3,1); % start of region where passband weight should be increased, to counter passband ripples at edge
    f2i = find(freqx > 290e3,1); % start of stopband
    target(1:f1i)=1.0;
    target(f2i:end) = 5e-7; % not 0, just for plotting reasons (log10 will blow up ...)
    target(f1i+1:f2i-1) = 5e-7; % don't-care band, weights in this region are set to 0
    errWeights(1:f1i) = 1; % passband error weights
    errWeights(f2i:end) = 10; % stopband error weights
    errWeights(f1i+1:f2i-1) = 0; % don't-care band error weights
    % compute a linear increase in weights near the passband edge to squash
    % down the ripples
    deltai = f1i-f1ai;
    rng = f1ai:f1i;
    errWeights(rng) = 1 + 4*(rng-f1ai)/deltai;
    
    % figure;
    % subplot(2,1,1), semilogx(freqx,target);
    % title('target response');
    % ylim([min(target)-0.5 max(target)+0.5]);
    % subplot(2,1,2), semilogx(freqx,errWeights);
    % ylim([min(errWeights)-0.5 max(errWeights)+0.5]);
    % title('errWeights');
end



%% Your circuit here, must set target and errWeights (vectors with size of freqx)
%% In the case of matchMode==3 (match both ampl and phase), target and errWeights must have 2x the length of freqx
% **************************************
% ***************************************
%
if example==3
    X = linspace(0, 80000, 21);  % every 4KHz from datasheet plot
    VdB = [0, 0.5, 1, 1.5, 2, 3, 4.8, 7, 11, 17, 22, 13, 8, 5, 6, 4.5, 3, 2.5, 3, 3, 4.5];  % from datasheet plot, freq resp in dB
    mic_dB = interp1(X,VdB,freqx);
    mic_lin = 10.^(mic_dB / 20);
    target = 30.0 ./ mic_lin;
    errWeights = sqrt(mic_lin);  % greatest weights where the target is smallest

end