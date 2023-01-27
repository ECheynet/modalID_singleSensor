function [fn,zeta] = modalID_singleSensor(acc,Sacc,f,Nmodes,fs,varargin)
% [fn,zeta] = modalID_singleSensor(ax,Sax,f,Nmodes,fs,varargin) estimates
% the eigen-frequencies and damping ratios from acceleration records
% collected by a single accelerometer mounted on a line-like structure
% (bridge or tower, for example). Note: it is not possible to estimate the mode
% shapes with a single sensor.
% 
% Inputs
%  - acc: [1 x N] double: accelration record at a given location
%  - Sacc: [1 x Nfreq] double: PSD of the acceleration record "acc"
%  - f: [1 x Nfreq] double: frequency associated wih "Sacc"
%  - Nmodes: [1 x 1] scalar: Number of modes to extract
%  - fs: [1 x 1] scalar: sampling frequency
%   varargin:
%       - 'PickingMethod':'auto' or 'manual'
%       - 'fnMin': empty [] or a [1 x Nmodes] vector for the lower cut-off
%       frequency for each mode.
%       - 'fnMax': empty [] or a [1 x Nmodes] vector for the higher cut-off
%       frequency for each mode
%       - 'Nperiod':[1 x 1] scalar: Number of period to include in the IRF
%       for each mode (30 by default)
%       - 'plotOpt': 1 or 0. If "1" is chosen, a figure is plotted to summarize
%       the fitting and filtering steps.
%   
% Outputs
%  - fn: [1 x Nmodes] double: vector of identified eigenfrequencies
%  - zeta: [1 x Nmodes] double: vector of identified damping ratios
% 
% 
% Author: E Cheynet - UiB - Last modified: 01.02.2021
% 

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('PickingMethod','auto'); % automated or manual peak picking ('auto or manual')
p.addOptional('fnMin',[]);
p.addOptional('fnMax',[]); 
p.addOptional('Nperiod',30); 
p.addOptional('plotOpt',0); 
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
fnMin = p.Results.fnMin ;
fnMax = p.Results.fnMax ;
PickingMethod = p.Results.PickingMethod;
Nperiod = p.Results.Nperiod;
plotOpt = p.Results.plotOpt;

if ~strcmpi(PickingMethod,'auto') && ~strcmpi(PickingMethod,'manual')
    error(' ''PickingMethod'' is not recognized. It must be a either ''auto'' or ''manual'' ')
end

dt = 1/fs;
zeta = nan(1,Nmodes);
%% Identify the eigen frequencies by peak-picking
if strcmpi(PickingMethod,'auto')
    [indP,~] = pickpeaks(Sacc,Nmodes,0);
elseif strcmpi(PickingMethod,'manual'),
    indP = manualPickPeaking(f,Sacc,Nmodes);
else
    error('pick-peaking method unknown. Please choose between ''auto'' and ''manual'' ');
end
fn = f(indP);

% sort the eigen frequencies    
[fn,~] = sort(fn);
if isempty(fnMax),    fnMax = fn.*1.05;end % limits of the eigenfrequencies for the filtering
if isempty(fnMin),    fnMin = fn.*0.95;end
%% Band pass filtering for each modes + compute the IRF + fit the IRF
if plotOpt==1
    tiledlayout(Nmodes,2,"Padding","compact","TileSpacing","none")
end

% For each mode selected, a band pass filter is applied + the IRF is
% computed and fitted by a modified exponential decay. 
for pp=1:Nmodes
    
    % Get the time-lag Ts for the autocorrelation function
    Ts = round(Nperiod./fn(pp));
    
    % Design the band-pass filter around the modes
    h1=fdesign.bandpass('N,F3dB1,F3dB2',8,fnMin(pp),fnMax(pp),fs);
    d1 = design(h1,'butter');
    acc_filtered = filtfilt(d1.sosMatrix,d1.ScaleValues, acc);
    
    % Get Impulse-response function 
    [IRF,t] = NExT(acc_filtered,dt,Ts,1);
    
    % Circular eigenfrequency of a lightly damped system
    wd = 2*pi*fn(pp);  
    
    % Fir the IRF
   [zeta(pp),rmse,y,newY,newWd] = expoFit(IRF,t,wd);
   
   % Check if the updated eigenfrquency changes signifciantly from the
   % initial guess by peak-picking
   dummyFn = [fn(pp),newWd./(2*pi)];
   fn(pp) = newWd./(2*pi); 
   if abs(diff(dummyFn)/dummyFn(2))*100>5% if difference larger than 5 %, it is announced
       fprintf(['Improved eigenfrequency estimate obtained \n'])
       fprintf(['Initial estimate: ',num2str(dummyFn(1)),' \n'])
       fprintf(['Final estimate: ',num2str(dummyFn(2)),' \n'])
   end
   
   
   % If plotting option enabled:
    if plotOpt==1
         nexttile
         [S0,f0] = pwelch(detrend(acc),[],[],[],fs);
         [S1,f1] = pwelch(acc_filtered,[],[],[],fs);
         loglog(f0,S0./max(S0),'b'); hold on;loglog(f1,S1./max(S0),'r');hold off; axis tight
         ylabel('PSD');
         if pp==1,legend('before filtering','After filtering','location','best');end
         if pp==Nmodes,xlabel('f (Hz)');end
        
        nexttile
        title(['RMSE = ',num2str(rmse),]);
        plot(t,y,'k','linewidth',1.2);hold on;plot(t,newY,'r','linewidth',0.75); 
        if pp==1,legend('IRF','Best fit','location','best');end
        axis tight;ylim([-1 1]); hold off;
        ylabel('IRF');
        if pp==Nmodes, xlabel('time (s)');end
        set(gcf,'color','w')
    end
end

%% NESTED FUNCTIONS %%
    function [Fp] = manualPickPeaking(f,S,Nmodes)
        % original author:  Mohammad Farshchin
        % FileExchange submission: https://se.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-/content/FDD.m
        %%
        display('Peak selection procedure')
        display('a: Draw rectangles around peaks while holding left click')
        display('b: Press "Space" key to continue the peak selection')
        display('c: Press "any other key" if you have selected a peak by mistake and want to ignore it')
        
        clf;close all;
        figure
        plot(f,mag2db(S))
        grid on
        ylim([min(mag2db(S)),max(mag2db(10*S))])
        xlim([f(2),f(end)])
        hold on
        xlabel('Frequency (Hz)')
        ylabel('1st Singular values of the PSD matrix (db)')
        Fp=[];% Frequencies related to selected peaks
        while numel(Fp)<Nmodes
            myRec=drawrectangle;                                                                          % Draw a rectangle around the peak
            [~,P1]=min(abs(f-myRec.Position(1)));
            [~,P2]=min(abs(f-(myRec.Position(1)+myRec.Position(3))));
            [~,P3]=max(S(P1:P2));
            indPeak=P3+P1-1;                                                                         % Frequency at the selected peak
            scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','b','MarkerFaceColor','b')         % Mark this peak
            pause;
            key=get(gcf,'CurrentKey');
            if strcmp(key,'space'),
                % Press space to continue peak selection
                Fp=[Fp,indPeak];
                scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','g','MarkerFaceColor','g')      % Mark this peak as green
            else
                % Press any other key to ignore this peak
                scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','r','MarkerFaceColor','r')      % Mark this peak as red
            end
        end
        % Number selected peaks, respectively
        Fp=sort(Fp);
        pause(0.01);
    end
    function [peaks,criterion] = pickpeaks(V,select,display)
        % -------------------------------------------------------------
        % Scale-space peak picking
        % ------------------------
        % This function looks for peaks in the data using scale-space theory.
        %
        % input :
        %   * V : data, a vector
        %   * select : either:
        %       - select >1 : the number of peaks to detect
        %       - 0<select<1 : the threshold to apply for finding peaks
        %         the closer to 1, the less peaks, the closer to 0, the more peaks
        %   * display : whether or not to display a figure for the results. 0 by
        %               default
        %   * ... and that's all ! that's the cool thing about the algorithm =)
        %
        % outputs :
        %   * peaks : indices of the peaks
        %   * criterion : the value of the computed criterion. Same
        %                 length as V and giving for each point a high value if
        %                 this point is likely to be a peak
        %
        % The algorithm goes as follows:
        % 1°) set a smoothing horizon, initially 1;
        % 2°) smooth the data using this horizon
        % 3°) find local extrema of this smoothed data
        % 4°) for each of these local extrema, link it to a local extremum found in
        %     the last iteration. (initially just keep them all) and increment the
        %     corresponding criterion using current scale. The
        %     rationale is that a trajectory surviving such smoothing is an important
        %     peak
        % 5°) Iterate to step 2°) using a larger horizon.
        %
        % At the end, we keep the points with the largest criterion as peaks.
        % I don't know if that kind of algorithm has already been published
        % somewhere, I coded it myself and it works pretty nice, so.. enjoy !
        % If you find it useful, please mention it in your studies by referencing
        % the following research report:
        %
        %@techreport{liutkus:hal-01103123,
        %  TITLE = {{Scale-Space Peak Picking}},
        %  AUTHOR = {Liutkus, Antoine},
        %  URL = {https://hal.inria.fr/hal-01103123},
        %  TYPE = {Research Report},
        %  INSTITUTION = {{Inria Nancy - Grand Est (Villers-l{\`e}s-Nancy, France)}},
        %  YEAR = {2015},
        %  MONTH = Jan,
        %  KEYWORDS = { scale-space ; peak detection},
        %  HAL_ID = {hal-01103123},
        %  HAL_VERSION = {v1},
        %}
        %
        %
        % running time should be decent, although intrinsically higher than
        % findpeaks. For vectors of length up to, say, 10 000, it should be nice.
        % Above, it may be worth it though.
        % ---------------------------------------------------------------------
        % Copyright (C) 2015, Inria, Antoine Liutkus
        %
        %Redistribution and use in source and binary forms, with or without
        %modification, are permitted provided that the following conditions are met:
        %    * Redistributions of source code must retain the above copyright
        %      notice, this list of conditions and the following disclaimer.
        %    * Redistributions in binary form must reproduce the above copyright
        %      notice, this list of conditions and the following disclaimer in the
        %       documentation and/or other materials provided with the distribution.
        %     * Neither the name of Inria nor the names of its contributors may
        %       be used to endorse or promote products derived from this software
        %       without specific prior written permission.
        %
        % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        % ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        % WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        % DISCLAIMED. IN NO EVENT SHALL INRIA BE LIABLE FOR ANY
        % DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        % (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        % LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        % ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        % (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        
        %data is a vector
        V = V(:)-min((V(:)));
        
        %input parsin
        if nargin < 3
            display=0;
        end
        if nargin < 2
            select= 0;
        end
        
        n = length(V);
        
        %definition of local variables
        buffer = zeros(n,1);
        criterion = zeros(n,1);
        if select < 1
            minDist = n/20;
        else
            minDist = n/select;
        end
        %horizons = round(linspace(1,ceil(n/20),50));
        horizons = unique(round(logspace(0,2,50)/100*ceil(n/20)));
        
        %horizons=1:2:50;
        Vorig = V;
        
        % all this tempMat stuff is to avoid calling findpeaks which is horribly
        % slow for our purpose
        tempMat = zeros(n,3);
        tempMat(1,1)=inf;
        tempMat(end,3)=inf;
        
        % loop over scales
        for is=1:length(horizons)
            
            %sooth data, using fft-based convolution with a half sinusoid
            horizon = horizons(is);
            if horizon > 1
                w=max(eps,sin(2*pi*(0:(horizon-1))/2/(horizon-1)));
                w=w/sum(w);
                %V=conv(V(:),w(:),'same');
                V = real(ifft(fft(V(:),n+horizon).*fft(w(:),n+horizon)));
                V = V(1+floor(horizon/2):end-ceil(horizon/2));
            end
            
            %find local maxima
            tempMat(2:end,1) = V(1:end-1);
            tempMat(:,2) = V(:);
            tempMat(1:end-1,3) = V(2:end);
            [useless,posMax] =max(tempMat,[],2);
            I = find(posMax==2);
            I = I(:)';
            
            %initialize buffer
            newBuffer = zeros(size(buffer));
            
            if is == 1
                % if first iteration, keep all local maxima
                newBuffer(I) = Vorig(I);
            else
                old = find(buffer);
                old = old(:)';
                if isempty(old)
                    continue;
                end
                
                %Now, for each element of I, find the closest element in
                %old along with its distance. The few nice lines below were
                %written by Roger Stafford in a forum post available here:
                %http://www.mathworks.fr/matlabcentral/newsreader/view_thread/24387
                [c,p] = sort(old);
                [useless,ic] = histc(I,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
                iOld = p(ic);
                d = abs(I-old(iOld));
                
                %done, now select only those that are sufficiently close
                neighbours = iOld(d<minDist);
                
                if ~isempty(neighbours)
                    newBuffer(old(neighbours)) = V(old(neighbours))*is^2;
                end
            end
            %update stuff
            buffer = newBuffer;
            criterion = criterion + newBuffer;
        end
        
        %normalize criterion
        criterion = criterion/max(criterion);
        
        %find peaks based on criterion
        if select<1
            peaks = find(criterion>select);
        else
            %     sorted = find(criterion>1E-3);
            %     [~,order] = sort(criterion(sorted),'descend');
            %     peaks = sorted(order(1:min(length(sorted),select)));
            [useless,order] = sort(criterion,'descend');
            peaks = order(1:select);
        end
        
        if display
            %display
            clf
            plot(Vorig,'LineWidth',2);
            hold on
            plot(criterion*max(Vorig),'r');
            hold on
            plot(peaks,Vorig(peaks),'ro','MarkerSize',10,'LineWidth',2)
            grid on
            title('Scale-space peak detection','FontSize',16);
            legend('data','computed criterion','selected peaks');
        end
    end
    function [IRF,t] = NExT(y,dt,Ts,method)
        %
        % [IRF] = NExT(y,ys,T,dt) implements the Natural Excitation Technique to
        % retrieve the Impulse Response FUnction (IRF) from the cross-correlation
        % of the measured output y.
        %
        % [IRF] = NExT(y,dt,Ts,1) calculate the IRF with cross-correlation
        % calculated by using the inverse fast fourier transform of the
        % cross-spectral power densities  (method = 1).
        %
        % [IRF] = NExT(y,dt,Ts,2) calculate the IRF with cross-correlation
        % calculated by using the unbiased cross-covariance function (method = 2)
        %
        %
        % y: time series of ambient vibrations: vector of size [1xN]
        % dt : Time step
        % method: 1 or 2 for the computation of cross-correlation functions
        % T: Duration of subsegments (T<dt*(numel(y)-1))
        % IRF: impusle response function
        % t: time vector asociated to IRF
        %%
        if nargin<4, method = 2; end % the fastest method is the default method
        if ~ismatrix(y), error('Error: y must be a vector or a matrix'),end
        
        
        [Nyy,N]=size(y);
        if Nyy>N
            y=y';
            [Nyy,N]=size(y);
        end
        
        % get the maximal segment length fixed by T
        M = round(Ts/dt);
        switch method
            case 1
                clear IRF
                for ii=1:Nyy
                    for jj=1:Nyy
                        y1 = fft(y(ii,:));
                        y2 = fft(y(jj,:));
                        h0 = ifft(y1.*conj(y2));
                        IRF(ii,jj,:) = h0(1:M);
                    end
                end
                % get time vector t associated to the IRF
                t = linspace(0,dt.*(size(IRF,3)-1),size(IRF,3));
                if Nyy==1
                    IRF = squeeze(IRF)'; % if Nyy=1
                end
            case 2
                IRF = zeros(Nyy,Nyy,M+1);
                for ii=1:Nyy
                    for jj=1:Nyy
                        [dummy,lag]=xcov(y(ii,:),y(jj,:),M,'unbiased');
                        IRF(ii,jj,:) = dummy(end-round(numel(dummy)/2)+1:end);
                    end
                end
                if Nyy==1
                    IRF = squeeze(IRF)'; % if Nyy=1
                end
                % get time vector t associated to the IRF
                t = dt.*lag(end-round(numel(lag)/2)+1:end);
        end
        % normalize the IRF
        if Nyy==1
            IRF = IRF./IRF(1);
        else
        end
        
        
        
    end
    function [zeta,rmse,y,newY,newWd] = expoFit(IRF,t,wd)
        % [zeta] = expoFit(y,t,wn) returns the damping ratio calcualted by fiting
        % an exponential decay to the envelop of the Impulse Response Function.
        %
        % y: envelop of the IRF: vector of size [1 x N]
        % t: time vector [ 1 x N]
        % wd: damped eigen frequencies (rad/Hz) :  [1 x 1]
        % wn: undamped eigen frequencies (rad/Hz) :  [1 x 1]
        % zeta: modal damping ratio:  [1 x 1]

        assert(license('test','optimization_Toolbox')==1,'The function expoFit requires Matlab Statistics Toolbox.')
        % Normalization of y
        y = IRF./IRF(1);
        options=optimset('Display','off');

        % Initialisation
        guess = [wd,1e-2,pi/2];
        % simple exponentiald ecay function
%         myFun = @(a,x) exp(-a(1).*x);
            
%         myFun = @(a,t) exp(-a(1).*a(2)).*cos(a(2).*sqrt(1-a(1).^2).*t+a(3));
        % application of nlinfit function
        myFun = @SDOF_dampedSystem;
        [coeff,mse] = lsqcurvefit(@(a,t) myFun(a,t),guess,t,y./y(1),[0.95*wd,0,0],[1.05*wd,0.99,2*pi],options);
        rmse = sqrt(mse);
        newY = myFun(coeff,t);
        % modal damping ratio:
        zeta = coeff(2);
        % damped eigen frequency
        newWd = coeff(1);
    end

    function [IRF] = SDOF_dampedSystem(coeff,t)
        % we have initially identified the eigenfrequency of a lightly damped system
        wd = coeff(1); % eigenfrequency of a lightly damped system
        zeta0 = coeff(2); % damping
        phi = coeff(3); % phase 
        wn0 = wd./sqrt(1-zeta0.^2);  % eigenfrequency of undamped system
        IRF = exp(-wn0.*zeta0.*t).*cos(wd.*t + phi);
    end
end

