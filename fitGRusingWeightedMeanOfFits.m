function [gest,gstd] = fitGRusingWeightedMeanOfFits(OD,dt,varargin)
% 
% 
% Required Inputs:
%   OD: vector double, optical density as a function of time (or other 
%       linear pop-size proxy)
%   dt: scalar double, timestep
%       This code assumes constant measurement spacing!!
%       This code is written assuming units of hours. You should be able to
%       adjust it to work 
%       
% Outputs:
%   gest: scalar double, estimated growth rate, same units as dt
%   gstd: scalar double, uncertainty on that estimate, same units as dt
%
% Optional Name-Value Inputs:
%   *NOTE: Planned update to the default values in the near future!*
%   'DontSubtractBG': logical. default false.
%                     if false the line 'OD = OD - min(OD);' is executed
%                     if true that line is not, so set to true if you
%                     prefer an alternative means of background subtraction
%   'MakeFigure': logical. default false. set to true if you want a figure
%                 summarizing the fit to be made
%                 Note: overwrites current figure!!
%   'TitleString': string or character array. adds additional text to the
%                  title of the figure.
%   'MinFitLength': scalar double. This code considers all possible time
%                   windows over which a growth rate could be fit, but has
%                   a minimum time window length. The default value is 3.
%                   This is in the same units as dt, assumed to be hours.
%                   This is done so that (i) noise over the scale of a few
%                   measurements doesn't unnecessarily increase the
%                   uncertainty on the fit and (ii) a couple measurements
%                   happening to perfectly fit the exponential growth model
%                   don't show up as an infinitely better fit than anything
%                   else. Not recommended to go below 3*dt due to (ii).
%   'MinFitFactor': scalar double. Likewise, after fitting growth rates to
%                   all possible time windows, growth rates fits which
%                   correspond to growth by less than a set factor are
%                   excluded. The default value is a factor of 4. Can be
%                   set to 1 to only exclude fits with negative slope or to
%                   0 to not exclude any fits.
%   'MinUsedOD': scalar double. This code excludes all data before the
%                optical density is permanently at or above some minimum
%                threshold below which all data is assume to be noise. The
%                default value is 5e-3. Can be set to 0, but better to use
%                min(OD(OD > 0)) to at least exclude population sizes of 0.
%   'ExcludedEndFactor': scalar double. This code excludes all data after
%                        the OD first reaches within a set factor of its
%                        maximum value. This is to (i) exclude everything
%                        happening after saturation, including slow
%                        residual growth and (ii) exclude any secondary
%                        growth phases. The default value is a factor of 2.
%                        Increase if data is expected to have diauxic
%                        shifts and you want to fit the growth rate to the
%                        initial exponential phase. Must be set to
%                        something greater than (NOT equal to) 1.
%   'NoGrowthThreshold': scalar double. Considers there to be no growth if
%                        the OD never crosses this threshold. Default value
%                        1e-2. Can be set to 0, but note that there are
%                        other conditions under which the code will
%                        consider there to be no growth, for example if
%                        very slow growth mean none of the fits satisfy the
%                        MinFitFactor or MinUsedOD requirement.
%
% Example usage:
%     dt = 10/60;
%     tmax = 24;
%     g = 0.6;
%     n0 = 1e-4;
%     nmax = 0.5;
%     noise = 5e-3;
%     n_t = min(n0*exp(g*(0:dt:tmax))',nmax);
%     noise_t = noise*randn(size(n_t));
%     OD_t = max(n_t + noise_t,0);
%     [gest,gstd] = fitGRusingWeightedMeanOfFits(OD_t,dt,'MakeFigure');
% 
% Copyright 2023 Blox Bloxham


% Default values 
DontSubtractBG = false;
MakeFigure = false;
TitleString = "";
MinFitLength = 3;
MinFitFactor = 4;
MinUsedOD = 5e-3;
ExcludedEndFactor = 2;
NoGrowthThreshold = 1e-2;

i = 1;
while i <= length(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        % Handle DontSubtractBG name-value pair
        if strcmp(varargin{i},'DontSubtractBG')
            if length(varargin) > i && (islogical(varargin{i+1}) || ...
                    (isnumeric(varargin{i+1}) && (varargin{i+1} == 0 || varargin{i+1} == 1)))
                DontSubtractBG = logical(varargin{i+1});
                i = i + 2;
            elseif length(varargin) == i || ischar(varargin{i+1}) || ...
                    isstring(varargin{i+1})
                DontSubtractBG = true;
                i = i + 1;
            else
                error('Invalid value for DontSubtractBG: must be logical')
            end
            
        % Handle MakeFigure name-value pair
        elseif strcmp(varargin{i},'MakeFigure')
            if length(varargin) > i && (islogical(varargin{i+1}) || ...
                    (isnumeric(varargin{i+1}) && (varargin{i+1} == 0 || varargin{i+1} == 1)))
                MakeFigure = logical(varargin{i+1});
                i = i + 2;
            elseif length(varargin) == i || ischar(varargin{i+1}) || ...
                    isstring(varargin{i+1})
                MakeFigure = true;
                i = i + 1;
            else
                error('Invalid value for MakeFigure: must be logical')
            end
            
        % Handle TitleString name-value pair
        elseif strcmp(varargin{i},'TitleString')
            if length(varargin) > i && (ischar(varargin{i+1}) || isstring(varargin{i+1}))
                TitleString = string(varargin{i+1});
                i = i + 2;
            else
                error('Invalid value for TitleString: must be string or char array')
            end
            
        % Handle MinFitLength name-value pair
        elseif strcmp(varargin{i},'MinFitLength')
            if length(varargin) > i && isnumeric(varargin{i+1})
                MinFitLength = varargin{i+1};
                i = i + 2;
            else
                error('Invalid value for MinFitLength: must be numeric')
            end
            
        % Handle MinFitFactor name-value pair
        elseif strcmp(varargin{i},'MinFitFactor')
            if length(varargin) > i && isnumeric(varargin{i+1})
                MinFitFactor = varargin{i+1};
                i = i + 2;
            else
                error('Invalid value for MinFitFactor: must be numeric')
            end
            
        % Handle MinUsedOD name-value pair
        elseif strcmp(varargin{i},'MinUsedOD')
            if length(varargin) > i && isnumeric(varargin{i+1})
                MinUsedOD = varargin{i+1};
                i = i + 2;
            else
                error('Invalid value for MinUsedOD: must be numeric')
            end
            
        % Handle ExcludedEndFactor name-value pair
        elseif strcmp(varargin{i},'ExcludedEndFactor')
            if length(varargin) > i && isnumeric(varargin{i+1})
                ExcludedEndFactor = varargin{i+1};
                i = i + 2;
            else
                error('Invalid value for ExcludedEndFactor: must be numeric')
            end
            
        % Handle NoGrowthThreshold name-value pair
        elseif strcmp(varargin{i},'NoGrowthThreshold')
            if length(varargin) > i && isnumeric(varargin{i+1})
                NoGrowthThreshold = varargin{i+1};
                i = i + 2;
            else
                error('Invalid value for NoGrowthThreshold: must be numeric')
            end
            
        else
            error('Unrecognized name-value pair name: %s',varargin{i})
        end
        
    else
        error('Non-string input as a name-value pair name.')
    end
end


if ~DontSubtractBG
    OD = OD - min(OD);
end

ts = (0:(length(OD)-1))*dt;


% Convert MinFitLength from hours to timepoints
MinFitLength = max(round(MinFitLength/dt),2);

N = size(OD,1);

B = nan(2,N,N);
SS = nan(N);
R2 = nan(N);

for startInd = (find(OD <= min(MinUsedOD,max(OD)/max(10,2*ExcludedEndFactor*MinFitFactor)),1,'last')+1):(N - 2)
    for endInd = (startInd + MinFitLength):find(...
            medfilt1(OD,round(1/dt+1)) > max(OD)/ExcludedEndFactor,1)%N
        % Calculate linear least-squares best fit using:
        %   B = (X' * W * X)^-1 * X' * W * Y
        %   where X = [ones(:,1),t];
        %   where W = diagonal matrix of OD^2 (propto 1/sigma ~1/(dlogOD/dt)^2)
        %   and where Y = log(OD)
        B(:,startInd,endInd) = ...
            ([ones(endInd-startInd+1,1),((startInd:endInd)'-1)*dt]' * ...   % (X'*
            (eye(endInd-startInd+1) .* OD(startInd:endInd).^2) * ...        %  W *
            [ones(endInd-startInd+1,1),((startInd:endInd)'-1)*dt])^-1 * ... %  X)^-1
            [ones(endInd-startInd+1,1),((startInd:endInd)'-1)*dt]' * ...    % * X'
            (eye(endInd-startInd+1) .* OD(startInd:endInd).^2) * ...        % * W
            log(OD(startInd:endInd));                                       % * Y
        
        % Calculate sum-of-squares and R2 value in (not-log) OD space
        SS(startInd,endInd) = sum( (OD(startInd:endInd) - ...
            exp( B(1,startInd,endInd) + B(2,startInd,endInd)*((startInd:endInd)'-1)*dt)).^2 );
        R2(startInd,endInd) = 1 - SS(startInd,endInd) ./ ...
            sum((OD(startInd:endInd) - mean(OD(startInd:endInd))).^2);
    end
end



ws_all = reshape(1./(1-R2),[],1);
gs_all = reshape(B(2,:,:),[],1);
keep = ~isnan(ws_all) & ~isnan(gs_all) & gs_all > 0 & ...
    reshape(squeeze(B(2,:,:)) .* (ts - ts') > log(MinFitFactor),[],1);
ws = ws_all(keep);
gs = gs_all(keep);

gest = sum(ws.*gs)/sum(ws);
gstd = std(gs,ws);


if sum(keep) == 1
    keep = ~isnan(ws_all) & ~isnan(gs_all) & gs_all > 0 & ...
        reshape(squeeze(B(2,:,:)) .* (ts - ts') > log(MinFitFactor/2),[],1);
    ws = ws_all(keep);
    gs = gs_all(keep);

    gest = sum(ws.*gs)/sum(ws);
    gstd = 2*std(gs,ws);

    
elseif isnan(gest)
    keep = ~isnan(ws_all) & ~isnan(gs_all) & gs_all > 0 & ...
    reshape(squeeze(B(2,:,:)) .* (ts - ts') > log(MinFitFactor/2),[],1);
    ws = ws_all(keep);
    gs = gs_all(keep);

    gest = sum(ws.*gs)/sum(ws)/2;
    gstd = gest;
end

if max(medfilt1(OD,round(1/dt+1))) < NoGrowthThreshold || isnan(gest)
    gest = 0;
    if median(OD(round(N/2):end)) <= 0
        gstd = 0;
    else
        One2N = 1:N;
        B0 = ([ones(sum(OD ~= 0),1),(One2N(OD ~= 0)'-1)*dt]' * ...  % (X'*
            (eye(sum(OD ~= 0)) .* OD(OD ~= 0).^2) * ...             %  W *
            [ones(sum(OD ~= 0),1),(One2N(OD ~= 0)'-1)*dt])^-1 * ... %  X)^-1
            [ones(sum(OD ~= 0),1),(One2N(OD ~= 0)'-1)*dt]' * ...    % * X'
            (eye(sum(OD ~= 0)) .* OD(OD ~= 0).^2) * ...             % * W
            log(OD(OD ~= 0));
        gstd = B0(2);
    end
end


if MakeFigure
    
    subplot(2,2,3)
    semilogx(ws,gs,'.')
    xlabel('Weight: 1/(1-R^2)')
    ylabel('Growth Rate (hr^{-1})')
    title('All Linear Fits')
    set(gca,'FontSize',16)
    drawnow
    yl = ylim;
    ylim(yl)

    kdf_gs = linspace(yl(1),yl(2),1000);
    kdfsig = 0.01;

    kdf_density = sum(ws.*exp(-0.5*( (gs - kdf_gs)/kdfsig ).^2)/sqrt(2*pi*kdfsig),1);

    kdf_localmax = find(islocalmax(kdf_density));

    subplot(2,2,4)
    hold off
    plot(kdf_density,kdf_gs,'LineWidth',2)
    hold on
    for i = 1:length(kdf_localmax)
        plot(kdf_density(kdf_localmax(i)),kdf_gs(kdf_localmax(i)),'k.',...
            'MarkerSize',12)
        text(kdf_density(kdf_localmax(i)),kdf_gs(kdf_localmax(i)),...
            sprintf('  %.3f hr^{-1}',kdf_gs(kdf_localmax(i))),'FontSize',12)
    end
    hold off
    xlabel('Weighted Density')
    ylabel('Growth Rate (hr^{-1})')
    title('Density of Fits')
    set(gca,'FontSize',16)
    ylim(yl)

    hold on
    plot(xlim,[1,1]*gest,'-','LineWidth',1,'Color',[0.65,0.65,0.65])
    plot(xlim,[1,1]*(gest+gstd),'--','LineWidth',1,'Color',[0.75,0.75,0.75])
    plot(xlim,[1,1]*(gest-gstd),'--','LineWidth',1,'Color',[0.75,0.75,0.75])
    text(sum([0.95,0.05].*xlim),gest,...
        sprintf('g_{est}=%.3f hr^{-1}',gest),'FontSize',12,...
        'VerticalAlignment','middle','HorizontalAlignment','left')
    hold off

    subplot(2,2,1:2)
    hold off
    semilogy(ts,OD,'LineWidth',2)

    if gest ~= 0
        centerOD = exp((log(min(MinUsedOD,max(OD)/10)) + log(max(OD)/ExcludedEndFactor))/2);
        centerT = ts(find(medfilt1(OD,round(1/dt+1)) >= centerOD,1));

        offsetDir = -0.5*(log(ylim)*[-1;1])/(xlim*[-1;1])/gest;
        offsetCenter = [centerT+0.15/offsetDir-0.3,...
            10^(log10(centerOD)+0.15-0.3*offsetDir)];

        fitPlotLength = ...
            min(ts(find(medfilt1(OD,round(1/dt+1)) >= max(OD)/ExcludedEndFactor,1)) - ...
            ts(find(medfilt1(OD,round(1/dt+1)) >= min(MinUsedOD,max(OD)/10),1)),...
            2*offsetCenter(1)-1);
    else
        offsetCenter = [mean(ts), 10^(log10(mean(OD))+0.5)];
        fitPlotLength = max(ts)/2;
    end

    hold on
    if gest ~= 0
        plot(offsetCenter(1)+[-0.5,0.5]*fitPlotLength,...
            offsetCenter(2)*exp([-0.5,0.5]*fitPlotLength*(gest-gstd)),...
            '--','LineWidth',1.5,'Color',[0.5,0.5,0.5])
    end
    plot(offsetCenter(1)+[-0.5,0.5]*fitPlotLength,...
        offsetCenter(2)*exp([-0.5,0.5]*fitPlotLength*(gest+gstd)),...
        '--','LineWidth',1.5,'Color',[0.5,0.5,0.5])

    plot(offsetCenter(1)+[-0.5,0.5]*fitPlotLength,...
        offsetCenter(2)*exp([-0.5,0.5]*fitPlotLength*gest),...
        'LineWidth',2.5,'Color',[0.2,0.2,0.2])
    hold off

    ylabel('Pop. Size (OD_{600})')
    xlabel('Time (hours)')
    if ~strcmp(TitleString,"")
        title(sprintf('%s:  g_{est} = %.3f +/- %.3f hr^{-1}',TitleString,gest,gstd))
    else
        title(sprintf('g_{est} = %.3f +/- %.3f hr^{-1}',gest,gstd))
    end
    set(gca,'FontSize',16)
end

end

%out = [nargin,length(varargin)];