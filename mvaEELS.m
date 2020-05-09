%MVAEELS 
% mvaEELS is a Matlab script that reads STEM-EELS datacubes from .dm3 
% files, conducts background subtraction, separates the signal components 
% by NNMF, and performs signal quantification.
% Author: Siyuan Zhang (<a
% href="mailto:siyuan.zhang@mpie.de">siyuan.zhang@mpie.de</a>)
% Please cite this paper: https://doi.org/10.1093/jmicro/dfx091

function output = mvaEELS(input)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',12)
if (nargin == 1) && isfield(input,'data') %have input data in proper format
    data = input.data;
    output = input;
else %new input
    [e,h,x,y] = readDualEELS;
    data.d3 = reshape(h',length(e),length(x),length(y));
    data.e = e;
    data.x = x;
    data.y = y;
    output.data = data;
end
%% Noise reduction
if isfield(output,'nComp') %have denoised input data
    si = output;
else
    sumS = sum(squeeze(sum(data.d3,2)),2);
    h = reshape(data.d3,length(data.e),[]);
    winEELS = selectWindow(data.e,sumS,'Select energy range with signal & background');
    e = data.e(winEELS);
    eels = h(winEELS,:);
    defaultE = {['E',num2str(floor(min(e))),'-',num2str(ceil(max(e)))]};
    str = inputdlg('Give a simple name','MVA analysis range',1,defaultE); %Name the energy range
    str = str{1};
    if ~isvarname(str)
        return
    end
    if size(eels,2) > 10
        [s,w,x,y,scree] = mva(eels,e,data.x,data.y,0,-2);%-2:PCA only
        data.d3 = reshape(s*w,length(e),length(x),[]);
        data.e = e;
        data.x = x;
        data.y = y;
        si.data = data;
        si.nComp = size(s,2);
        si.pca.spect = s;%chang to spect(e,nComp);
        w = reshape(w,si.nComp,length(x),length(y));%change to weight(nComp,x,y);
        si.pca.weight = permute(w,[3,2,1]);%change to weight(y,x,nComp);
        si.pca.scree = scree;
    end
    output.(str) = si;
end
%% Background model and removal
if ~isfield(si,'Sig')
    [si.BG,si.nmf] = eelsBG(si);
end
if isempty(si.nmf)
    return
end
%% Process signals
iSig = 1;
while menu(sprintf('Evaluate signal %d',iSig),{'Yes','Quit'}) == 1
    si.(['Sig',num2str(iSig)]) = eelsSignal(si.nmf);
    iSig = iSig + 1;
end
if exist('str','var')
    output.(str) = si;
else 
    output = si;
end
end

%% Background modelling
function [bg,signal] = eelsBG(input)
nComp = input.nComp;
E = input.data.e;
x = input.data.x;
y = input.data.y;
Nx = length(x);
Ny = length(y);
nPix = Nx*Ny;
H = input.data.d3;
if length(size(H)) == 3
    H = reshape(H,length(E),[]);
end
if nComp > 0
    sOption = menu('Background removal for each',{'Component','Factor','Pixel','Quit'});
else
    sOption = menu('Background removal',{'Proceed','Quit'})+2;
end
if sOption == 4
    bg = [];
    signal = {};
    return
end
if sOption == 1
    sStr = 'Comp';
    defBG = 'power';
    [nmfS,nmfW] = mva(H,E,x,y,0,nComp);
elseif sOption == 2
    sStr = 'Fact';
    defBG = 'logx';
    logH = log(H);
    [nmfS,nmfW] = mva(logH,E,x,y,0,nComp);
end
if sOption < 3
    nComp = size(nmfS,2);
    nmf.nComp = nComp;
    nmf.method = sStr;
    nmf.spect = nmfS;
    nmf.weight = reshape(nmfW,nComp,Nx,Ny);
    nmf.weight = permute(nmf.weight,[3,2,1]);
    nmfBG = nmfS;
    nmfSig = nmfS;
    winBGs = cell(1,nComp);
    typeBG = cell(1,nComp);
    aBG = zeros(nComp,1);
    bBG = zeros(nComp,1);
    checkSig = 1;
    while checkSig
        for iComp = 1:nComp
            if isempty(winBGs{iComp})
                winBGs{iComp} = winBGs{1};
            end
            [bg,aBG(iComp),bBG(iComp),sig,winBGs{iComp},typeBG{iComp}] = ...
                eelsBackground(E,nmfS(:,iComp),sprintf('Component %d background',iComp),winBGs{iComp},defBG);
            nmfBG(:,iComp) = bg;
            nmfSig(:,iComp) = sig;
        end
        if sOption == 1
            eelsSig = H - nmfBG*nmfW;
        elseif sOption == 2
            eelsSig = H - exp(nmfBG*nmfW);
        end
        figID = figure('Name','Check background subtraction','KeyPressFcn',@(obj,evt) 0);
        nPix = size(eelsSig,2);
        if nPix > 10
            ySig = eelsSig(:,ceil(rand(10,1)*nPix));
        else
            ySig = eelsSig;
        end
        plot(E,ySig);
        axis tight
        set(gca,'lineWidth',2)
        title('[Q] Proceed, [Other keys] Adjust background model');
        waitfor(gcf,'CurrentCharacter');
        key = uint8(get(gcf,'CurrentCharacter'));
        set(gcf,'CurrentCharacter',char(0))
        if key == 113;
            checkSig = 0;
        end
        close(figID)
    end
    winSG = 1:length(E);
    nmf.sig = nmfSig;
    nmf.backg = nmfBG;
    for iComp = 1:nComp
        nmf.(['BG',num2str(iComp)]).window = E(winBGs{iComp});
        winSG = winSG(~ismember(winSG,winBGs{iComp}));
        nmf.(['BG',num2str(iComp)]).type = typeBG{iComp};
        nmf.(['BG',num2str(iComp)]).r = aBG(iComp);
        nmf.(['BG',num2str(iComp)]).A = bBG(iComp);
    end
    if sOption == 2
        nmf.BG.r = reshape(aBG'*nmfW,Nx,Ny)';
        nmf.BG.A = reshape(bBG'*nmfW,Nx,Ny)';
    end
    bg = nmf;
    eelsSig = eelsSig(winSG,:);
    data.d3 = reshape(eelsSig,length(winSG),Nx,Ny);
    data.e = E(winSG);
    data.x = x;
    data.y = y;
    signal.data = data;
elseif sOption == 3
    sumH = sum(H,2);
    [~,aBG,bBG,sig,winBG,typeBG] = eelsBackground(E,sumH,'Select background energy window');
    eelsSig = sig';
    if nPix > 1
        aBG = zeros(nPix,1);
        bBG = zeros(nPix,1);
        eelsSig = zeros(length(E),nPix);
        for iPix = 1:nPix
            [~,aBG(iPix),bBG(iPix),sig] = eelsBackground(E,H(:,iPix),0,winBG,typeBG);
            eelsSig(:,iPix) = sig';
        end
    end
    bg.window = E(winBG);
    winSG = 1:length(E);
    winSG = winSG(~ismember(winSG,winBG));
    bg.type = typeBG;
    bg.r = reshape(aBG,Nx,Ny)';
    bg.A = reshape(bBG,Nx,Ny)';
    eelsSig = eelsSig(winSG,:);
    data.d3 = reshape(eelsSig,length(winSG),Nx,Ny);
    data.e = E(winSG);
    data.x = x;
    data.y = y;
    signal.data = data;
end
%% NNMF after background removal
eelsSig(eelsSig<0) = 0; %Negative counts removed for NNMF!
[nmfS,nmfW,nmfX,nmfY,nmfScree] = mva(eelsSig,data.e,x,y,0,-1); %-1:PCA&NNMF
nComp = size(nmfW,1);
signal.nComp = nComp;
signal.spect = nmfS;
signal.weight = reshape(nmfW,nComp,length(nmfX),length(nmfY));
signal.weight = permute(signal.weight,[3,2,1]);
signal.scree = nmfScree;
end

%% Background removal of individual EELS spectra
function [BG,a,b,sig,win,sType] = eelsBackground(e,spect,winName,win,sType)
sig = 'N/A';
if nargin < 2
    readEELS = 1;
    while readEELS
        [spect,xyz,~,expara] = readDM;
        e = xyz{end};
        readEELS = isempty(expara);
        spect = spect/expara(3)/expara(4);
    end
end
if nargin < 3
    winName = '';
end
if nargin < 4
    win = selectWindow(e,spect,'Background window');
end
if isempty(win)
    win = selectWindow(e,spect,'Background window');
end
if nargin < 5
    sType = 'power';
end
if size(e,1) > size(e,2)
    e = e';
end
if size(e,1) ~= 1
    e = e(1,:);
end
if size(spect,2) ~= size(e,2)
    if size(spect,1) == size(e,2)
        spect = spect';
    else
        error('Input dimension of energy and spectrum does not match.')
    end
end

if winName == 0 % Automatic
    BG = spect;
    a = zeros(1,size(spect,1));
    b = zeros(1,size(spect,1));    
    for i = 1:size(spect,1)
        [BG(i,:),a(i),b(i)] = fitEELS(e,spect(i,:),win,'p');
        sig = spect - BG;
    end
else % Interactive
if size(spect,1) ~= 1
    spect = sum(spect);
end
screenSize = get(0,'ScreenSize');
figID = figure('Name',winName,'KeyPressFcn',@(obj,evt) 0,'Position',[1 screenSize(4)/2.5 screenSize(3) screenSize(4)/2]);
sSwitch = 1;
BG = NaN;
a = 0;
while sSwitch
    figure(figID)
    switch sType
        case 'power'
            [BG,a,b] = fitEELS(e,spect,win,'p',a);
            sBG = sprintf('x^%.3g*exp(%.3g)',a,b);
        case 'logx'
            [BG,a,b] = fitEELS(e,spect,win,'x',a);
            sBG = sprintf('%.3g*log(x)+%.3g',a,b);
        case 'logy'
            [BG,a,b] = fitEELS(e,spect,win,'y',a);
            sBG = sprintf('exp(%.3g*x+%.3g)',a,b);
        case 'linear'
            [BG,a,b] = fitEELS(e,spect,win,'l',a);
            sBG = sprintf('%.3g*x+%.3g',a,b);
        case 7
            if isnan(sum(BG))
                BG = min(spect(win))*ones(1,length(e));
            else
                BG = (9*BG+min(spect(win))*ones(1,length(e)))/10;
            end
            sBG = sprintf('%.2g',BG(1));
        case 'constant'
            BG = b*ones(1,length(e));
            sBG = sprintf('%.2g',BG(1));
        case 'signal'
            BG = spect;
            sBG = sprintf('no Signal (BG = Signal)');
    end
    sig = spect - BG;
    winBG = zeros(1,length(e));
    winBG(win) = max(sig);
    subplot(1,2,1)
    plot(e,spect,e,BG,'lineWidth',2)
    axis tight
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    iStr = sprintf('\nBackground: [S] signal, [P] power, [L] linear, [X] logx, [Y] logy\nConstant background: [0] no, [M] maximum, [N] minimum, [C] mean\nBackground windows: [W] choose windows, [Left/Right Arrow] adjust\n[Up/Down Arrow] adjust background\n[Q] quit');
    title(sprintf(['Background [%.2f,%.2f]: ',sBG,iStr],e(min(win)),e(max(win))));
    legend('Spectrum','Background')
    subplot(1,2,2)
    plot(e,sig,e,winBG,'lineWidth',2)
    axis tight
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    iStr = sprintf('\nBackground windows: [W] choose windows, [Left/Right Arrow] adjust\n[Up/Down Arrow] adjust background\n[Q] quit');
    title(sprintf(['Background [%.2f,%.2f]: ',sBG,iStr],e(min(win)),e(max(win))));
    legend('Signal','Window range')
    waitfor(gcf,'CurrentCharacter');
    key = uint8(get(gcf,'CurrentCharacter'));
    set(gcf,'CurrentCharacter',char(0))
    if key == 30 % Up arrow: enhance background
        if strcmp(sType,'power')||strcmp(sType,'logx')||strcmp(sType,'logy')||strcmp(sType,'linear')
            a = ceil(a*100+1)/100;
        elseif strcmp(sType,'constant')
            b = b + (max(spect(win))-min(spect(win)))/100;
        end
    elseif key == 31 % Down arrow: reduce background
        if strcmp(sType,'power')||strcmp(sType,'logx')||strcmp(sType,'logy')||strcmp(sType,'linear')
            a = floor(a*100-1)/100;
        elseif strcmp(sType,'constant')
            b = b - (max(spect(win))-min(spect(win)))/100;
        end
    elseif key == 113 % Q: quit
        sSwitch = 0;
    elseif key == 115 % S: signal
        sType = 'signal';
    else
        a = 0;
    end
    if key == 28 % Left arrow: adjust right window boundary
        if length(win) > 1
            win(end) = [];
        else
            key = 119;
        end
    elseif key == 29 % Right arrow: adjust left window boundary
        if length(win) > 1
            win(1) = [];
        else
            key = 119;
        end
    end
    if key == 119 % W: window selection
        win = selectWindow(e,spect,'Background window');
    elseif key == 112 % P: power
        sType = 'power';
    elseif key == 120 % X: log(x)
        sType = 'logx';
    elseif key == 121 % Y: log(y)
        sType = 'logy';
    elseif key == 108 % L: linear
        sType = 'linear';
    elseif key == 110 % N: minimum (constant)
        sType = 'constant';
        b = min(spect(win));
    elseif key == 99 % C: mean (constant)
        sType = 'constant';
        b = mean(spect(win));
    elseif key == 109 % M: maximum (constant)
        sType = 'constant';
        b = max(spect(win));
    elseif key == 48 % 0
        sType = 'constant';
        b = 0;
    elseif key == 117 % U: undefined
        sig = 'N/A';
        sType = 'undefined';
        BG = NaN;
        b = 0;
        sSwitch = 0;
    end
end
close(figID)
end
end

function [BG,a,b,sig] = fitEELS(x0,y0,xWin,typeBG,a,b)
if nargin < 4
    typeBG = 'p';
end
x = x0(xWin);
y = y0(xWin);
if typeBG == 'p'
    y = log(y);
    x = log(x);
elseif typeBG == 'x'
    x = log(x);
elseif typeBG == 'y'
    y = log(y);
end
if nargin < 5
    p = real(polyfit(x,y,1));
    a = p(1);
    b = p(2);
elseif nargin < 6
    if a == 0
        p = real(polyfit(x,y,1));
        a = p(1);
        b = p(2);
    else
        b = y(end)-a*x(end);
    end
end
if typeBG == 'p'
    BG = x0.^a*exp(b);
elseif typeBG == 'l'
    BG = a*x0+b;
elseif typeBG == 'x'
    BG = a*log(x0)+b;
elseif typeBG == 'y'
    BG = exp(a*x0+b);
end
sig = y0-BG;
end

%% EELS signal analysis
function output = eelsSignal(input)
eSig = input.data.e;
Nx = length(input.data.x);
Ny = length(input.data.y);
nComp = input.nComp;
nmfSig = input.spect;
nmfW = permute(input.weight,[3,2,1]);%weight(y,x,nComp) to (nComp,x,y)
nmfW = reshape(nmfW,nComp,Nx*Ny);
%% Signal processing
winSig = cell(nComp,1);
for iComp = 1:nComp
    winSigA = 1;
    winSigB = length(eSig);
    if iComp > 1
        winSigA = min(winSig{1});
        winSigB = max(winSig{1});
    end
    winSig{iComp} = selectWindow(eSig,nmfSig(:,iComp),sprintf('Select signal windows for component %d',iComp),winSigA,winSigB);
end
%% Signal integration
nmfI = zeros(1,nComp);
nmfIx = zeros(1,nComp);
nmfIx2 = zeros(1,nComp);
for iComp = 1:nComp
    win = winSig{iComp};
    nmfSig(win,iComp) = nmfSig(win,iComp)-nmfSig(win(1),iComp); %% Remove constant background!
    nmfI(iComp) = sum(nmfSig(winSig{iComp},iComp))';
    nmfIx(iComp) = sum(nmfSig(winSig{iComp},iComp).*eSig(winSig{iComp})')';
    nmfIx2(iComp) = sum(nmfSig(winSig{iComp},iComp).*eSig(winSig{iComp})'.^2)';
end
output.sum = reshape(nmfI*nmfW,Nx,Ny)';
mean = (nmfIx*nmfW)./(nmfI*nmfW);
output.mean = reshape(mean,Nx,Ny)';
output.var = reshape((nmfIx2*nmfW)./(nmfI*nmfW)-mean.^2,Nx,Ny)';
output.sum1 = reshape(nmfIx*nmfW,Nx,Ny)';
output.sum2 = reshape(nmfIx2*nmfW,Nx,Ny)';

comp.sum = nmfI;
mean = nmfIx./nmfI;
comp.mean = mean;
comp.var = nmfIx2./nmfI-mean.^2;
comp.sum1 = nmfIx;
comp.sum2 = nmfIx2;
for iComp = 1:nComp
    comp.(['window',num2str(iComp)]) = eSig(winSig{iComp});
end
output.comp = comp;
end

function [winAll,win] = selectWindow(x,y,winName,win1,win2)
if nargin < 2 %% original
    y = x;
    x = 1:length(x);
end
if size(x,1) == 1
    x = x';
end
if nargin < 3
    winName = '';
end
if nargin < 4
    win1 = 1;
end
if nargin < 5
    win2 = length(x);
end
if length(x) == size(y,1)
elseif length(x) == size(y,2)
    y = y';
else
    error('Input vectors should have the same length!')
end
if win2 == 0
    winAll = [];
    win = {};
else
figID = figure('Name',winName,'WindowKeyPressFcn',@(obj,evt) 0);
sSwitch = 1;
sStep = 1;
sX = 0;
sY = 0;
win = {};
iWin = 1;
while sSwitch
    figure(figID)
    if sY
        semilogy(x(win1:win2),y(win1:win2,:),'lineWidth',2)
    else
        plot(x(win1:win2),y(win1:win2,:),'lineWidth',2)
    end
    if sX
        I = sum(y(win1:win2,:));
        Ix = sum(y(win1:win2,:).*repmat(x(win1:win2),1,1));
        Ix2 = sum(y(win1:win2,:).*repmat(x(win1:win2).^2,1,1));
        xM = Ix./I;
        xQ = real(sqrt(Ix2./I-xM.^2));
        [~,xm] = min(abs(x-xM));
        [~,xm1] = min(abs(x-xM+xQ));
        [~,xm2] = min(abs(x-xM-xQ));
        hold on
        stem([xM,xM-xQ,xM+xQ],y([xm,xm1,xm2]))
        hold off
    end
    axis tight
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    if sSwitch == 1
        tempStr1 = 'left';
        tempStr2 = 'right';
    elseif sSwitch == 2
        tempStr1 = 'right';
        tempStr2 = 'left';
    end
    title([sprintf('window %d: [%.2f, %.2f]\n',iWin,x(win1),x(win2)),...
        '[\leftarrow],[\rightarrow]:Move the ',tempStr1,sprintf(' boundary\n'),...
        '[Enter]:Control the ',tempStr2,sprintf(' boundary\n'),...
        '[\uparrow],[\downarrow]:Change step size ',sprintf('%d\n',sStep),...
        '[Q]:Proceed, ',sprintf('[N]:Select window %d',iWin+1)]);
    waitfor(gcf,'CurrentCharacter');
    key = uint8(get(gcf,'CurrentCharacter'));
    set(gcf,'CurrentCharacter',char(0))
    if key == 13 % Enter: select/deselect
        sSwitch = 3 - sSwitch;
    elseif key == 28 % Left arrow
        if sSwitch == 1
            win1 = max(win1-sStep,1);
        elseif sSwitch == 2
            win2 = max(win2-sStep,win1+1);
        end
    elseif key == 29 % Right arrow
        if sSwitch == 1
            win1 = min(win1+sStep,win2-1);
        elseif sSwitch == 2
            win2 = min(win2+sStep,length(x));
        end
    elseif key == 30 % Up arrow
        sStep = min(sStep*10,ceil(length(x)/2));
    elseif key == 31 % Down arrow
        sStep = ceil(sStep/10);
    elseif key == 120 % X: X-axis mean
        if size(y,2) == 1
            sX = ~sX;
        end
    elseif key == 121 % Y: Semilogy
        sY = ~sY;
    elseif key == 113 % Q: quit
        win{iWin} = win1:win2;
        iWin = iWin+1;
        sSwitch = 0;
    elseif key == 110 % N: next window selection loop
        win{iWin} = win1:win2;
        iWin = iWin+1;
    end
end
close(figID)
nWin = iWin-1;
winAll = [];
for iWin = 1:nWin
    winAll = [winAll,win{iWin}];
end
end
end