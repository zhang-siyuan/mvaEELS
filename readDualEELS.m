%READDUALEELS 
% readDualEELS is a Matlab script that reads and combines high and low loss
% EELS datacubes from dm3, and centres the zero-loss channels for output. 
% Author: Siyuan Zhang (<a
% href="mailto:siyuan.zhang@mpie.de">siyuan.zhang@mpie.de</a>)
% Please cite this paper: https://doi.org/10.1093/jmicro/dfx091

function [eSum,eelsSum,x,y,zlCor] = readDualEELS
filePath = pwd;
readHL = 1;
while readHL
    disp('read high/low loss data')
    [eels1,axes,filename,expara] = readDM(filePath);
    fileSep = max(strfind(filename,'\'));
    filePath = filename(1:fileSep);
    fileName = filename(fileSep+1:end);
    e1 = axes{end};
    readHL = isempty(expara);%||~expara(1)||expara(2);
end
eels1 = eels1/expara(3)/expara(4);
if ~expara(1)||expara(2) 
    disp('low loss data read: self calibration')
    e0 = e1;
    eels0 = eels1;
else
    readLL = 1;
    while readLL
        fileName = strrep(fileName,'high','low');
        if exist([filePath,fileName],'file')
            [eels0,e,~,expara] = readDM(filePath,fileName);
        else
            disp('high loss data read: read corresponding zero loss data')
            [eels0,e,~,expara] = readDM(filePath);            
        end
        e0 = e{end};
        readLL = isempty(expara)||~expara(1)||~expara(2);
    end
    eels0 = eels0/expara(3)/expara(4);
end

if length(axes) == 1
    x = 0;
    y = 0;
elseif length(axes) == 2
    x = axes{1};
    y = 0;
elseif length(axes) == 3
    x = axes{1};
    y = axes{2};
end

nP = size(eels0,1);
zlP = zeros(1,nP);
for i = 1:nP
    zlP(i) = eelsZL(e0,eels0(i,:));
end
maxZLP = max(zlP);
minZLP = min(zlP);
oZLP = find(e0>=0,1);
zlCor = oZLP - zlP;
eelsL = eels0;
eelsC = eels1;
for i = 1:nP
    eelsL(i,:) = circshift(eels0(i,:),[0,zlCor(i)]);
    eelsC(i,:) = circshift(eels1(i,:),[0,zlCor(i)]);
end
eL = e0;
eC = e1;
eL([1:oZLP-minZLP,end-maxZLP+oZLP:end]) = [];
eC([1:oZLP-minZLP,end-maxZLP+oZLP:end]) = [];
eelsL(:,[1:oZLP-minZLP,end-maxZLP+oZLP:end]) = [];
eelsC(:,[1:oZLP-minZLP,end-maxZLP+oZLP:end]) = [];
[~,eCut] = find(eL>=min(eC));
eL(eCut) = [];
eelsL(:,eCut) = [];
eSum = [eL,eC];
eelsSum = [eelsL,eelsC];
end

function nZL = eelsZL(e,eels)
dE = e(2)-e(1); %new
nE= length(e);
minC = min(eels);
eels = eels - minC;
[~,eMax] = max(eels); %new
nLeft = max(1,floor(eMax-2/dE)); %new 
if isempty(nLeft)
    nZL = 1;
else
    w = 5;
    iW = 1;
    f = eels(nLeft:nE);
    while sum(f(iW:iW+w-1)) < sum(f(iW+w:iW+2*w-1))
        iW = iW+1;
    end
    [~,wM] = max(f(iW:iW+2*w-1));
    nZL = nLeft+wM+iW-2;
end
end

function [dmCube,dmAxis,filename,eels,dmData] = readDM(filePath,fileName,logfile)
% function [dmCube,dmAxis,filePath,eels,dmData] = readDM(filePath,fileName) 
% read a Digital Micrograph file version 3 or 4 and return the first image
% in its native data format (e.g. int16), along with its pixel scale 
% information. If a logfile name is specified, a description of the entire 
% tree will be written to the console and also to the file.

% Siyuan Zhang version to include dualEELS readouts, October 2016
% Based on ReadDMFile F. Sigworth, August 2013
% Based on ReadDM3 F. Sigworth, July 2009
% Code was based on the description by Greg Jefferis at
% http://rsb.info.nih.gov/ij/plugins/DM3Format.gj.html
% Version 4 differs from version 3 mainly in that it uses uint64 instead of
% int32 or uint32 for various count entries.  The function GetLLong is used
% to get these counts in a version-dependent manner.
% 
% This function has been written assuming that it will run on a
% little-endian (e.g. Intel) machine, reading a file written by a
% little-endan machine. Otherwise the use of the swapbytes function will
% have to be made conditional on the local machine type and the byteorder
% variable. Also, swapbytes may have to be added to the GetData function.
% 
% Note that the code here allows any fields from the tree to be extracted. 
% Here is where we define the fields that we will extract. We grab the data 
% value at the first match of each of these tags.  Here numerals represent 
% unnamed fields. To see what all the tags are, specify a logfile to 
% receive the hundreds of lines of information!

if nargin < 1
    [fileName,filePath] = uigetfile('*.dm3');
elseif nargin == 1
    [fileName,filePath] = uigetfile([filePath,'/.dm3']);
end
filename = [filePath,fileName];

celltags={'ImageList 2 ImageData Calibrations Dimension 1 Scale'
    'ImageList 2 ImageData Calibrations Dimension 1 Origin'
    'ImageList 2 ImageData Dimensions 1'
    'ImageList 2 ImageData Calibrations Dimension 2 Scale'
    'ImageList 2 ImageData Calibrations Dimension 2 Origin'
    'ImageList 2 ImageData Dimensions 2'
    'ImageList 2 ImageData Calibrations Dimension 3 Scale'
    'ImageList 2 ImageData Calibrations Dimension 3 Origin'
    'ImageList 2 ImageData Dimensions 3'
    'ImageList 2 ImageData Data'
    'ImageList 2 ImageTags EELS Acquisition Dual_acquire_enabled'
    'ImageList 2 ImageTags EELS Acquisition Is_dual_acquire_low-loss'
    'ImageList 2 ImageTags EELS Acquisition Exposure_(s)'
    'ImageList 2 ImageTags EELS Acquisition Number_of_frames'};

% Strings to write out data types
dataTypeString={'1?'      'int16'   'int32' 'uint16' 'uint32'...
                'float32' 'float64' 'bool'  'int8'   'uint8'...
                'int64'   'uint64'  '13?'   '14?'    'struct'...
                '16?'     '17?'     'string' '20?'   'array'};

found=zeros(size(celltags));
output=cell(size(celltags));

% Set up the log file.  We use my mprintf function to allow printing to the
% console as well, and suppressing printing when the handle is zero.
if nargin > 2
    flog=fopen(logfile,'w');
    hs = [1 flog];  % log file handles
else
    flog=0;
    hs = 0;
end
tabstring='| ';
level=0;
maxprint=4;
OutputOn=1;

% Read the whole file into memory as a byte array.
fb=fopen(filename,'r');
tic;
d=fread(fb,inf,'*uint8');  % read the whole file as bytes
dt=toc;
fclose(fb);

p=uint64(1);  % byte pointer--also a global variable
Tags=cell(1,10); % Keeps track of the sequence of tags, for matching with the tag strings.

% Pick up the header
version=GetLong;
if (version<3) || (version>4)
    error(['ReadDM34: Wrong file type.  Version = ' num2str(version)]);
end;
mprintf(hs,'Version %d\n',version);

nbytes=GetLLong;
mprintf(hs,'Total size %d MB\n',nbytes/2^20);
mprintf(hs,'Read in %d s\n',dt);
% Handle little- and big-endian files and machines
dle=GetLong;  % boolean to tell whether data are little endian
[~,~,endian] = computer; % [str,maxsize,endian]
mle= (endian=='L');  % machine is little endian: we'll have to swap bytes in reading the tree.
dswap=(dle~=mle);  % swap byte-order when reading data

% Traverse the tree
GetTagGroup;

% Close the logfile
if flog>0
    fclose(flog);
end

% Extract the output parameters
if isempty(output{1})
else
    dx=output{1};
    x0=output{2};
    x0=-x0*dx;
    xdim=output{3};
    xAxis = x0:dx:x0+dx*double(xdim-1);
    dmData.xAxis = double(xAxis);
    dmAxis{1} = dmData.xAxis;
    if isempty(output{4})
        dmData.cube = double(output{10});
        dmCube = double(output{10})';
    else
        dy=output{4};
        y0=output{5};
        y0=-y0*dy;
        ydim=output{6};
        yAxis = y0:dy:y0+dy*double(ydim-1);
        dmData.yAxis = double(yAxis);
        dmAxis{2} = dmData.yAxis;
        if isempty(output{7})
            dmData.cube = reshape(double(output{10}),[xdim ydim]);
            dmCube = reshape(double(output{10}),[xdim ydim]);
        else
            dz=output{7};
            z0=output{8};
            z0=-z0*dz;
            zdim=output{9};         % no. of images in a stack
            zAxis = z0:dz:z0+dz*double(zdim-1);
            dmData.zAxis = double(zAxis);
            dmAxis{3} = dmData.zAxis;
            dmData.cube = reshape(double(output{10}),[xdim ydim zdim]);
            dmCube = reshape(double(output{10}),[xdim*ydim,zdim]);
        end
    end
end
if isempty(output{13})
    eels = [];
else
    eels = zeros(1,4);
    eels(1) = output{11};
    eels(2) = output{12};
    eels(3) = output{13};
    eels(4) = output{14};
end
% end of main function
% -------------------------------

% ---- here are all the local functions, called recursively ----

    function GetTagGroup
        sorted=GetByte;
        open=GetByte;
        NumTags=GetLLong;
        for i=1:NumTags
            GetTagEntry(i);
        end
    end

    function GetTagEntry(MemberIndex)
        level=level+1;
        PutNewline;
        PutTabs;
        isdata=GetByte;
        labelsize=GetInt;
        labelstring=GetString(labelsize);
        PutStr('-');
        PutStr([labelstring ':']);
        if numel(labelstring)<1
            labelstring=num2str(MemberIndex);
        end
        Tags{level}=labelstring;
        if version==4
            totalBytes=GetLLong;
        end
        if isdata==21
            GetTagType
        elseif isdata==20
            GetTagGroup
        else
            error(['Unknown TagEntry type ' num2str(isdata)]);
        end
        Tags{level}=[];
        level=level-1;
    end

    function GetTagType
        dum=GetLong;
        if dum ~= 623191333
            disp(['Illegal TagType value ' num2str(dum)]);
        end
        deflen=GetLLong;  % number of items in encoding array
        EncType=GetLLong;
%         for i=1:deflen
%             EncType(i)=GetLLong;  % Don't know how to use the array...
%         end;
        x=GetData(EncType);
        index=CheckTags;
        if index>0
            output{index}=x;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is the function that slows everything down, which checks the
    % entire tag list for each tagged type.  It is inefficient, but simple.
    function r=CheckTags
        for i=1:numel(celltags)
            ok=~found(i);
            c=celltags{i};
            j=1;
            while ok && (numel(c)>0)
                [s c]=strtok(c);
                ok=strcmp(s,'*')||strcmp(s,Tags{j});
                j=j+1;
            end
            if ok
                r=i;
                return
            end
        end
        r=0;
    end


    function x=GetData(ftype,num)
        if nargin < 2
            num=1;
        end
        num=uint64(num);
        x=[];
        %         disp(['GetData ' num2str(ftype)]);
        switch ftype
            case 2  % short
                x=typecast(d(p:p+num*2-1),'int16');
                p=p+2*num;
            case 3  % long
                x=typecast(d(p:p+num*4-1),'int32');
                p=p+4*num;
            case 4  % ushort
                x=typecast(d(p:p+num*2-1),'uint16');
                p=p+2*num;
            case 5  % ulong
                x=typecast(d(p:p+num*4-1),'uint32');
                p=p+4*num;
            case 6  % float
                x=typecast(d(p:p+num*4-1),'single');  % Takes 4 s for 142 M elements
                p=p+4*num;
            case 7  % double
                x=typecast(d(p:p+num*8-1),'double');
                p=p+8*num;
            case 8  % boolean
                x=d(p:p+num-1);
                p=p+num;
            case 9  % char
                x=char(d(p:p+num-1));
                p=p+num;
            case 10  % octet
                x=(d(p:p+num-1));
                p=p+num;
            case 11  % int64
                x=typecast(d(p:p+num*8-1),'int64');
                p=p+8*num;
            case 12  % uint64
                x=typecast(d(p:p+num*8-1),'uint64');
                p=p+8*num;
            case 15  % Struct
                PutStr('struct');
                StructNameLength=GetLLong;
                NumFields=GetLLong;
                x=[];
                for i=1:NumFields
                    FieldNameLength(i)=GetLLong;
                    FieldType(i)=GetLLong;
                end;
                StructName=GetString(StructNameLength);
                PutStr(StructName);
                PutNewline; PutTabs;
                for i=1:NumFields
                    %                     FieldNameLen=FieldNameLength(i);
                    FieldName=GetString(FieldNameLength(i));
                    PutStr(FieldName);
                    x(i)=GetData(FieldType(i));
                    PutNewline; PutTabs;
                end;
            case 18 % string
                length=GetLong;
                x=char(d(p:p+length-1)');
                PutVal(x); PutNewline;
                p=p+uint64(length);
                
            case 20  % Array
                ArrayType=GetLLong;
                if ArrayType==15  % Struct is special case
                    StructNameLength=GetLLong;
                    NumFields=GetLLong;
                    x=[];
                    for i=1:NumFields
                        FieldNameLength(i)=GetLLong;
                        FieldType(i)=GetLLong;
                    end;
                end;
                ArrayLength=GetLLong;
                
                if ArrayType ~=4
                    PutStr('array of');
                    PutVal(ArrayLength);
                    PutStr(' --type: ');
                    PutStr(dataTypeString{ArrayType});
%                     PutVal(ArrayType);
                end;
                
                if ArrayType==15
                    PutStr('structs');
                    PutNewline;
                    for j=1:ArrayLength
                        OutputOn=j<=maxprint;
                        for i=1:NumFields
                            FieldNameLen=FieldNameLength(i);
                            FieldName=GetString(FieldNameLength(i));
                            FieldTy=FieldType(i);
                            PutTabs;
                            PutStr(FieldName);
                            x(i)=GetData(FieldType(i));
                            PutNewline;
                        end;
                        OutputOn=1;
                    end;
                elseif ArrayType==4
                    OutputOn=0;
                    for j=1:ArrayLength
                        x(j)=GetData(ArrayType);
                    end;
                    OutputOn=1;
                    PutVal(char(x'));
                else
                    % Might be long data
                    if (ArrayLength > 1000)  % try to handle a long array
                        OutputOn=0;
                        x=GetData(ArrayType,ArrayLength);
                        OutputOn=1;
                    else
                        PutNewline;
                        for j=1:ArrayLength
                            OutputOn=j<=maxprint;
                            PutTabs;
                            x(j)=GetData(ArrayType);
                            PutNewline;
                        end;
                        OutputOn=1;
                    end; % long data
                end;
            otherwise
                x=0;
                disp(['Unrecognized data type ' num2str(ftype)]);
        end; % switch
        if (ftype < 15) && OutputOn
            PutVal(x);
        end;
    end % GetData


    function PutStr(s)
        if OutputOn
            mprintf(hs,'%s ',s);
        end
    end

    function PutTabs
        if OutputOn
            for i=1:level-1
                mprintf(hs,'%s',tabstring);
            end
        end
    end

    function PutVal(x)
        if OutputOn
            if isa(x,'char')
                mprintf(hs,'%s ',x);
            else
                mprintf(hs,'%d ',x);
            end
        end;
    end

    function PutNewline
        if OutputOn
            mprintf(hs,'\n');
        end
    end

    function s=GetString(len)
        len=uint64(len);
        if len<1
            s=[];
        else
            stemp = d(p:p+len-1)';
            stemp(stemp<=32) = 95;
            s = char(stemp);
            %s=char(d(p:p+len-1)');
            p=p+len;
        end;
    end

    function x=GetLong  
        x=typecast(d(p:p+3),'int32');
        x=swapbytes(x);
        p=p+4;        
    end

    function x=GetLLong % uint32 in version 3, uint64 in version 4
        if version == 3
            x=GetLong;
        else  % version 4 code:
            x=typecast(d(p:p+7),'uint64');
            x=swapbytes(x);
            p=p+8;
        end;
    end

    function x=GetWord
        x=typecast(d(p:p+1),'int16');
        x=swapbytes(x);
        p=p+2;
    end

    function x=GetInt
        x=typecast(d(p:p+1),'int16');
        x=swapbytes(x);
        p=p+2;
    end
    function x=GetByte
        x=d(p);
        p=p+1;
    end

    function mprintf(handles,varargin)
    % function mprintf(handles,varargin) % copy of my utility function to
    % make ReadDM3 self-contained.
    % Write the same formatted text to multiple files.  handles is an array of
    % file handles.  The function fprintf is called multiple times, once for
    % each handle number.  Handles of 0 are ignored.
    for i=1:numel(handles)
        if handles(i)>0
            fprintf(handles(i),varargin{:});
        end;
    end
    end % mprintf
end