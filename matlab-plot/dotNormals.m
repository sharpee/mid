function result = dotNormals(modelV1, modelV2, MIDVectorPrefix, nParts, ...
                             Nh, Nv, nlags)
%
% Calculate the dot product of the normal vectors of two planes.
% must be run from directory with .dat files in it
%
% USAGE: result = dotNormals(modelV1, modelV2, MIDVectorPrefix, nParts,
%                            Nh, Nv, nlags)
%
% modelV1 : Model Vector 1 file name (e.g. 'mv_model_v1_1110_1.dat'). 
%           A PNG image is generated in the output folder.
%
% modelV2 : Model Vector 2 file name (e.g. 'mv_model_v2_1110_1.dat'). 
%           A PNG image is generated in the output folder.
%
% MIDVectorPrefix : The prefix string (e.g. 'V1model-1D-n2') from which
%                   the file names of the individual data vector files can be
%                   generated as:
%                     sprintf('%s-v%u-p%u.dat',MIDVectorPrefix,vector,part);
%                   These files are averaged together, and a png is generated.
%
% nParts: Number of parts in the analysis ("number of parts" in
%         params.xml Typically 4).
%
% Nh: Horizontal size ("sta width"/"x downsample" in params.xml. Typically 16)
%
% Nv: Vertical size. ("sta height"/"y downsample" in params.xml. Typically 16)
%
% nlags: number of frames. ("sta duration" in params.xml. Typically 3)
%

mkdir './output';  
%load model vectors
fp = fopen(modelV1);
v1 = fread(fp, inf, 'double');
fclose(fp);

parsed = regexp(modelV1, '\.', 'split');
imageName = ['./output/' char(parsed(1)) '.png'];
saveVector(v1, imageName, Nh, Nv, nlags);

fp = fopen(modelV2);
v2 = fread(fp, inf, 'double');
fclose(fp);

parsed = regexp(modelV2, '\.', 'split');
imageName = ['./output/' char(parsed(1)) '.png'];
saveVector(v2, imageName, Nh, Nv, nlags);

fileSize = size(v1);
fileSize = fileSize(1);

u1 = zeros(fileSize,1);
u2 = zeros(fileSize,1);
vector = 1;
%load data vectors and average
for part=1:nParts
    fileName = sprintf('%s-v%u-p%u.dat', MIDVectorPrefix, vector, part);
    fp= fopen(fileName);
    u1d = fread(fp,inf,'double');
    if(dot(u1d, u1) < 0)
        u1d = -u1d;
    end
    u1 = u1 + u1d;
    fclose(fp);     
end
u1 = u1/norm(u1);
imageName = sprintf('./output/%s-v%u.png', MIDVectorPrefix, vector);
saveVector(u1, imageName, Nh, Nv, nlags);


vector = 2;
for part=1:nParts
    fileName = sprintf('%s-v%u-p%u.dat', MIDVectorPrefix, vector, part);
    fp= fopen(fileName);
    u2d = fread(fp, inf, 'double');
    if(dot(u2d, u2) < 0)
       u2d = -u2d; 
    end
    u2 = u2 + u2d;
    fclose(fp);
end
u2 = u2/norm(u2);
imageName = sprintf('./output/%s-v%u.png', MIDVectorPrefix, vector);
saveVector(u2, imageName, Nh, Nv, nlags);

%calculate dot product of normal vectors (1 for matching, 0 for orthogonal)
numerator = det([dot(v1,u1) dot(v1, u2); dot(v2,u1) dot(v2, u2)]);
denominator = sqrt(det([dot(u1,u1) dot(u1, u2); dot(u1, u2) dot(u2, u2)]));

result = abs(numerator/denominator);

