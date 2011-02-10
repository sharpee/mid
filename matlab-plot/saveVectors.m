function saveVectors(Nh, Nv, nlags)
%
% Generates figures of all output data.
% The active directory must contain the output data.
%
% Usage: saveVectors(Nh, Nv, nlags)
%

mkdir('./output');
list = dir('*.dat');
for i= 1:length(list)
    parsed = regexp(list(i).name, '\.', 'split');
    if(length(parsed)==2)
        fileName = list(i).name;
        outputName = ['./output/' char(parsed(1)) '.png'];
        fp=fopen(fileName, 'rb');

        Nn = Nh*Nv*nlags;
        if ~(fp==-1)
            vect=fread(fp,Nn,'double');
            minv = min(vect);
            maxv = max(vect);
            vect = reshape(vect, Nh, Nv*3);
            fig1 = figure;
            set(fig1, 'OuterPosition', [100 500 1000 400]); %left bottom width height
            for lag=1:nlags
                subplot(1,3,lag);
                imagesc(vect(1:Nh, 1+Nv*(lag-1):Nv*lag), [minv, maxv]);     
                colormap('hot');           
                colorbar;
                axis square;
       
            end
            saveas(fig1, outputName, 'png');
            close(fig1);
        end
    end
end

