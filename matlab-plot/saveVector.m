function saveVector(vector, fileName, Nh, Nv, nlags)
%
% Save vector to disk as a PNG image.
%
% Usage: saveVector(vector, fileName, Nh, Nv, nlags)
%

Nn = Nh*Nv*nlags;           
minv = min(vector);
maxv = max(vector);
vect = reshape(vector, Nh, Nv*3);
fig1 = figure;
set(fig1, 'OuterPosition', [100 500 1000 400]); %left bottom width height
for lag=1:nlags
    subplot(1,3,lag);
    imagesc(vect(1:Nh, 1+Nv*(lag-1):Nv*lag), [minv, maxv]);     
    colormap('hot');           
    colorbar;
    axis square;       
end
saveas(fig1, fileName, 'png');
close(fig1);

