Z = peaks; surf(Z)
axis tight
set(gca,'NextPlot','replacechildren');
for j = 1:20
    surf(sin(2*pi*j/20)*Z,Z)
    F(j) = getframe;
end
movie(F,20) % Play the movie twenty times
%movie2avi(F,'myavifile.avi','Compression','Cinepak')


