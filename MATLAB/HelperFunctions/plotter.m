function s = plotter(h, X, Y)

s = surf(X,Y,h,'FaceAlpha',0.5);
s.EdgeColor = 'none';
end