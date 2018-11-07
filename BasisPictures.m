%% Basis Pictures

clearvars
InitializeDisplay
set(0,'defaultLineLineWidth',6) %thick lines

n1plot = 1e4;
x1plot = -1:2/(n1plot-1):1;
n2plot = 500;
x2plot = -1:2/(n2plot-1):1;
[xx, yy] = meshgrid(x2plot);

s_max = 4;
d = 2;
twoD = [1 1; 1 2; 1 3; 1 4; 2 2; 2 3];
nTwoD = size(twoD,1);

basisHandles = {@legendreBasis,@chebyshevBasis};
basisNames = ["Legendre","Chebyshev"];
nBasis = length(basisNames);
for nB = 1:nBasis
   % One dimensional basis plots
   horizOff = -0.7;
   vertOff = 1.2;
   for s = 0:s_max
      y = basisHandles{nB}(s,x1plot);
      figure
      plot(x1plot,y,'-')
      axis([-1 1 -1.2 1.2])
      set(gca,'Visible','off')
      text(horizOff,vertOff,['\boldmath\(j\,\) \unboldmath \(= (' int2str(s) ', 0, \ldots, 0)\)'], ...
         'interpreter','latex','fontsize',48)
       eval(['print(''' char(basisNames(nB)) '_Degree_' int2str(s) '.png'',' ...
          '''-dpng'')'])
      eval(['print(''' char(basisNames(nB)) '_Degree_' int2str(s) '.eps'',' ...
         '''-depsc'', ''-opengl'')'])
   end
   xOff = -1;
   yOff = 0.5;
   zOff = 2.2;
   for j = 1:nTwoD
      z = basisHandles{nB}(twoD(j,1),xx).*basisHandles{nB}(twoD(j,2),yy);
      figure
      surf(xx,yy,z); shading interp
      axis([-1 1 -1 1 -1.2 1.2])
      set(gca,'Visible','off')
      text(xOff,yOff,zOff,['\boldmath\(j\,\) \unboldmath \(= (' int2str(twoD(j,1)), ...
         ',' int2str(twoD(j,2)) ', 0, \ldots, 0)\)'], ...
         'interpreter','latex','fontsize',48)
       eval(['print -dpng ' char(basisNames(nB)) '_Degree_' int2str(twoD(j,1)) ...
          '_' int2str(twoD(j,2)) '.png'])
      eval(['print -depsc ' char(basisNames(nB)) '_Degree_' int2str(twoD(j,1)) ...
         '_' int2str(twoD(j,2)) '.eps'])
   end
end


