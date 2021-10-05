function h=Layout(width,height)

h = figure('PaperType','<custom>','PaperSize',[width height],'PaperPosition', [0 0 width height], 'color','w');
set(0,'DefaultAxesLineStyleOrder','-|-.', 'DefaultTextInterpreter','latex', 'DefaultAxesFontName', 'Times', 'DefaultTextFontName','Times','DefaultAxesFontSize',9,'DefaultTextFontSize',9)

end