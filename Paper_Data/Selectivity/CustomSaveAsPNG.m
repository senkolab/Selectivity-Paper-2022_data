h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
set(h,'Color','white');
set(h, 'InvertHardCopy', 'off')
if exist(pngFilepath) == 0
    mkdir(pngFilepath)
end
print([pngFilepath pngfilename],'-dpng', '-r0');