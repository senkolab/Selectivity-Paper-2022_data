h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
set(h,'Color','none');
set(h, 'InvertHardCopy', 'off')
if exist(pdfFilepath) == 0
    mkdir(pdfFilepath)
end
print([pdfFilepath pdffilename],'-dpdf','-r0')