%This script is meant specifically for processing the Barium ion
%chain selectivity data. Define the variable Filepath, which should be a
%string indicating the folder containing the files to be processed, before
%running this script. 
%e.g. Filepath = 'Z:\Lab Data\Sessions\2021\2021_09\2021_09_03\Ion_Chain_Exp3_332\';

files = dir([Filepath 'Ion-Image_*']);

Repump_num = 100;
Set_num = numel(files)/Repump_num;

files_underscore_ind = find(files(1).name == '_');
A = double(imread([Filepath files(1).name]));
Image_array = nan(size(A,1),size(A,2),Set_num,Repump_num);

for h = 1:numel(files)
    A = double(imread([Filepath files(h).name]));
    exp_ind = str2num(files(h).name(files_underscore_ind(1)+1:files_underscore_ind(1)+2))+1;
    repump_ind = str2num(files(h).name(files_underscore_ind(2)+1:files_underscore_ind(2)+2))+1;
    Image_array(:,:,exp_ind,repump_ind) = A;
end

number_of_ions = nan(Set_num,Repump_num);

pngFilepath = [Filepath 'Processed_images\'];
if exist(pngFilepath) == 0
    mkdir(pngFilepath)
end

for exp_ind = 1:Set_num
    pngFilepath = [Filepath 'Processed_images\' 'Individual_PNG\'];
    for repump_ind = 1:Repump_num
        [peaks,locx,locy,C] = findpeaks2D(Image_array(:,:,exp_ind,repump_ind));
        %Frame = Image_array(:,:,exp_ind,repump_ind);
        number_of_ions(exp_ind,repump_ind) = numel(locx(peaks > 100));
        A = Image_array(:,:,exp_ind,repump_ind);
        p = pcolor(A);
        daspect([1 1 1]);
        cc = colorbar;
        xlabel('x pixel');
        ylabel('y pixel');
        ylabel(cc,'Photon count');
        set(p,'linestyle','none');
        pngfilename = ['Pcolor_image_' num2str(exp_ind-1,'%02.f') '_' num2str(repump_ind-1,'%02.f') '.png'];
        CustomSaveAsPNG;
    end
    B = sum(Image_array(:,:,exp_ind,:),4);
    pdfFilepath = [Filepath 'Processed_images\' 'Aggregated_PDF\'];
    p = pcolor(B);
    daspect([1 1 1]);
    cc = colorbar;
    xlabel('x pixel');
    ylabel('y pixel');
    ylabel(cc,'Photon count');
    set(p,'linestyle','none');
    pdffilename = ['Pcolor_image_' num2str(exp_ind-1,'%02.f') '.pdf'];
    CustomSaveAsPDF;
end