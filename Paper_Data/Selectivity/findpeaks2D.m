function [peaks,locx,locy,C] = findpeaks2D(InputMatrix)

FX = nan(size(InputMatrix));
FY = nan(size(InputMatrix));

FX(:,1:end-1) = InputMatrix(:,2:end) - InputMatrix(:,1:end-1);
FY(1:end-1,:) = InputMatrix(2:end,:) - InputMatrix(1:end-1,:);

FX = sign(FX);
FY = sign(FY);

FXX = nan(size(FX));
FYY = nan(size(FY));

FXX(:,2:end-1) = FX(:,2:end-1) - FX(:,1:end-2);
FYY(2:end-1,:) = FY(2:end-1,:) - FY(1:end-2,:);

B = FXX + FYY;
C = zeros(size(B));
C(B == -4) = 1;

peaks = InputMatrix(C==1);
[locy,locx] = find(C==1);