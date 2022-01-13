function threshold = OutlierUpper(InputMatrix,Percentile,Multiplier)

threshold = median(InputMatrix(:)) + Multiplier*(prctile(InputMatrix(:),Percentile) - median(InputMatrix(:)));