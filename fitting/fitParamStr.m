function strParams = fitParamStr(params)

params = round(params,2);

strParams = ['' ...
    '1. T: ' num2str(params(1)) ' K' newline ...
    '2. dnuGL: ' num2str(params(2)) ' 1/cm' newline ...
    '3. dnuLL:' num2str(params(3)) ' 1/cm' newline ...
    '4. Scale:' num2str(params(4)) newline ...
    '5. dnuSh:' num2str(params(5)) ' 1/cm' newline ...
    '6. Offset:' num2str(params(6)) newline ...
    '7. Slope:' num2str(params(7))];


end
