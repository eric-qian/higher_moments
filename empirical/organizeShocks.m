%% organizeShocks.m  Input results .mat file and output Excel spreadsheet

initSetting = 'GlobalSearch';  % Initialization setting


outFile = ['figures2019/initSetting=' initSetting '_shocks.xls'];
delete(outFile)  % Delete file if exists (prevent duplicate sheets)

load(['figures2019/initSetting=' initSetting '_Results'])


for jSpec = 1:nSpec
    shocks = Spec(jSpec).shocks;
    Time   = Spec(jSpec).dfEst.Properties.RowTimes;
    p      = Spec(jSpec).p;
    
    % Make table
    df        = array2table(Time(p+1:end), 'VariableNames', {'Time'});
    df.shock1 = shocks(:, 1);
    df.shock2 = shocks(:, 2);
    df.shock3 = shocks(:, 3);
    
    % Reformat sheet name
    sheetName = strrep(Spec(jSpec).nameDate, ':', '');
    sheetName = strrep(sheetName, ', ', ' ');
    sheetName = strrep(sheetName, ' (', ' ');
    sheetName = strrep(sheetName, ')' , '' );
    

   writetable(df, outFile, 'Sheet', sheetName)
end


