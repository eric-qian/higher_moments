%% Preliminaries
clear
clc

pathData = 'data/';

% FRED mnemonics
mnems  = {'GDPC1', ...  Real GDP
    'GDPPOT', ...       Potential GDP
    'UNRATE', ...       Unemployment rate
    'NROU', ...         Natural rate
    'GDPDEF', ...       GDP deflator
    'TB3MS', ...        3M
    'WTISPLC'};  %      WTI spot crude


% Final estimation sample will start 1959:Q4. Fetch additional quarters for padding.
tStart = datetime(1956, 1, 1);  
tEnd   = datetime(2019, 12, 31);   


%% Fetch raw data

n      = length(mnems);
Raw    = timetable('RowTimes', tStart:calquarters(1):tEnd);
T      = height(Raw);


% Fetch from FRED
url = 'https://fred.stlouisfed.org/';
c   = fred(url);


for jSeries = 1:n
    Resj  = fetch(c, mnems{jSeries}, tStart, tEnd);
    ttj   = timetable(datetime(Resj.Data(:, 1), 'ConvertFrom', 'datenum'),...
        Resj.Data(:, 2));
    ttjq  = retime(ttj, 'quarterly', 'mean') ;  % Take quarterly average
    assert(all(ttjq.Time == Raw.Time))          % Make sure time indices match
    Raw.(mnems{jSeries}) = ttjq.Var1;           % Add to raw data
            
    
    % Plot data, verify no issues with aggregation
    figure()
    plot(ttj.Time, ttj.Var1)
    hold on
    plot(ttjq.Time, ttjq.Var1)
    legend({'original', 'aggregated'})
    title(mnems{jSeries})
end


%% Add transformations
% See footnote 20 of Gourieroux, Monfort, Renne (2017) for details on 
% real economic activity series and oil price transformation.

df      = Raw;
df.ygap = log(df.GDPC1) - log(df.GDPPOT);                                % Output gap
df.ugap = df.UNRATE - df.NROU;                                           % Unemployment gap
%df.pi   = 100*[NaN; (df.GDPDEF(2:end)./ df.GDPDEF(1:end-1)).^4 - 1];     % CAGR
df.pi   = 400*[NaN; log(df.GDPDEF(2:end)) - log(df.GDPDEF(1:end-1))];     
df.doil = 100*[NaN; log(df.WTISPLC(2:end)) - log(df.WTISPLC(1:end-1))];  % Growth in oil price   
df.R    = df.TB3MS;
df      = df(2:end, :);                                                  % Drop first period

writetimetable(df, [pathData 'data.xls'])
