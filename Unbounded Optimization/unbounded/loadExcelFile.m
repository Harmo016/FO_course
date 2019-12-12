function [dates, ric, prices] = loadExcelFile(fileName, sheetName)

[data,txt,raw] = xlsread(fileName, sheetName, '', 'basic');
ric = raw(1,2:end);

dates = cell2mat(raw(3:end,1));
prices = cell2mat(raw(3:end,2:end));
dates = dates + datenum(1899,12,30);

% Remove dates which include NaN values for any asset
ind = ~any(isnan(prices),2);
dates= dates(ind);
prices = prices(ind,:);

% Change order (oldest value first)

dates = dates(end:-1:1);
prices = prices(end:-1:1,:);
