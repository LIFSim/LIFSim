function ts= getTimeStamp
% Returns current time as string formatted: yyyymmdd-HHMMSS
ts = datestr(now, 'yyyymmdd-HHMMSS', 'local');
end