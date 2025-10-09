function tsolve_values = parse_tsolve(filename)
    fid = fopen(filename, 'r');
    tsolve_values = [];
    capture_tsolve = false;
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'Tsolve =')
            capture_tsolve = true;
        elseif capture_tsolve
            value = str2double(strtrim(line));
            if ~isnan(value)
                tsolve_values = [tsolve_values; value];
                capture_tsolve = false;
            end
        end
    end
    fclose(fid);
end

% % Example usage
% filename = 'data.txt';
% tsolve_values = parse_tsolve(filename);
% disp(tsolve_values);
