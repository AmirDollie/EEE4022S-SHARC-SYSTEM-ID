function file = writeStudyCsv(run, name, header, rows)
%WRITESTUDYCSV  Write one deliverable table as CSV (readable in Excel and by pgfplotstable in LaTeX).
%   file = writeStudyCsv(run, 'gamma_tolerance', {'sensor', 'eps', 'dSys'}, rows)
%
%   rows  cell array, one row per line; numbers are written with %.6g, NaN and [] as empty,
%         logicals as 0/1, text as is. Text containing a comma, quote or line break is quoted
%         ("...") with embedded quotes doubled (RFC 4180), so it stays in one column.
%   Writes <Results>/<ID>/<ID>_<name>.csv (no timestamp: the latest run is the one the thesis uses).
%   finishStudy lists every table and figure written during the run in its manifest.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
file = fullfile(run.dir, sprintf('%s_%s.csv', run.id, name));
fid = fopen(file, 'w');
if fid < 0, error('writeStudyCsv:open', 'Cannot write %s (open in Excel?)', file); end
fprintf(fid, '%s\n', strjoin(cellfun(@csvText, header, 'UniformOutput', false), ','));
for i = 1:size(rows, 1)
    cells = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if ischar(v), cells{j} = csvText(v);
        elseif islogical(v), cells{j} = sprintf('%d', v);
        elseif isempty(v) || (isnumeric(v) && isnan(v)), cells{j} = '';
        else, cells{j} = sprintf('%.6g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(cells, ','));
end
fclose(fid);
fprintf('  table saved: %s\n', file);
end

function t = csvText(t)
if any(t == ',' | t == '"' | t == sprintf('\n') | t == sprintf('\r'))
    t = ['"' strrep(t, '"', '""') '"'];
end
end
