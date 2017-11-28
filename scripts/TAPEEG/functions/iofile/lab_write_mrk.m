% Write Cartool marker file .mrk
%
% lab_write_mrk(filename,header)
%
% written by F. Hatz 2012

function lab_write_mrk(filename,header)

skipprocessing = 0;

if isstruct(header) & isfield(header,'events') & isfield(header.events,'POS') & ~isempty(header.events.POS)
    events = header.events;
elseif isstruct(header) & isfield(header,'POS') & ~isempty(header.POS)
    events = header;
else
    skipprocessing = 1;
end

if skipprocessing == 0
    if length(events.POS) > 500
        save([filename 's'],'events');
    else
        fid=fopen(filename,'w');
        fprintf(fid,'TL02');
        for i = 1:size(events.POS,2)
            fprintf(fid,native2unicode([13 10]));
            fprintf(fid,sprintf('% 12.0f',events.POS(1,i)-1));
            fprintf(fid,'\t');
            fprintf(fid,sprintf('% 12.0f',events.DUR(1,i) + events.POS(1,i) - 1));
            fprintf(fid,'\t');
            fprintf(fid,['"' events.TYP{1,i} '"']);
        end
        fclose(fid);
    end
end