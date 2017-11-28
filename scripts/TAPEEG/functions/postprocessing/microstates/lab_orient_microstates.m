function Template = lab_orient_microstates(Template,header)
    
if isstruct(Template) & isfield(Template,'Template') & isfield(Template,'header')
    for i = 1:length(Template)
        T = Template(i).Template;
        IdxL = Template(i).header.locs.y < 0;
        IdxH = Template(i).header.locs.y > 0;
        for j = 1:size(T,2)
            if mean(T(IdxL,j)) > mean(T(IdxH,j))
                T(:,j) = -T(:,j);
            end
        end
        Template(i).Template = T;
    end
elseif exist('header','var') & isfield(header,'locs') & ~isempty(header.locs)
    IdxL = header.locs.y < 0;
    IdxH = header.locs.y > 0;
    for j = 1:size(Template,2)
        if mean(Template(IdxL,j)) > mean(Template(IdxH,j))
            Template(:,j) = -Template(:,j);
        end
    end
else
    return
end
