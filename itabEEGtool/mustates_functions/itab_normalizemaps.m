function maps = itab_normalizemaps(maps, modelmaps_normalization_type, gfp_method)
    if  strcmp(modelmaps_normalization_type,'gfp1')
        for i = 1:size(maps,1)
            normalization_factor = compute_gfp(maps(i,:), gfp_method);
            maps(i,:) = maps(i,:)./ normalization_factor;
        end
    elseif strcmp(modelmaps_normalization_type,'vector_norm_1')
        for i = 1:size(maps,1)
            normalization_factor = norm(maps(i,:),2);
            maps(i,:) = maps(i,:)./normalization_factor;
        end
    end