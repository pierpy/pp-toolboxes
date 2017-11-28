function lab_set_X_Angle(labels,angle)
   text_h = findobj(gca,'Type','text');
   if length(text_h) < length(labels)
       return
   end
   for cnt = 1:length(labels)
       if angle == 0
           set(text_h(length(text_h)-cnt+1),'Rotation',0, ...
               'String',labels{cnt},'HorizontalAlignment','center');
       elseif angle > 0
           set(text_h(length(text_h)-cnt+1),'Rotation',angle, ...
               'String',labels{cnt},'HorizontalAlignment','right');
       elseif angle < 0
           set(text_h(length(text_h)-cnt+1),'Rotation',angle, ...
               'String',labels{cnt},'HorizontalAlignment','left');
       end
   end
end