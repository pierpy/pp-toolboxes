function output = lab_color2html(color,text)

if iscell(text)
    tmp = char(text{end});
    text = tmp;
end
if size(color,2) == 3 & size(color,1) > 1
    color = color(end,:);
end
if isnumeric(color)
    output = ['<html><table border=0 width=400 bgcolor=rgb(',num2str(color(1)*255),',',...
        num2str(color(2)*255),',',num2str(color(3)*255),')><TR><TD>',text,'</TD></TR> </table></html>'];
else
    output = ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
end

end