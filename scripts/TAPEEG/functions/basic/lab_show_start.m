% Shows initial window of TAPEEG
%
% written by F. Hatz 2012

function lab_show_start

global VERSION

disp('')

fig = figure('renderer','zbuffer','Color',[0 0 0],'Menubar','none');
set(0,'units','pixels');
pos = get(0,'ScreenSize');
pos = [100 (pos(4)-300) 460 225];
if pos(2) < 0
    pos(2) = 100;
end
set(fig,'position',pos,'MenuBar','None','Name','TAPEEG','NumberTitle','off');
axes('units', 'pixels','position', [180 25 270 190],'Visible','off');
mylogo = lab_logo;
image(mylogo);
axis off
%uicontrol('Style', 'text','String','written by F. Hatz 2013','fontsize',8,'Position', ...
%    [280 4 200 16],'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],'HorizontalAlignment','right');
pos= get(gca,'Position');
text(pos(3),pos(4)+20,'written by F. Hatz','FontSize',8,'BackgroundColor',[0 0 0], ...
    'Color',[1 1 1],'HorizontalAlignment','right','VerticalAlignment','bottom');
text(pos(3),0,VERSION,'FontSize',8,'BackgroundColor',[0 0 0], ...
    'Color',[1 1 1],'HorizontalAlignment','right','VerticalAlignment','top');
btnw = 160;
btnh = 30;
btnx = 20;
btny = 180;
btn1 = uicontrol('style', 'pushbutton', 'units', 'pixels',...
    'position', [btnx, btny, btnw, btnh],...
    'string','Process data','fontsize',13,'Callback','close;lab_tapeeg;lab_show_start;',...
    'TooltipString','Process files or folder with eeg/meg data'); %#ok<NASGU>
btnx = 20;
btny = 140;
btn2 = uicontrol('style', 'pushbutton', 'units', 'pixels',...
    'position', [btnx, btny, btnw, btnh],...
    'string','Edit settings','fontsize',13,'Callback','close;lab_select_settings;lab_show_start;',...
    'TooltipString','Edit settings.mat (created when processing files)'); %#ok<NASGU>
btnx = 20;
btny = 100;
btn3 = uicontrol('style', 'pushbutton', 'units', 'pixels',...
    'position', [btnx, btny, btnw, btnh],...
    'string','Collect data','fontsize',13,'Callback','close;lab_select_collect;lab_show_start;',...
    'TooltipString','Collect results of frequency analysis, bad channels or eeg/meg data'); %#ok<NASGU>
btnx = 20;
btny = 60;
btn4 = uicontrol('style', 'pushbutton', 'units', 'pixels',...
    'position', [btnx, btny, btnw, btnh],...
    'string','Statistics','fontsize',13,'Callback','close;lab_select_statistics;lab_show_start;',...
    'TooltipString','Statistics'); %#ok<NASGU>
btnx = 20;
btny = 20;
btn5 = uicontrol('style', 'pushbutton', 'units', 'pixels',...
    'position', [btnx, btny, btnw, btnh],...
    'string','Plot','fontsize',13,'Callback','close;lab_select_plot;lab_show_start;',...
    'TooltipString','Plot results'); %#ok<NASGU>

% Create menus
lab_select_main(1);
lab_select_settings(1);
lab_select_collect(1);
lab_select_statistics(1);
lab_select_plot(1);
lab_select_extra(1);
lab_select_system(1);

return
