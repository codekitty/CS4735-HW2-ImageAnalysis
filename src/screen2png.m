% SCREEN2JPEG Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
function screen2png(h,filename)

% create directory if needed
[pathstr, ~, ~] = fileparts(filename);
if ~exist(pathstr, 'dir')
   mkdir(pathstr); 
end

set(gcf,'Units','pixels');

pos = get(h,'Position');
newpos = pos/100;

% sets the position\size of the current graphical object before printing
set(h,'PaperUnits','inches',...
     'PaperPosition',newpos)

% print
print('-dpng', filename, '-r100');
drawnow;
end
