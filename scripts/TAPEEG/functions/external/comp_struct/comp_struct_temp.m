function [df, match, er1, er2, erc, erv] = comp_struct(s1,s2,prt,pse,tol,n1,n2,wbf)
% check two structures for differances - i.e. see if strucutre s1 == structure s2
% function [match, er1, er2, erc, erv] = comp_struct(s1,s2,prt,pse,tol,n1,n2,wbf)
%
% inputs  8 - 7 optional
% s1      structure one                              class structure
% s2      structure two                              class structure - optional
% prt     print test results (0 / 1 / 2 / 3)         class integer - optional
% pse     pause flag (0 / 1 / 2)                     class integer - optional
% tol     tol default tolerance (real numbers)       class integer - optional
% n1      first structure name (variable name)       class char - optional
% n2      second structure name (variable name)      class char - optional
% wbf     waitbar flag (0 / 1) default is 1          class integer - optional
%
% outputs 6 - 6 optional
% df      mis-matched fields with contents           class cell - optional
% match   matching fields                            class cell - optional
% er1     non-matching feilds for structure one      class cell - optional
% er2     non-matching feilds for structure two      class cell - optional
% erc     common error (both structures listed)      class cell - optional
% erv     common error values (both structures)      class cell - optional
%
% prt:
%	0 --> no print
%	1 --> print major erros
%	2 --> print all errors
%	3 --> print errors and matches
% pse:
%	1 --> pause for major erros
%	2 --> pause for all errors
%
% example:	[match, er1, er2] = comp_struct(data1,data2,1,1,1e-6,'data1','data2')
% michael arant - may 27, 2013
%
% updated - aug 22, 2013
%
% hint:
% passing just one structure causes the program to copy the structure
% and compare the two.  This is an easy way to list the structure

if nargin < 1; help comp_struct; error('I / O error'); end
if nargin < 2; s2 = s1; prt = 3; end
if nargin < 3 || isempty(prt); prt = 1; end
if nargin < 4 || isempty(pse); pse = 0; elseif pse ~= 1 && prt == 0; pse = 0; end
if nargin < 5 || isempty(tol); tol = 1e-6; end
if nargin < 6 || isempty(s1); n1 = 's1'; end
if nargin < 7 || isempty(s2); n2 = 's2'; end
if nargin < 8 || isempty(wbf); wbf = 1; end
if pse > prt, pse = prt; end

% solve
[match, er1, er2, erc, erv] = comp_struct_loop(s1,s2,prt,pse,tol,n1,n2,wbf);
keyboard; stop


%% recursive loop
function [match, er1, er2, erc, erv] = comp_struct_loop(s1,s2,prt,pse,tol,n1,n2,wbf)

% init outputs
match = {}; er1 = {}; er2 = {}; erc = {}; erv = cell(0,2);

% test to see if both are structures
if isstruct(s1) && isstruct(s2)
	% both structures - get the field names for each structure
	fn1 = fieldnames(s1);
	fn2 = fieldnames(s2);
	% missing fields? get the common fields
	temp1 = ismember(fn1,fn2);
	temp2 = ismember(fn2,fn1);
	% missing fileds in set 1
	for ii = find(~temp2)'
		er1{end+1} = sprintf('%s is missing field %s',n1,fn2{ii});
		if prt; fprintf('%s\n',er1{end}); end; if pse; pause; end
		erc(end+1) = er1(end);
		erv(end+1,2) = {s2(1).(fn2{ii})};
	end
	% missing fields in set 2
	for ii = find(~temp1)'
		er2{end+1} = sprintf('%s is missing field %s',n2,fn1{ii});
		if prt; fprintf('%s\n',er2{end}); end; if pse; pause; end
		erc(end+1) = er2(end);
		erv(end+1,1) = {s1(1).(fn1{ii})};
	end
	% index sizes match?  i.e. do both structures have the same # of indexes?
	inda = numel(s1); indb = numel(s2); inder = inda-indb;
	if inder < 0
		% struct 1 is smaller
		for ii = inda+1:indb
			er1{end+1} = sprintf('%s(%g) is missing',n1,ii);
			if prt; fprintf('%s\n',er1{end}); end; if pse; pause; end
			erc(end+1) = er1(end);
			erv(end+1,2) = {s2(ii).(fn1{ii})};
		end
	elseif inder > 0
		% index 2 is smaller
		for ii = indb+1:inda
			er2{end+1} = sprintf('%s(%g) is missing',n2,ii);
			if prt; fprintf('%s\n',er2{end}); end; if pse; pause; end
			erc(end+1) = er2(end);
			erv(end+1,1) = {s1(ii).(fn1{ii})};
		end
	end
	% get common fields
	fn = fn1(temp1); fnn = numel(fn); 
	% loop through structure 1 and match to structure 2
	ind = min([inda indb]); cnt = 0; 
	if wbf; wb = waitbar(0,'Comparing ....'); end
	for ii = 1:ind
		% loop each index
		for jj = 1:fnn
			% loop common field names
			if wbf; cnt = cnt + 1; waitbar(cnt/(ind*fnn),wb); drawnow; end
			% add index and field name to the structure name
			n1p = sprintf('%s(%g).%s',n1,ii,fn{jj});
			n2p = sprintf('%s(%g).%s',n2,ii,fn{jj});
			% recurse - run the program again on the sub-set of the structure
			[m e1 e2 ec ev] = comp_struct_loop(s1(ii).(fn{jj}),s2(ii).(fn{jj}),prt,pse, ...
				tol,n1p,n2p,wbf);
			% add the sub-set (field name) results to the total results
			match = [match m']; 
			if ~isempty(e1) || ~isempty(e2)
				er1 = [er1 e1']; er2 = [er2 e2'];
				if numel(erc) > 1 && size(erc,1) == 1; erc = erc'; end
				if numel(ec) > 1 && size(ec,1) == 1; ec = ec'; end
				erc = [erc; ec]; 
				erv = [erv; ev];
			end
		end
	end
	if wbf;	close(wb); end
else
	% both are non-structures - compare
	% get the varable class and test
	c1 = class(s1); c2 = class(s2);
	if strcmp(c1,c2);
		% both are the same class
		if isequal(s1,s2)
			% results are equal
			match{end+1} = sprintf('%s and %s match',n1,n2);
			if prt == 3; fprintf('%s\n',match{end}); end
		else
			% same class but not equal
			% calculate error if type is single or double
			% test for function type match if function handle
			switch c1
				case {'single', 'double'}, 
					if numel(s1) ~= numel(s2); er = 1; else; er = abs(s1 - s2); end
				case {'function_handle'},
					s1f = functions(s1); s2f = functions(s2);
					if strcmp(s1f.function,s2f.function)
						% same function with different values - record deviation and exit
						er = 0;
						er1{end+1} = sprintf('%s and %s are both %s but have different values', ...
							n1,n2,char(s1));
						er2{end+1} = er1{end};
						if prt > 1; fprintf('%s\n',er1{end}); end;
						if pse > 1; pause; end
						erc(end+1) = er1(end);
						erv(end+1,1) = {s1}; erv(end,2) = {s2};
					else
						er = 1;
					end
				otherwise, er = 1;
			end
			% test error - error will be 0 (no error) or 1 (error) for all
			% classes except double and single.  double and single are the 
			% actual error which is tested against the tolerance
			% this was done for cases where structures are run on different 
			% platforms and numerical precision errors are observed
			if er > tol
				% sets do not match
				er1{end+1} = sprintf('%s and %s do not match',n1,n2);
				er2{end+1} = sprintf('%s and %s do not match',n1,n2);
				if prt > 1; fprintf('%s\n',er1{end}); end;
				if pse > 1; pause; end
				erv(end+1,1) = {s1}; erv(end,2) = {s2}; erc(end+1) = er1(end);
			else
				% sets are a tolerance match
				match{end+1} = sprintf('%s and %s are tolerance match',n1,n2);
				if prt > 2; fprintf('%s\n',match{end}); end
			end
		end
	else
		% fields are different classes
		er1{end+1} = sprintf('%s is class %s, %s is class %s',n1,c1,n2,c2);
		er2{end+1} = sprintf('%s is class %s, %s is class %s',n1,c1,n2,c2);
		if prt; fprintf('%s\n',er1{end}); end
		if pse; pause; end
		erc(end+1) = er1(end);
		erv(end+1,1) = {c1}; erv(end,2) = {c2};
	end
end

% transpose outputs
match = match'; er1 = er1'; er2 = er2';


%% old code for historical purposes....
% % common error list?
% if nargout > 3
% 	% combine errors into single list
% 	erc = cell(0); erv = cell(0,2); 
% 	% compare errors
% 	temp1 = ismember(er1,er2); temp2 = ismember(er2,er1);
% 	% matching errors
% 	idx1 = find(temp1);
% 	% unique errors
% 	inx1 = find(~temp1); inx2 = find(~temp2);
% 	% record unique errors
% 	if ~isempty(inx1); erc = [erc; er1(inx1)]; end 
% 	if ~isempty(inx2); erc = [erc; er2(inx2)]; end 
% 	% recod matching errors
% 	if ~isempty(idx1); erc = [erc; er1(idx1)]; end 
% end

% % list error contents?
% if nargout > 4
% 	% display contents from erc - data that did nto match
% 	erv = cell(numel(erc),2);
% 	for ii = 1:numel(erc)
% 		% get the error string
% 		temp = erc{ii};
% 		% lookfor structure listed
% 		n1temp = strfind(temp,n1);
% 		n2temp = strfind(temp,n2);
% 		% collect data
% 		if isempty(n1temp)
% 			junk = findstr(' is missing field ',temp);
% 			erv(ii,1) = {eval([n1 temp(3:junk-1) '.'  temp(junk+18:end)])};
% 		elseif isempty(n2temp)
% 			junk = findstr(' is missing field ',temp);
% 			erv(ii,2) = {eval([n2 temp(3:junk-1) '.'  temp(junk+18:end)])};
% 		else
% 			% both fields exist
% 			% class error?
% 			junk = strfind(temp,' is class ');
% 			if ~isempty(junk)
% 				trash = strfind(temp,', '); 
% 				erv(ii,1) = {temp(junk(1)+4:trash-1)};
% 				erv(ii,2) = {temp(junk(2)+4:end)};
% 			end
% 			% differnet content?
% 			junk = strfind(temp,' do not match');
% 			if ~isempty(junk)
% 				trash = strfind(temp,' and ');
% 				erv(ii,1) = {eval(temp(1:trash-1))};
% 				erv(ii,2) = {eval(temp(trash+5:junk-1))};
% 			end
% 		end
% 		
% 		
% 		
% % 		% get field definitions
% % 		temp = strfind(erc{ii},' ');
% % 		% parse fields
% % 		f1 = erc{ii}(1:temp(1)-1);
% % 		f2 = erc{ii}(temp(2)+1:temp(3)-1);
% % 		% get contents
% % 		f1 = regexprep(f1,n1,'s1','once');
% % 		f2 = regexprep(f2,n2,'s2','once');
% % 		% evaluate
% % % 		erv(ii,1) = {eval(f1)};
% % % 		erv(ii,2) = {eval(f2)};
% % 		erv(ii,1) = {eval(f1)};
% % 		erv(ii,2) = {eval(f2)};
% 	end
% end
