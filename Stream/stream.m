function h = stream(d, m, list) % perform a series of tasks on a data table

% Inputs:
%   d = data table: table
%   m = metadata: structure
%		list = list of tasks to execute: table
%	Output:
%		h = handles of functions in this file
% Description:
%   Rows in the data table represent observations and columns represent
%   variables. Data are analysed in two sequences. First, the tasks listed
%   in m.tasks are performed on the table in sequence. Second, for each
%   task, data in the table are separated into groups and each group is
%   analysed in turn. Grouping variables are contained in m.task.group.
%   Plotting tasks differ in that there are two grouping variables:
%   m.plot.group indicates variables that change between axes and m.plot.line
%   indicates variables that change between lines on an axes.
	
	%	Return handles for task functions in this file
	if nargout > 0 % handles are required
		h = localfunctions; % handles of functions in this file: fs x 1, cell
		return % nothing else to do
	end

	%	Set the name of the dependent variable
	if ~ isfield(m, "z") && isfield(m, "y") % m.z doesn't exist: give it a value
		m.z = m.y; % % several functions use m.z as the name of the dependent var.
	end
	
	% Perform tasks
	task = list.task; % task sequence: ts x 1, string
	for i = 1: length(task) % loop over tasks
		taskC = task(i); % current task
		m.task = taskC; % store task name
		if i > 1 && isempty(d) % data table is empty for tasks after first
			error("Data table input to task ""%s"" is empty", taskC);
		end
		[d, m] = anaGroups(d, m, list(i, :)); % execute task
	end

function [d, m] = anaGroups(d, m, list) % analyse data group by group

% Input:
%   list = task name, handle, and grouping variables: 1 x 3, table

	% Find grouping variable(s) and function to execute
	task = list.task; % task name: string
	fun = list.hand{:}; % handle of function to execute
	var = list.group{:}; % grouping variables: string
	m.(task).fun = fun; % define at least one field for m.(task),
		% as assumed by following code *** fix ***

	% Set independent and dependent variables, x, y, and z
	for v = "xyz" % loop over variable names
		if isfield(m.(task), v) % m.task.v is defined, so use it to set m.v
			if isfield(m, v) % m.v is defined
				prev.(v) = m.(v); % save m.v for later restoration
			else % m.v is undefined
				prev.(v) = ""; % flag m.v as undefined
			end
			m.(v) = m.(task).(v); % set m.v
		end
	end

	% Split d into groups and apply function to each group
	if isempty(var) % no grouping
		g = ones(size(d, 1), 1); % all rows
		gs = 1; % one group
	else
		%	[g, dU] = findgroups(d(:, var)); % find groups specified by var
		[dU, ~, g] = unique(d(:, var), "rows"); %	find groups specified by var
		gs = size(dU, 1); % number of groups
	end
	m.groups = gs; % store number of groups
	dOut = cell(gs, 1); % allocate storage
	for j = 1: gs % loop over groups
		m.group = j; % current group number		
		%	dIn = d(g == j, :); % current group: takes too much RAM for big	files
		%	[dOut{j}, m] = fun(dIn, m); % analyse group
		[dOut{j}, m] = fun(d(g == j, :), m); % analyse group
	end
	
	% Finalise
	d = vertcat(dOut{:}); % combine data tables over groups
	for v = "xyz" % loop over independent and dependent variables
		if exist("prev", "var") && isfield(prev, v) % previous values were changed
			if isempty(prev.(v)) % m.v didn't previously exist
				m = rmfield(m, v); % remove this variable
			else % m.v previously existed
				m.(v) = prev.(v); % restore the previous value
			end
		end
	end

function [d, m] = add(d, m)

% Add data to existing axes
%
% Inputs:
%   m.add.arg (optional) = arguments for plotting function: cell
%   m.add.funFun (optional) = plotting function: function handle
%   m.line (optional) = line grouping variable(s): char or cell
%   m.x = independent variable: char
%   m.y = dependent variable: char

	% Initialise plotting function handle
	if isfield(m.add, 'funFun') % set plotting function
		fun = m.add.funFun; % user-provided function
	else
		if isfield(m.add, 'arg') % plot arguments
			arg = m.add.arg; % user-supplied arguments
		else
			arg = {}; % default
		end
		fun =  @(d, m)plot(d.(m.x), d.(m.y), arg{:}); % default
	end

	%	Initialise line grouping variables
	if isfield(m.add, 'line') % m.add.line is specified
		line = m.add.line; % names of line grouping variables
	elseif isfield(m.plot, 'line') % m.plot.line is specified
		line = m.plot.line; % use same line grouping variables as for plot task
	else % nothing specified
		line = {}; % default
	end
	line = string(line); % string array

	% Loop over lines
	hC = m.handle.axes(m.group); % axes handle
	axes(hC); % make axes current
	aC = d; % current data table
	p = group(aC, line); % group into lines
	ps = length(p); % number of lines
	set(hC, 'nextPlot', 'add'); % add lines, don't replace
	for i = 1: ps % loop over lines
		pC = p{i}; % current data to plot
		hP = fun(pC, m); % plot: pre-R2014b
		labelLine(pC, hP, line); % label line by adding display name
	end
	
function [d, m] = doAnova(d, m)

% Perform and display analyses of variance on the model(s)
%
% Input:
%   m.model = model(s) produced by task fitlm
% Output:
%   ANOVA tables, one for each model

	% Loop over models
	for i = 1: length(m.model)
		modelC = m.model{i}; % current model
		dC = modelC.Variables; % data used to obtain model
		name = modelC.PredictorNames; % predictor variables
		
		% Ensure that all predictor variables are multi-valued
		for j = 1: length(name) % loop over predictor variables
			nameC = name{j}; % current name
			val = unique(dC.(nameC)); % unique values
			if size(val, 1) == 1
				error('Predictor variable %s is single-valued', nameC);
			end
		end

		% Perform ANOVA
		disp(anova(modelC));
		
	end

function [d, m] = doExport(d, m) % export graphics object to a file

	% Set outputs
	if isfield(m.export, 'obj') % set type of graphics object to export
		t = m.export.obj; % user-specified: axes, figure, ...
	else
		t = 'figure'; % default
	end
	if isfield(m.export, 'name') % set name and type of export file
		n = m.export.name; % user-specified name
	else
		n = 'Export.eps'; % default
	end
	if isfield(m.export, 'arg') % exportgraphics arguments
		arg = m.export.arg; % user-supplied arguments
	else
		arg = {}; % default
	end

	% Export
	switch t % obtain graphics object handle
		case 'axes', h = gca; % handle of current axes
		case 'figure', h = gcf; % handle of current figure
	end
	exportgraphics(h, n, arg{:}); % export

function [d, m] = doFitglm(d, m) % fit generalised linear model

% Inputs:
%   m.fitglm.arg (optional) = list of arguments to fitlm: cell
%   m.group (internal) = group number: double
% Output:
%   m.model{m.group} = model generated by fitglm: cell member

	if isfield(m.fitglm, 'arg') % m.fitglm.arg specified
		arg = m.fitglm.arg; % user-specified arguments
	else
		arg = {}; % default
	end
	m.model{m.group} = fitglm(d, arg{:});

function [d, m] = doFitlm(d, m) % fit linear model

% Inputs:
%   m.fitlm.arg (optional) = list of arguments to fitlm: cell
%   m.group (internal) = group number: double
% Output:
%   m.model = models generated by fitlm, one for each group: cell

	% Initialise
	if isfield(m.fitlm, 'arg') % m.fitlm.arg specified
		arg = m.fitlm.arg; % user-specified arguments
	else
		arg = {}; % default
	end
	if isfield(m.fitlm, 'weight') % regression weights are specified
		w = m.fitlm.weight; % name of weight variable
		w = {'weights', d.(w)}; % arguments for fitlm
	else % no regression weights specified
		w = {}; % default
	end
		
	% Fit
	m.model{m.group} = fitlm(d, arg{:}, w{:}); % call fitlm

function [d, m] = doFitnlm(d, m) % fit nonlinear model

% Inputs:
%   m.fitnlm.arg (optional) = list of arguments to fitnlm: cell
%   m.fitnlm.funFun = model function: function handle
%   m.fitnml.beta0 = starting values for model parameters: double
%   m.group (internal) = group number: double
% Output:
%   m.model = models generated by fitnlm, one for each group: cell

	% Initialise
	if isfield(m.fitnlm, 'arg') % m.fitnlm.arg specified
		arg = m.fitnlm.arg; % user-specified arguments
	else
		arg = {}; % default
	end
	
	% Fit
	m.model{m.group} = fitnlm(d, m.fitnlm.funFun, m.fitnlm.beta0, arg{:});

function [d, m] = limit(d, m) % accept variable values within specified limits

%	Purpose:
%		Select rows in which specified variables fall within specified limits
% Input:
%		Variables are specified by fields in the metadata struture, m.limit.
%		Allowable values are specified as a closed range.

	% Initialise
	var = d.Properties.VariableNames; % names of variables in d
	field = fieldnames(m.limit); % names of fields in m.limit
	name = intersect(var, field); % names in both var and field
	names = length(name); % number of names
	rows = size(d, 1); % number of rows in d
	row = zeros(rows, names); % rows to select

	% For each name, find rows to select
	for i = 1: names % loop over names
		nameC = name{i}; % current name
		valV = d.(nameC); % values of variable
		valF = m.limit.(nameC); % value of field
		classC = class(valV); % class of variable
		switch classC % check that the variable is double
			case 'double'
				rowC = valV >= valF(1) & valV <= valF(2);
			otherwise
				warning('Row selection: class %s not supported', classC);
		end
		row(:, i) = rowC; % store result for this name
	end

	% Select rows
	row = all(row, 2); % find rows that satisfy all tests
	d = d(row, :); % select rows
		
function [d, m] = doLoad(~, m) % load data table from file

	if isfield(m, 'folder') % m.folder exists
		f = m.folder + filesep; % folder from which to load
	else
		f = ""; % default
	end
	if isfield(m.load, 'name')
		n = f + m.load.name; % user-specified name
	else
		n = ""; % default
	end
	load(n, 'd'); % load
	
function [d, m] = doPlot(d, m) % plot data in a single axes

	% Set plotting function
	if isfield(m.plot, 'funFun') % set plotting function
		hFun = m.plot.funFun; % user-provided function handle
	else % default
		if isfield(m.plot, 'arg') % plot arguments
			arg = m.plot.arg; % user-supplied arguments
		else
			arg = {}; % default
		end
		hFun =  @(d, m)plot(d.(m.x), d.(m.y), arg{:}); % function handle
	end
	fun = func2str(hFun); % convert handle to char
	i = strfind(fun, ')'); i = i(1) + 1; % starting index of function name
	fun = fun(i: end); % trim leading arguments
	i = strfind(fun, '('); i = i(1) - 1; % end index of function name
	fun = fun(1: i); % function name

	%	Set axes grouping
	if ~ isfield(m.plot, 'group') % m.plot.group is undefined
		m.plot.group = {}; % set it to empty *** better: don't call label functions
	end
	m.plot.group = string(m.plot.group); % convert char to string

	% Set line grouping
	if isfield(m.plot, 'line') % m.plot.line is defined
		line = string(m.plot.line); % internal reference is to line
	else
		line = missing; % default
	end
	
	% Initialise axes
	if strcmp(fun, 'polarplot') % check axes type
		[ax, i] = iniAxes(m, 'polar'); % polar axes
	else
		[ax, i] = iniAxes(m); % everything else
	end
	labelFigure(i, m.plot.group, line); % label figure with variable names
	labelAxes(d, m.plot.group); % label axes with variable values
	hold on; % add plots, don't replace
	
	% Loop over lines
	if ismissing(line) % no line grouping
		iLine = ones(height(d), 1); % all rows
		gs = 1; % one group
	else % line grouping is required
		[dLine, ~, iLine] = unique(d(:, line), 'rows'); % group into lines
		gs = height(dLine); % number of groups
	end
	for i = 1: gs % loop over lines
		pC = d(iLine == i, :); % data for this line
		hFun(pC, m); % plot
	end
	
	% Add legend
	if ~ ismissing(line) % there is line grouping
		t = makeLegend(dLine); % make the legend text
		hL = legend(t, 'location', 'best'); % show legend
		legend('boxOff'); % hide box and background
		if m.group == 1, m.handle.legend = []; end % initialise list of handles
		m.handle.legend = [m.handle.legend, hL]; % add to list
	end
	
	% Store handles
	if m.group == 1, m.handle.axes = []; end
	m.handle.axes = [m.handle.axes, ax]; % store axes handle(s) in m
	
	% Finalise
	if ~ strcmp(fun, 'polarplot') % can't label polar axes
		labelAxis(d, m); % label axes
	end
	hold off; % revert to original state
	prettyPlot(ax); % improve plot readability

function [d, m] = doPred(d, m)

% Evaluate model

	if isfield(m.pred, 'x'), xName = m.pred.x; else, xName = m.x; end
	if isfield(m.pred, 'y'), yName = m.pred.y; else, yName = m.y; end
	mC = m.model{m.group}; % current model
	if isfield(m.pred, 'row') % data are contained in a single row
		x = d.(xName); % predictor values: 1 x os x vs
		x = shiftdim(x, 1); % predict requires observations x variables: os x vs
		y = predict(mC, x); % predict
		d.(yName) = shiftdim(y, -1); % 1 x os
	else % data are contained in columns
		if all(size(mC.Variables) == size(d)) % fit is based on table d
			d.(yName) = predict(mC, d); % predict
		else % fit is based on x, y
			d.(yName) = predict(mC, d.(xName)); % predict
		end
	end

function [d, m] = pub(d, m) % format graphic for publication, and export

	%	Initialise
	drawnow; % finish figure updates to prevent confusing exportgraphics
	f = gcf; % handle of current figure
	a = m.handle.axes; % handles of all the axes in the current figure
	p = fieldnames(m.pub); % names of graphics properties to be set
	if isfield(m.pub, 'dest') % export destination
		dest = m.pub.dest; % user setting
	else
		dest = 'pdf'; % default
	end

	%	Format for publication	
	for pC = string(p)' % loop over properties
		switch pC % choose property
			case 'axisLength' % length of axes
				for aC = a % loop over axes
					pos = aC.Position; % axes position: [left, bottom, width, height]
					aC.Position = [pos(1: 2), m.pub.axisLength]; % set axes lengths
				end
			case 'lineWidth' % width of axes, lines
				set(a, 'LineWidth', m.pub.lineWidth); % set axes width
				l = findobj(f, 'type', 'line'); % handles of all lines in figure
				set(l, 'LineWidth', m.pub.lineWidth); % set line widths
			case 'units' % unit of distance
				set(a, 'Units', m.pub.units); % set units for lengths and widths
		end
	end
	
	%	Export graphic
	switch dest % choose export destination
		case 'copy' % copy figure to clipboard
			copygraphics(a(end), contentType = 'vector', preserveAspectRatio = 'on');
		case 'eps' % write to EPS file
			exportgraphics(f, 'Export.eps', contentType = 'vector', ...
				preserveAspectRatio = 'on', colorspace = 'cmyk');
		case 'pdf' % write to PDF file
			exportgraphics(f, 'Export.pdf', contentType = 'vector', ...
				preserveAspectRatio = 'on');
	end

function [d, m] = reduce(d, m) % reduce the d.resp dimensions and/or elements

% Inputs:
%		m.reduce: structure containing names and values of variables to be retained,
%			for example m.reduce.loc = [0, 0; 1, 1; 2, 2]
%		d.var: variables to be reduced must have size 1 x rows x cols, where
%			cols can be 1
%	Optional inputs:
%		m.reduce.squeeze = 1: remove response dimensions with a single elemebt
%	Outputs:
%		d.resp: response dimensions named in m.reduce are reduced to the required
%			elements
%		d.var: only required values of named variables are retained

	% Find variables to reduce
	dim = d.Properties.CustomProperties.RespDim; % response dimensions: 1 x ds
	i = ismember(dim, fieldnames(m.reduce)); % truth of resp. dim. to reduce
	v = dim(i); % names of variables to reduce
	vs = length(v); % number of variables to reduce
	i = find(i); % indices of response dimensions to reduce
	
	% Reduce
	r = d.(m.z); % response: ds dimensions
	s = size(r); % response array size: 1 x ds
	for j = 1: vs % loop over variables
		
		% Prepare to remove unwanted variable values
		vC = v{j}; % name of current variable
		val = d.(vC); % values of current variable
		val = shiftdim(val, 1); % remove leading single-element dimension
		valC = m.reduce.(vC); % values to keep
		
		%	Remove unwanted variable values
		if isnumeric(val) % values are numerical
			if size(valC, 2) ~= size(val, 2) % knnsearch requires same number of col.
				valC = valC'; % assume that this is a match
			end
			k = knnsearch(val, valC); % indices of values to keep
			ks = length(k); % number of values to keep
		else % values are not numerical
			k = contains(val, valC); % indices of values to keep
			ks = length(find(k)); % number of values to keep
		end
		val = val(k, :); % reduce variable values
		d.(vC) = shiftdim(val, -1); % store
		
		% Remove unwanted elements in d.resp
		l = i(j); % dimension to be reduced
		sPrec = s(1: (l - 1)); % size of preceding dimensions
		sFoll = s((l + 1): end); % size of following dimensions
		r = reshape(r, [prod(sPrec), s(l), prod(sFoll)]); % isolate reduction dim.
		r = r(:, k, :); % reduce
		s = [sPrec, ks, sFoll]; % reduced size
		
	end
	
	% Finish up
	if isfield(m.reduce, 'squeeze') % remove single-element dim. if required
		i = ~ (s == 1); i(1) = true; % first and multi-element dimensions
		s = s(i); dim = dim(i); % new size and response dimensions
	end
	if numel(s) > 1 % reshape response
		r = reshape(r, s); % match response r with reduced size, s
	end
	d.(m.z) = r; % store
	d.Properties.CustomProperties.RespDim = dim; % update dimensions

function [d, m] = doSave(d, m) % save data table to file

	if isfield(m.save, 'folder') % m.save.folder exists
		f = [m.save.folder, filesep]; % folder in which to save
	else
		f = []; % default
	end
	if isfield(m.save, 'name') % m.save.name exists
		n = [f, m.save.name]; % user-specified name
	else
		n = [f, 'd.mat']; % default name
	end
	save(n, 'd'); % save data table

function [d, m] = select(d, m) % select specified values for specified variables
		
% Input:
%		Values are specified by fields in m.select. Fields and
%		variables are matched by having the same name. Allowable values are:
%			char-based: categorical, cell or char
%			double: single values or a closed range

	% Initialise
	var = d.Properties.VariableNames; % names of variables in d
	field = fieldnames(m.select); % names of fields in m.select
	name = intersect(var, field); % names in both var and field
	names = length(name); % number of names
	rows = size(d, 1); % number of rows in d
	row = zeros(rows, names); % rows to select

	% For each name, find rows to select
	for i = 1: names % loop over names
		nameC = name{i}; % current name
		valV = d.(nameC); % values of variable
		valF = m.select.(nameC); % value of field
		switch class(valV) % class of variable
			case {'categorical', 'cell', 'char', 'string'}
				rowC = ismember(valV, valF);
			case 'double'
				if isvector(valV) && isvector(valF) % variable and field both vectors
					if size(valV, 2) == 1 % variable is a column vector
						valF = valF(:); % make field a column vector
					end
				end
				rowC = ismember(valV, valF, 'rows');
			otherwise
				warning('Row selection: class %s not supported', class(valV));
		end
		row(:, i) = rowC; % store result for this name
	end

	% Select rows
	row = all(row, 2); % find rows that satisfy all tests
	d = d(row, :); % select rows

function [d, m] = doSet(d, m) % set the properties of graphics objects

	drawnow; % update graphics, ready for setting them
	type = fieldnames(m.set); % object types
	a = m.handle.axes; % axes handles
	as = length(a); % number of axes
	for j = type' % loop over types
		switch char(j) % current type
			case 'axes' % set axes properties
				set(m.handle.axes, m.set.axes{:});
			case 'axesLength' % set length of axes
				for i = 1: as
					hC = m.handle.axes(i); % current handle
					set(hC, 'units', 'centimeters');
					pos = get(hC, 'position'); % current position and size
					pos = [pos(1: 2), m.set.axesLength]; % current position, new size
					set(hC, 'position', pos); % location and size
				end
			case 'caxis' % set colormap limits
				for i = 1: as
					clim(m.handle.axes(i), m.set.caxis);
				end
			case 'colorbar' % display a colorbar for each axes
				for i = 1: as
					c = colorbar(m.handle.axes(i)); % create colorbar
					arg = m.set.colorbar; % colorbar arguments
					if ~ isempty(arg) % there are arguments
						p = reshape(arg, 2, [])'; % name-value pairs: ps x 2, cell
						ps = size(p, 1); % number of pairs
						for k = 1: ps % loop over pairs
							n = p{k, 1}; v = p{k, 2}; % current name and value
							switch n % which type of text?
								case 'label'
									c.Label.String = v; % label string
								case 'title'
									c.Title.String = v; % title string
								otherwise
									set(c, n, v); % everything else
							end
						end
					end
				end
			case 'colormap' % set colormap identity
				for i = 1: as
					colormap(m.handle.axes(i), m.set.colormap);
				end
			case 'legend' % set legend properties
				for i = 1: as
					legend(m.handle.axes(i), m.set.legend{:});
				end
			case 'line' % set line properties
				l = findobj(a, 'type', 'line'); % line handles in all axes
				switch m.set.line % check for special cases
					case 'fill' % fill markers with same colour as marker edges
						set(l, {'markerFaceColor'}, {l.Color}');
					otherwise
						set(l, m.set.line{:}); % standard case
				end
			case 'lineWidth' % set line width
				h = m.handle.axes; % axes handles
				wid = m.set.lineWidth; % line width (pt)
				set(h, 'lineWidth', wid); % line width (pt)
				hLine = findobj(h, 'type', 'line'); % handles of all lines
				set(hLine, 'lineWidth', wid); % set line width
			case 'title' % set axes title properties
				for i = 1: as % loop over axes
					aC = a(i); % current axes
					set(aC.Title, m.set.title{:}); % set properties
				end
		end
	end
	
function [d, m] = show(d, m) % display models

	for i = 1: length(m.model)
		disp(m.model{i});
	end

function [d, m] = unpack(d, m) % unpack specified variables

% Inputs:
%		m.unpack.var: names of variables to be unpacked, for example
%			m.unpack.var = {'time', 'loc'}
%		d.var: variables to be unpacked, for example d.time, d.loc:
%			must have size 1 x rows x cols, where cols can be 1
%		d.(m.z): response variable
% Optional inputs:
%		m.unpack.varName: new variables to split variable varName, for example
%			m.unpack.loc = {'x', 'y'}
%	Outputs:
%		d.var: changed from 1 x rows x cols to rows x cols
%		d.x, d.y: these are added if m.unpack.var includes 'loc'
%		d.(m.z): unpacked dimensions are removed
%	Method:
%		Multi-column variables are replaced with single-column indices and
%		unpacked. Then the multi-column values are replaced.

	% Find variables to unpack
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	if ~ isfield(m.unpack, 'var'), return, end % there is no variable to unpack
	i = ismember(dim, m.unpack.var); % true for response dimensions to unpack
	v = dim(i); % names of variables to unpack
	vs = length(v); % number of variables to unpack
	if vs == 0, return, end % nothing to do
	iUn = find(i); % indices of response dimensions to unpack
	iRet = find(~ i); % indices of response dimensions to retain
	iRet(1) = []; % leading dimension, which is single-element, is not retained
	val = cell(1, vs); % allocate storage for values of variables: 1 x vs
	valN = val; % allocate storage for values of multi-column variables: 1 x vs
	
	% Unpack variables
	for i = 1: vs % loop over variables to unpack
		vC = v{i}; % name of current variable
		valC = d.(vC); % values of current variable: 1 x rs x cs
		valC = shiftdim(valC, 1); % numerical values: rs x cs
		s = size(valC); % size of values
		if s(2) == 1 % variable is single column
			val{i} = valC; % store: rs x 1
		else % variable is multi-column
			val{i} = (1: s(1))'; % replace values with linear index: rs x 1
			valN{i} = valC; % store values: rs x cs
		end
	end
	[val{:}] = ndgrid(val{:}); % expand values
	rows = numel(val{1}); % number of rows in unpacked table
	
	% Unpack response
	r = d.(m.z); % response: first dimension is assumed to have a single element
	s = ones(1, length(dim)); % initialise size for trailing single elements
	sC = size(r); ss = length(sC); % response size
	s(1: ss) = sC; % response size, including trailing single elements
	s = s(iRet); % size of retained dimensions
	r = permute(r, [iUn, iRet, 1]); % shift dimensions to left for unpacking
	r = reshape(r, [rows, s]); % unpack response
	dim(iUn) = []; d.Properties.CustomProperties.RespDim = dim; % update
	
	% Store
	d = repmat(d, [rows, 1]); % new number of rows equals number of values
	for i = 1: vs % loop over variables
		vC = v{i}; % name of current variable
		valC = val{i}; % values or indices of current variable
		valC = valC(:); % expand values or indices: rows x 1
		if isempty(valN{i}) % variable is single column
			d.(vC) = valC; % store as column vector: rows x 1
		else % variable is multi-column
			valC = valN{i}(valC, :); % restore the values
			d.(vC) = valC; % store as multi-column variable: rows x cols
			if isfield(m.unpack, vC) % split variable into single-column var.
				d = splitvars(d, vC, 'newVariableNames', m.unpack.(vC)); % split
			end
		end
	end
	d.(m.z) = r; % store response

function g = group(d, col) % group data table into rows with unique parameters
%
% Input:
%   d = data table
%   col = cell array of names of parameters that are constant within a group
% Output:
%   g = d split into a cell array, with one member for each group

	if isempty(col) % no grouping
		g{1} = d;
	else
		u = unique(d(:, col)); % unique rows
		gs = size(u, 1); % number of groups
		g = cell(1, gs); % allocate storage
		for i = 1: gs % loop over groups
			r = ismember(d(:, col), u(i, :)); % find all rows in this group
			g{i} = d(r, :); % assign them to this group
		end
	end
	
function [h, i] = iniAxes(m, p) % initialise new axes

% Input:
%   m = metadata structure
%		p = 'polar' for polar axes, rather than Cartesian: optional
% Output:
%   h = axes handle
%   i = axes index within figure
%	Notes:
%		plot formats are 'basic', 'paper', 'poster', 'talk', or 'thesis'

	% Default input arguments
	if isfield(m.plot, 'format') % plot format
		form = m.plot.format;
	else % default
		form = 'basic';
	end
	if isfield(m.plot, 'subplot') % subplot layout
		layout = m.plot.subplot; % user-defined
	else
		layout = [2, 2]; % default
	end
	if exist('p', 'var') % axes type is specified
		a = p; % polar axes
	else
		a = 'cart'; % Cartesian axes
	end
	axesPerFig = prod(layout); % number of axes per figure
	
	% Create axes
	i = m.group; % current group and, therefore, axes number
	i = mod(i - 1, axesPerFig) + 1; % axes index on this figure
	if i == 1
		figure('windowStyle', 'docked'); % new figure
	end
	switch a % check axes type
		case 'cart' % Cartesian axes
			h = subplot(layout(1), layout(2), i); % create Cartesian axes
		case 'polar' % polar axes
			h = polaraxes; % create polar axes
			subplot(layout(1), layout(2), i, h); % add axes to figure
	end

	% Set formatting parameters
	switch form
		case 'paper'
			set(h, 'fontSize', 16); % axis label size (pt)
		case 'poster'
			set(h, 'fontSize', 28);
		case 'talk'
			set(h, 'fontSize', 16);
		case 'thesis'
			set(h, 'fontSize', 12);
	end

function labelAxes(d, gAxes) % label axes with values of axes variables

	vs = length(gAxes); % number of variables
	if vs == 0 % no label required
		return
	end
	s = cell(1, vs); % one string for each variable
	for i = 1: vs
		val = d.(gAxes{i})(1, :); % value of current axes variable
		s{i} = makeString(val); % convert to char
	end
	s = sprintf('%s, ', s{:}); % concatenate
	s = s(1: end - 2); % remove trailing punctuation
	title(s);

function labelAxis(d, m) % label axes

	n = d.Properties.VariableNames; % variable names
	des = d.Properties.VariableDescriptions; % variable descriptions
	for a = ["x", "y", "z"] % loop over possible axes
		if isfield(m, a) % axis exists
			l = m.(a); % default label
			i = strcmp(l, n); % index of axis in variable names
			if any(i) % axis has a variable
				l = n{i}; % label axis with variable name
				if ~ isempty(des) % description property exists
					dC = des{i}; % description of this variable
					if ~ isempty(dC) % description exists
						l = des{i}; % use description for axis label
					end
				end
			end
			f = a + "label"; % name of labelling function
			feval(f, l); % label axis
		end
	end

function labelFigure(i, gAxes, gLine) % axes variables followed by line var.

	if i > 1 % figure is already labelled
		return
	end
	sAxes = []; sLine = []; % initialise
	if ~ isempty(gAxes) % there is at least one axes variable
		sAxes = sprintf('%s, ', gAxes{:}); % list variables
		sAxes = sAxes(1: end - 2); % remove trailing punctuation
	end
	if ~ ismissing(gLine) % there is at least one line variable
		sLine = sprintf('%s, ', gLine{:}); % list variables
		sLine = sLine(1: end - 2); % remove trailing punctuation
	end
	if ~ isempty(sAxes) && ~ isempty(sLine) % there are both axes and line var.
		s = [sAxes '; ' sLine]; % concatenate
	else
		s = [sAxes, sLine]; % no need for separator
	end
	if ~ isempty(s) % label figure
		sgtitle(s); % label
	end

function labelLine(d, h, gLine) % add line's display name

	vs = length(gLine); % number of variables
	if vs == 0 % no label required
		return
	end
	s = cell(1, vs); % one string for each variable
	for i = 1: vs
		val = d.(gLine{i})(1, :); % value of current line variable
		s{i} = makeString(val); % convert to char
	end
	s = sprintf('%s, ', s{:}); % concatenate
	s = s(1: end - 2); % remove trailing punctuation
	set(h, 'displayName', s); % display name

function l = makeLegend(d) % make text for legend

% Input: d = table containing values for each row of the legend

	[rs, vs] = size(d); % number of rows and variables
	v = d.Properties.VariableNames; % variable names
	l = cell(rs, 1); % one member for each row
	for i = 1: rs % loop over rows
		s = cell(1, vs); % one member for each variable
		for j = 1: vs % loop over variables
			vC = v{j}; % name of current variable
			val = d.(vC)(i, :); % current value
			s{j} = makeString(val); % current text
		end
		s = sprintf('%s, ', s{:}); % concatenate
		s = s(1: end - 2); % remove trailing punctuation
		l{i} = s; % store
	end
	
function s = makeString(val) % convert value, numeric or not, to char

	% *** Change to string()? ***
	if isnumeric(val) || islogical(val) % numeric value, including row vector
		s = num2str(val, '%.3g '); % convert to char
	else % char, string, cell or categorical
		s = char(val); % convert to char if necessary
	end

function prettyPlot(h) % improve the readability of an existing plot

% Input:
%   h = axes handle(s)

	if isprop(h, 'plotBoxAspectRatio') % set aspect ratio if possible
		h.PlotBoxAspectRatio = [1, 1, 1]; % square aspect ratio
		h.Box = 'off'; % no enclosing box around axes
	end
	h.TickLength = [.02, .02]; % make ticks more visible by enlarging them
	h.FontName = 'Helvetica'; % easily read font