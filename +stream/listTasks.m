function list = listTasks(namespace, taskFile, m)

% listTasks: list handles of functions in task sequence
%
%	Method: task functions can appear in three libraries:
%	*	namespace = folder containing functions specific to a run file. For example,
%		the namespace for runPrim is +prim. Specify namespace without the +.
%	*	taskFile = file containing task functions specfic to a run file.
%		For example, the file for runPrim is runPrimTask.
%	*	streamFile = file containing task functions run by stream, e.g., select
%	Input:
%		namespace = namespace in which to find task functions
%		taskFile = file in which to find task functions
%		m = metadata, including the task sequence m.tasks
%	Output:
%		list = table with variables:
%			task = task name, e.g. "max": string
%			hand = task function handle, e.g. @doMax: function handle
%			group = grouping variable: string

	%	List the handles of namespace functions
	namespace = string(namespace); % ensure that namespace is a string
	d = userpath + "/Drive/+" + namespace + "/*.m"; % namespace directory
	f = dir(d); % functions in the directory: fs x 1, structure
	fs = length(f); % number of functions
	nNamespace(fs) = ""; % allocate storage: 1 x fs, string
	hNamespace = cell(fs, 1); % allocate storage: fs x 1, cell
	for i = 1: fs % loop over functions
    	[~, n] = fileparts(f(i).name); % function file name without extension
			nNamespace(i) = string(n); % function name
    	n = namespace + "." + n; % Matlab-accessible name: string
    	hNamespace{i} = str2func(n); % function handle
	end

	%	List the handles of task functions
	if ~ isempty(taskFile) % task file specified
		h = str2func(taskFile); % task file as a handle
		hTask = h(); % list the handles of functions in task file: hs x 1, cell
	else % no task file
		hTask = {}; % empty cell
	end
	hStream = stream(); % list handles of func. in stream file: hs x 1, cell
	handLib = [hTask; hStream]; % handle library: hs x 1, cell
	hs = length(handLib); nameLib(hs) = ""; % allocate storage
	for i = 1: hs % loop over handles
		nameLib(i) = func2str(handLib{i}); % names of library functions: 1 x hs
	end

	%	Combine name and handle lists
	nameLib = [nNamespace, nameLib]; % combine name lists: 1 x fs, string
	handLib = [hNamespace; handLib]; % combine handle lists: fs x 1, cell

	%	Initialise the task sequence
	t = m.tasks; % task sequence
	t = string(t); % ensure that the task sequence is a string
	t = split(t); % split the sequence into a list of tasks: ts x 1, string
	i = matches(t, "break"); % point at which to break the sequence
	if any(i) % there is a break task
		i = find(i); % index of break task
		t = t(1: i - 1); % discard the remaining tasks
	end
	nameTask = t; % task list: ts x 1, string

	%	Create a table of tasks
	ts = length(nameTask); % number of tasks
	hand = cell(ts, 1); % task handles
	group = cell(ts, 1); % grouping variables
	list = table(nameTask, hand, group, ...
		variableNames = ["task", "hand", "group"]); % ts x 2

	%	Fill in the handles of functions in the task sequence
	for i = 1: ts % loop over tasks

		% Manage task names that duplicate Matlab function names
		tC = nameTask(i); % name of current task: string
		tHand = tC; % name of function to be executed: string
		tDup = [ ... % tasks that duplicate Matlab functions
			"anova", "export", "fitglm", "fitlm", "fitnlm", "load", "max", "mod", ...
			"pca", "plot", "pred", "print", "save", "set"];
		if ismember(tC, tDup)
			s1 = extract(tC, 1); % first character of task name
			s2 = extractAfter(tC, 1); % remaining characters
			tHand = "do" + upper(s1) + s2; % e.g. "print" becomes "doPrint"
		end

		%	Fill in handles
		j = matches(nameLib, tHand); % find matching function name in library
		if any(j) % there is a match
			h = handLib{j}; % required handle
			list.hand(i) = {h}; % fill in task handle
		else
			error("Task ""%s"" is undefined", tC); % throw error
		end

		%	Fill in grouping variables
		if isfield(m, tC) && isfield(m.(tC), "group") % this task uses grouping
			list.group{i} = string(m.(tC).group); % store in list
		end
		if tC == "add" % special case: "add" uses the same grouping as "plot"
			if ~ isfield(m, "plot")
				error("Task ""add"" requires previous task ""plot"""); % throw error
			elseif isfield(m.plot, "group")
				list.group{i} = string(m.plot.group); % add grouping for task "add"
			end
		end

	end