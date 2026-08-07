function m = setHandle(m) % set handles for task functions

	%	List task functions
	taskC = string(m.tasks); % sequence of current task functions: 1 x 1, string
	taskC = split(taskC)'; % individual task names: 1 x ts, string
	hand = runPrimTask'; % handles of all task functions: 1 x hs, cell
	taskA = strings; % allocate storage: string
	for i = 1: length(hand) % loop over handles
		hC = hand{i}; % current handle
		tC = func2str(hC); % name of task
		taskA(i) = tC; % store task name
	end

	%	Assign handles to current tasks
	for tC = taskC % loop over current tasks
		i = ismember(taskA, tC); % index of current task in list of all tasks
		if any(i) % task found
			m.(tC).fun = hand{i}; % store task handle in metadata
		elseif contains(tC, ["max", "mod", "pca"]) % can't overlay Matlab functions
			s = char(tC); s1 = upper(s(1)); s2 = s(2: end); % capitalise first char.
			tDoName = "do" + s1 + s2; % e.g., replace max with doMax
			m.(tC).fun = hand{taskA == tDoName}; % set handle
			if isfield(m.(tC), "group") % fix grouping
				m.(tDoName).group = m.(tC).group; % grouping
			end
		end
	end
