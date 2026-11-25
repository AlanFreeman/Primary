function showStim(stim) % show stimulus, both still and moving

%	Input:
%		stim = structure, with fields:
%			dur = frame duration (s)
%			stim = stimulus: xs x ys x cs x ts
%			x = x limits (deg): 1 x 2

	%	Prepare for image
	s = stim.stim; % stimulus: xs x ys x cs x ts
	s = permute(s, [2, 1, 3, 4]); % swap x and y: ys x xs x cs x ts
	s = flipud(s); % flip up for down: ys x xs x cs x ts

	% Show movie
	fs = size(s, 4); % number of movie frames
	x = stim.x; y = x; % x and y limits (deg): 1 x 2
	dur = stim.dur; % frame duration (s)
	for i = 1: fs % loop over frames
		sC = s(:, :, :, i); % current stimulus: ys x xs x cs
		image(x, y, sC); % draw image
		pause(dur); % pause between frames (s)
	end
