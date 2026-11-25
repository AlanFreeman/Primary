function [d, m] = bpi(~, m)

%	Plot "bandpass index": Lee (12)

	%	Obtain data
	f = m.folder + filesep + "Lee (12)" + filesep + "Lee"; % data file
	load(f, "d"); % load data file
	x = d.bpi; % bandpass index: bs x 1
	%	y = d.gen; % geniculate data

	%	Adjust data
	x = x - 0.05; % bin starting points: bs x 1
	d.bpi = x; % store