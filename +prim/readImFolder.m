function [folder, name] = readImFolder(file) % list names of image files

	%	Find name of folder containing image files
	if computer('arch') == "maca64" % local computer
		folder = [userpath, '/Data/Prim/Macaque n02487347']; % folder
	else % remote
		folder = './Macaque n02487347'; % folder
	end
	folder = string(folder); % name of folder containing image files: string

	% List files in folder, remove those with names starting with dot
	f = dir(folder); % files in folder: structure
	n = {f.name}; % file names: 1 x fs, cell
	i = startsWith(n, '.'); % file names starting with .
	f(i) = []; n(i) = []; % remove them

	%	List image files
	if isempty(file) % set number of files
		fs = length(f); % all files
		file = 1: fs;
	end
	n = n(file); % file name: 1 x fs, cell
	n = string(n); % file name: 1 x fs, string

	%	Sort the file names
	switch 'sort' % sort method
		case 'num' % numerical, as in Mac OS directory listings
			num = extractBetween(n, "_", "."); % extract file number
			num = double(num); % convert from string to numerical
  		[~, i] = sort(num); % indices resulting from numerical sort
  		name = n(i); % return names in numerical order
		case 'sort' % alphanumeric, as in Matlab's sort function
			name = sort(n); % ensure consistency across operating systems
	end
