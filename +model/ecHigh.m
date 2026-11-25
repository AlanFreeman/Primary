function [d, m] = ecHigh(~, m) % func. v. anat. ec., high range: Schein (88)

	%	Read empirical data
	folder = m.folder + filesep + "Schein (88)"; % data folder
	d = prim.readFolder(folder); % create table from folder
	d = renamevars(d, "label", "source"); % rename label
	d.source = categorical(d.source); % make it categorical
	i = d.source == "total"; % rows containing total offset
	d = d(i, :); % use total offset (mm): es x 1
	d.source(:) = "schein"; % indicate author
	d = sortrows(d, "eccMm"); % sort by source, eccentricity

	%	Calculate anatomical eccentricity
	ecFunc = m.p.magRet * d.eccMm; % functional eccentricity (deg): es x 1
	offset = m.p.magRet * d.offsetMm; % offset (deg)
	ecAna = ecFunc + offset; % anatomical eccentricity (deg)
	
	%	Store
	d.ecAna = ecAna; % anatomical eccentricity (deg): es x 1
	d.ecFunc = ecFunc; % functional eccentricity (deg)
	d.Properties.VariableDescriptions("source") = "Source";
	d.Properties.VariableDescriptions("ecAna") = "Anatomical eccentricity (deg)";
	d.Properties.VariableDescriptions("ecFunc") = "Functional eccentricity (deg)";