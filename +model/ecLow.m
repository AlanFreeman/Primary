function [d, m] = ecLow(~, m) % func. v. anat. ec., low range: McGregor (18)

	%	Read empirical data
	folder = m.folder + filesep + "McGregor (18)"; % data folder
	d = prim.readFolder(folder); % create table from folder
	d = renamevars(d, "label", "source"); % rename label
	d.source = categorical(d.source); % make it categorical
	d.ecAnaMm = d.eccUm / 1000; % anatomical eccentricity (mm)
	d.ecFuncMm = d.eccFunUm / 1000; % functional eccentricity (mm)
	d = sortrows(d, ["source", "ecAnaMm"]); % sort by source, anatomical ec.

	%	Store
	d.ecAna = m.p.magRet * d.ecAnaMm; % anatomical eccentricity (deg)
	d.ecFunc = m.p.magRet * d.ecFuncMm; % functional eccentricity (deg)
	d.Properties.VariableDescriptions("source") = "Source";
	d.Properties.VariableDescriptions("ecAna") = "Anatomical eccentricity (deg)";
	d.Properties.VariableDescriptions("ecFunc") = "Functional eccentricity (deg)";