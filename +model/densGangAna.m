function [d, m] = densGangAna(~, m) % gang. cell dens. v. ana. ec.: Wässle (89)

	% Load the empirical data, and convert mm to degrees
	folder = m.folder + filesep + "Wässle (89)"; % data folder
	d = prim.readFolder(folder); % create table from folder
	d = renamevars(d, "label", "source"); % rename label
	d.source = categorical(d.source); % for better listing
	i = d.source == "func"; % functional eccentricity
	d = d(~ i, :); % keep only anatomical eccentricity
	ecAna = m.p.magRet * d.ecc; % anatomical eccentricity (deg)
	dens = 10 .^ d.densLog; % density (mm^-2)
	dens = (m.p.magRet) ^ -2 * dens; % density (deg^-2)

	% Correct for shrinkage, and store
	i = d.source == "anatNarr"; % anatWide data are already corrected
	s = .9; % linear shrinkage
	ecAna(i) = (1 / s) * ecAna(i); % shift away from fovea
	dens(i) = s ^ 2 * dens(i); % reduce density
	d.ecAna = ecAna; d.dens = dens; % store
	d = sortrows(d, "ecAna"); % sort by eccentricity
	d.Properties.VariableDescriptions("ecAna") = "Anatomical eccentricity (deg)";
	d.Properties.VariableDescriptions("dens") = "Ganglion cell density (deg^-^2)";