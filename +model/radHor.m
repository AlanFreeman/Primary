function [d, m] = radHor(~, m) % fit horizontal cell dendritic radius

	% Obtain and condition the horizontal cell dendritic data
	folder = m.folder + filesep + "Wässle (89) Horizontal"; % data folder
	d = prim.readFolder(folder); % create table from folder
	d = renamevars(d, "label", "type"); % rename label
	d.type = categorical(d.type); % make it categorical
	ecAna = m.p.magRet * d.ecc; % anatomical eccentricity (deg)
	ecFun = model.readTab("ecFunc", ecAna); % convert ana. to func. ec. (deg)
	d.ec = ecFun; % store functional eccentricity (deg): es x 1
	d = sortrows(d, "ec"); % sort on eccentricity

	%	Calculate radius and store
	r = 10 .^ d.areaLog; % dendrite area (mm^2)
	r = sqrt(r / pi); % radius (mm): es x 1
	d.radius = m.p.magRet * r; % radius (deg): es x 1
	d.Properties.VariableDescriptions("type") = "Cell type"; % describe variable
	d.Properties.VariableDescriptions("ec") = "Functional eccentricity (deg)";
	d.Properties.VariableDescriptions("radius") = "Radius (deg)";