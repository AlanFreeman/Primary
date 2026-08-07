function index = indexStruct(d, m, name) % calculate structural index

%	Input:
%		name of index: "cdi", "sdi"
%	Output:
%		index: structure with fields L, M, disp: bs x 1, complex

	%	Initialise
	locC = m.p.cone.loc; % cone location (deg): cs x 2
	typeC = m.p.cone.type; % cone type: cs x 1
	locB = shiftdim(d.loc(1, :, :), 1); % cortical cell location (deg): bs x 2
	
	switch string(name) % select name of index
		case "cdi" % centre-of-mass displacement index

			%	Loop over cone types
			indexC = zeros(size(locB, 1), 2); % allocate storage: bs x 2
			for i = 1: 2 % L cones then M cones

				%	Calculate cone weights from cone-beta displacements
				locCC = locC(typeC == i, :); % current cone location (deg): cs x 2
				%	w = m.calConv(locCC, locB, m.p.radBetaCen); % weight matrix: bs x cs
				w = m.calConv(locCC, locB, m.p.radCen); % weight matrix: bs x cs
				%	prim.see(m, scatter = {locCC, locB, m.p.c.b}, weight = w); % debug

				%	Calculate cone-beta displacements
				locCP = permute(locCC, [3, 1, 2]); % prepare for subtraction: 1 x cs x 2
				locBP = permute(locB, [1, 3, 2]); % bs x 1 x 2
				disp = locCP - locBP; % disp. of cone from beta cell (deg): bs x cs x 2
				disp = disp(:, :, 1) + 1i * disp(:, :, 2); % complex disp.: bs x cs

				%	Calculate index
				indC = dot(w, disp, 2); % weighted displacement (deg): bs x 1
				indexC(:, i) = indC ./ sum(w, 2); % norm. by weights (deg): bs x 1

			end

			%	Store
			index.L = indexC(:, 1); % L-cone cdi (deg): bs x 1, complex
			index.M = indexC(:, 2); % M-cone cdi (deg): bs x 1, complex
			index.disp = abs(index.L - index.M); % difference cdi: bs x 1
			index.dim = 'loc'; % dimensions

		case "sdi" % S-cone distance index
			locS = locC(typeC == 3, :); % S-cone location (deg): ss x 2
			i = knnsearch(locS, locB); % index of nearest S-cone: ls x 1
			dist = locB - locS(i, :); % displacement from nearest S-cone: ls x 2
			dist = sqrt(sum(dist .^ 2, 2)); % dist. to nearest S-cone (deg): ls x 1
			index = dist / max(dist); % index: ls x 1
	end