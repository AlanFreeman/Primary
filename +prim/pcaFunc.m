function out = pcaFunc(resp, m) % principal components analysis

%	Input:
%		resp = neuronal response indices: structure with one field per index
%		index = name of response indices: char or string
%	Ouput:
%		out = PCA outputs: structure with one field per output

	%	Choose method
	p = resp.pca.val; % PCA variables
	switch m.pca.method % choose method
		case 'pca' % principal components analysis

			%	Perform PCA
			[vec, comp, varC, tSq, frac, mu] = pca(p); % perform PCA
			out.vec = vec; % eigenvectors: ps x ps
			out.comp = comp; % principal components: ls x ps
			out.var = varC; % variance accounted for by principal components: ps x 1
			out.tSq = tSq; % Hotelling's T-squared statistic: ls x 1
			out.frac = frac; % fraction of variance accounted for: ps x 1
			out.mu = mu'; % estimated mean of each property: ps x 1
		
			%	Normalise principal components
			p = comp(:, 1); q = comp(:, 2); % first and second principal comp.: ls x 1
			pMin = min(p); pMax = max(p); % minimum and maximum of first component
			qMin = min(q); qMax = max(q); % minimum and maximum of second component
			sP = 1 / (pMax - pMin); % scaling factor so that p is in range [0, 1]
			sQ = sP; % same scaling for q
			oP = - sP * pMin; % offset for p
			oQ = .5 - sQ * mean([qMin, qMax]); % offset for q
			switch m.pca.reverse % switch axis directions?
				case 0 % no
				case 1 % reverse horizontal axis
					oP = sP * pMax; sP = - sP;
				case {2, 3} % reverse vertical, both axes (add code if required)
			end
			x = oP + sP * p; % normalised first component: ls x 1
			y = oQ + sQ * q; % normalised second component: ls x 1
			out.pc1 = x; % principal component 1: ls x 1
			out.pc2 = y; % principal component 2: ls x 1
			out.pc3 = comp(:, 3); % principal component 3: ls x 1

		case 'tsne' % t-distributed stochastic neighbour embedding

			comp = tsne(p); % first two components: ls x 2
			out.pc1 = comp(:, 1); % first component: ls x 1
			out.pc2 = comp(:, 2); % second component: ls x 1

	end