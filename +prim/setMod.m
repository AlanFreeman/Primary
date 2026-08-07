function [d, m] = setMod(iFile, ecc, m) % set modulations or eccentricity

	import prim.readTab % find function
	if ~ isempty(iFile) % file number is defined: read it

		%	Read modulation file
		file = "Modulation " + iFile + ".mat"; % modulation file
		m.p.file = file; % store
		[d, m] = m.readFile(m); % set the data folder and read the data file
	
		%	Set constants associated with modulation file
		switch iFile
			case 7 % eccentricity = 10 deg
			case 16 % ratio of off-centre to all midget ganglion cells = .546
				m.p.c.big = [.03, -.03]; % location of cell with big response (deg): 1 x 2
				m.p.ratSign = .546; % ratio of off-centre to all midget ganglion cells
			case 17 % ratio of off-centre to all midget ganglion cells = .52
				m.p.c.big = [0, -.015]; % location of cell with big response (deg): 1 x 2
				m.p.c.dirPref = 90; m.p.c.freqSPref = 9; % preferred stimuli
				m.p.ratDensBeta = 1; % ratio of inhibitory to excitatory cortical cells
				m.p.ratSign = .52; % ratio of off-centre to all midget ganglion cells
			case 18 % ratio of inhibitory to excitatory cortical cells = .25
				m.p.c.big = [-.015, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 36; m.p.c.freqSPref = 8.3; % preferred stimuli
				m.p.ratDensBeta = .25; % ratio of inhibitory to excitatory cortical cells
			case 19 %	m.p.ratDensBeta replaced by m.p.ratBetaIn
				m.p.c.big = [-.015, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 36; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 20 %	modulation matrix replaced by synaptic structure: gen_betaEx
				m.p.c.big = [-.015, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 36; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 21 %	mod'n matrix replaced by syn. structure: gen_betaEx, gen_betaIn
				m.p.c.big = [-.015, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 36; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 22 % modulations calc. for both gen_betaEx and gen_betaIn synapses
				m.p.c.big = [0, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -144; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 23 % gain calculated for betaIn_betaEx synapse, kIncIn = 0
				m.p.c.big = [0, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -144; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 24 % gain calculated for betaIn_betaEx synapse, kIncIn > 0
				m.p.c.big = [0, -.015]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -144; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 25 %	second cone stage added, time constant recalculated
				m.p.c.big = [0, -.03]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 144; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 26 %	set inhibitory-excitatory convergence radius from Packer (11)
				m.p.c.big = [0, -.03]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = 144; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 27 %	visual field width = .5 deg
				m.p.c.big = [.015, -.06]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -90; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 28 % cleaned drive: driveF, driveS, fixed rgb2lms to include gamma
				m.p.c.big = [.015, -.06]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -90; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 29 % reduced ganglion cell resting impulse rate to match Troy (94)
				m.p.c.big = [.015, -.06]; % location of cell with big resp. (deg): 1 x 2
				m.p.c.dirPref = -90; m.p.c.freqSPref = 8.3; % preferred stimuli
			case 30 % development phase 2, drifting equiluminant gratings
				% Development phase 1: Modulation 29
				m.p.c.big = [-.15, .015]; % big response to M, equiluminant, achromatic
				m.p.c.dirPref = -36; m.p.c.freqSPref = 7.7; % L, M, achromatic
				m.p.c.dirPrefEqui = 36; m.p.c.freqSPrefEqui = 6.5; % equiluminant
			case 31 % reduced ganglion cell dendritic radius at low eccentricity
				m.p.c.big = [.015, -.075]; % big response to L, achromatic
				m.p.c.dirPref = 90; m.p.c.freqSPref = 9.0; % L, achromatic
			case 32 % at low ecc.: equated g.c. with cone dens.; zeroed g.c. std. dev.
				%	updated default spatial frequency range
				m.p.c.big = [.015, -.075]; % big response to L, achromatic			
				m.p.c.dirPref = -108; m.p.c.freqSPref = 8.3; % achromatic
			case 33 % eccentricity = 1 deg; updated default spatial frequency range:
				%	no orientation selectivity for achromatic stimulus
				m.p.c.big = [.015, -.015]; % big response to achromatic			
				%	m.p.c.dirPref = -18; m.p.c.freqSPref = 7.3; % achromatic
			case 34 % randomly varied geniculocortical gain; drifting ach. grating:
				%	orientation selectivity, maximum response amplitude about 8 Hz
				m.p.c.big = [.015, 0]; % big response to achromatic			
				m.p.c.dirPref = 180; m.p.c.freqSPref = 11; % achromatic
				m.p.c.dirPrefEqui = 36; m.p.c.freqSPrefEqui = 16.5; % equiluminant
			case 35 % pulsed natural images: negligible orientation selectivity
				m.p.c.big = [.008, .008]; % big response to achromatic			
				m.p.c.dirPref = 72; m.p.c.freqSPref = 7.3; % achromatic
				%	m.p.c.dirPrefEqui = 36; m.p.c.freqSPrefEqui = 16.5; % equiluminant
			case 36 % drifting achromatic grating:
				%	orientation selective, max. resp. amplitude about 8 Hz
				m.p.c.big = [.015, 0]; % big response to achromatic			
				m.p.c.dirPref = 180; m.p.c.freqSPref = 11; % achromatic
			case 37 % drifting achromatic grating then drifting equiluminant grating
				m.p.c.bigAch = [.015, -.023]; % high pref. spatial freq. to achromatic
				m.p.c.big = [.023, -.008]; % big oriented response to equiluminance
				%	m.p.c.big = [-.015, .015]; % big oriented response to equiluminance
				m.p.c.dirPref = 36; m.p.c.freqSPref = 16.5; % achromatic
				m.p.c.dirPrefEqui = -90; m.p.c.freqSPrefEqui = 9.2; % equiluminance
			case 38 % update of 37
				m.p.c.big = [.015, -.015]; % 8.7 Hz response to achromatic
				%	m.p.c.big = [-.008, .008]; % 11.7 Hz response to equiluminance
				m.p.c.dirPref = -18; m.p.c.freqSPref = 8.2; % 8.5 Hz, achromatic
				%	m.p.c.dirPref = 108; m.p.c.freqSPref = 4.9; % 11.7, equiluminance
			case 39 % drifting achromatic grating
				m.p.c.big = [.015, -.015]; % 8.7 Hz response to achromatic
				%	m.p.c.big = [-.008, .008]; % 11.7 Hz response to equiluminance
				m.p.c.dirPref = -18; m.p.c.freqSPref = 8.2; % 8.5 Hz, achromatic
				%	m.p.c.dirPref = 108; m.p.c.freqSPref = 4.9; % 11.7, equiluminance
			case 40 % drifting achromatic grating then pulsed natural images
				m.p.c.big = [.023, .031]; % 3.7 Hz response to equiluminant grating
				m.p.c.dirPref = -144; m.p.c.freqSPref = 8.2; % 9.2 Hz, achromatic
			case 41 % achromatic development, kAdapt = .9, kGainOff = 1.1
				switch 'equi' % achromatic grating
					case 'ach'
						m.p.c.big = [-.038, -.023]; % medium resp. amp., good orient. sel.
						m.p.c.dirPref = 72; m.p.c.freqSPref = 11.4; % 4.1 Hz, achromatic
						%	m.p.c.dirPref = -162; m.p.c.freqSPref = 9.8; % 3.4 Hz, equilum't
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .023];
						%	m.p.c.dirPref = -54; m.p.c.freqSPref = 11.4; % 4.7 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 8.2; % 3.9 Hz, equiluminant
				end
			case 42 % achromatic development, kAdapt = .9, kGainOff = 1.05
				switch 'equi' % achromatic grating
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015];
						%	m.p.c.dirPref = -54; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 6.5; % 3.9 Hz, equiluminant
				end
			case 43 % achromatic development, kAdapt = .95, kGainOff = 1.05
				switch 'equi' % achromatic grating
					case 'ach'
						%	m.p.c.big = [.008, -.015]; % medium resp. amp., good orient. sel.
						%	m.p.c.dirPref = 72; m.p.c.freqSPref = 11.4; % 4.1 Hz, achromatic
						%	m.p.c.dirPref = -162; m.p.c.freqSPref = 9.8; % 3.4 Hz, equilum't
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015];
						m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 5.2 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 5.4 Hz, equiluminant
				end
			case 44 % achromatic development, kAdapt = .95, kGainOff = 1.025
				switch 'equi' % achromatic grating
					case 'ach'
						%	m.p.c.big = [.008, -.015]; % medium resp. amp., good orient. sel.
						%	m.p.c.dirPref = 72; m.p.c.freqSPref = 11.4; % 4.1 Hz, achromatic
						%	m.p.c.dirPref = -162; m.p.c.freqSPref = 9.8; % 3.4 Hz, equilum't
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015];
						%	m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 5 Hz, achromatic
						%	m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 4.9 Hz, equiluminant
				end
			case 45 % achromatic development, kAdapt = .975, kGainOff = 1.025
				switch 'equi' % achromatic grating
					case 'ach'
						%	m.p.c.big = [.008, -.015]; % medium resp. amp., good orient. sel.
						%	m.p.c.dirPref = 72; m.p.c.freqSPref = 11.4; % 4.1 Hz, achromatic
						%	m.p.c.dirPref = -162; m.p.c.freqSPref = 9.8; % 3.4 Hz, equilum't
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015];
						%	m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 5.2 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 5.5 Hz, equiluminant
				end
			case 46 % achromatic development, m.p.adapt = 0, kGainOff = 1.025
				switch 'equi' % achromatic grating
					case 'ach'
						m.p.c.big = [0, .023]; % medium resp. amp., good orient. sel.
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						%	m.p.c.dirPref = 144; m.p.c.freqSPref = 6.5; % 4.6 Hz, equilum't
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015]; % good resp. amp., poor orient. sel.
						%	m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 5.3 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 5.5 Hz, equiluminant
				end
			case 47 % phase 2 development, m.p.adapt = 0, kGainOff = 1.025
				switch 'equi' % achromatic grating
					case 'ach'
						m.p.c.big = [0, .023]; % medium resp. amp., good orient. sel.
						m.p.c.dirPref = 162; m.p.c.freqSPref = 9.8; % 4.7 Hz, achromatic
						%	m.p.c.dirPref = 162; m.p.c.freqSPref = 9.8; % 2.3 Hz, equiluminant
					case 'equi' % equiluminant grating
						m.p.c.big = [.023, .015]; % good resp. amp., poor orient. sel.
						%	m.p.c.dirPref = 144; m.p.c.freqSPref = 9.8; % 4.8 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 4.7 Hz, equiluminant
				end
			case 48 % achromatic development, m.p.adapt = 0, m.p.varyGain = 0
				% phase 1 has poorer contrast sensitivity than does file 46
				switch 'equi' % achromatic grating
					case 'ach'
						m.p.c.big = [0, .023]; % medium resp. amp., good orient. sel.
						m.p.c.dirPref = 162; m.p.c.freqSPref = 9.8; % 4.7 Hz, achromatic
						%	m.p.c.dirPref = 162; m.p.c.freqSPref = 9.8; % 2.3 Hz, equiluminant
					case 'equi' % equiluminant grating
						%	m.p.c.big = [-.015, -.046]; % good resp. amp., poor orient. sel.
						%	m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						%	m.p.c.dirPref = 54; m.p.c.freqSPref = 4.9; % 5.9 Hz, equiluminant
						m.p.c.big = [.023, .008]; % med. resp. amp., poor orient. sel.
						%	m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 4.4 Hz, achromatic
						m.p.c.dirPref = 126; m.p.c.freqSPref = 4.9; % 4.3 Hz, equiluminant
				end
			case 49 % achromatic development, m.p.kAdapt = .975, m.p.varyGain = 0
				switch 'equi' % achromatic grating
					case 'ach' % optimal locations for achromatic stimulus
						m.p.c.big = [.023, .023]; % good OS for ach., medium equi., dO = 162
						m.p.c.dirPref = 18; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 4.9; % 2.7 Hz, equiluminant
						m.p.c.big = [-.008, -.031]; % good OS for ach., equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.9 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 9.8; % 4.1 Hz, equiluminant
						m.p.c.big = [.023, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = -144; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 9.8; % 3.8 Hz, equiluminant
						m.p.c.big = [-.008, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = -18; m.p.c.freqSPref = 13.1; % 3.4 Hz, equiluminant
					case 'equi' % optimal locations for equiluminant stimulus
						m.p.c.big = [-.023, -.046]; % good OS for ach., equi., dO = 54
						m.p.c.dirPref = 180; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 8.2; % 4.8 Hz, equiluminant
						%{
						m.p.c.big = [-.008, -.038]; % good OS for ach. and equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.8 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 8.2; % 4.8 Hz, equiluminant
						m.p.c.big = [.023, .008]; % good OS for ach., poor equi., dO = 162
						m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 4.3 Hz, achromatic
						m.p.c.dirPref = 126; m.p.c.freqSPref = 4.9; % 4.1 Hz, equiluminant
						m.p.c.big = [0, .023]; % good OS for ach., medium for equi., dO = 18
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.3 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 6.5; % 4.3 Hz, equiluminant
						%}
				end
			case 50 % development with images, m.p.kAdapt = .975, m.p.varyGain = 0
				% basic properties similar to Mod 49, but resp. amp. are much reduced
				switch 'equi' % achromatic grating
					case 'ach' % optimal locations for achromatic stimulus
						m.p.c.big = [.023, .023]; % good OS for ach., medium equi., dO = 162
						m.p.c.dirPref = 18; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 4.9; % 2.7 Hz, equiluminant
						m.p.c.big = [-.008, -.031]; % good OS for ach., equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.9 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 9.8; % 4.1 Hz, equiluminant
						m.p.c.big = [.023, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = -144; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 9.8; % 3.8 Hz, equiluminant
						m.p.c.big = [-.008, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = -18; m.p.c.freqSPref = 13.1; % 3.4 Hz, equiluminant
					case 'equi' % optimal locations for equiluminant stimulus
						m.p.c.big = [-.015, -.046]; % good OS for ach., poor equi., dO = 108
						m.p.c.dirPref = -180; m.p.c.freqSPref = 9.8; % 4.8 Hz, achromatic
						m.p.c.dirPref = 72; m.p.c.freqSPref = 4.9; % 4.2 Hz, equiluminant
						m.p.c.big = [.008, .053]; % poor OS for both ach. and equi., dO = 18
						m.p.c.dirPref = -90; m.p.c.freqSPref = 6.5; % 3.8 Hz, achromatic
						m.p.c.dirPref = 72; m.p.c.freqSPref = 4.9; % 1.44 Hz, equiluminant
						m.p.c.big = [.023, .008]; % good OS for ach., poor equi., dO = 36
						m.p.c.dirPref = -36; m.p.c.freqSPref = 9.8; % .77 Hz, achromatic
						m.p.c.dirPref = -72; m.p.c.freqSPref = 4.9; % .17 Hz, equiluminant
						m.p.c.big = [0, .023]; % good OS for ach., medium for equi., dO = 18
						m.p.c.dirPref = -36; m.p.c.freqSPref = 9.8; % 4.4 Hz, achromatic
						m.p.c.dirPref = 162; m.p.c.freqSPref = 9.8; % 2.3 Hz, equiluminant
				end
			case 51 % achromatic development, m.p.kAdapt = .975, m.p.GainOff = 1.025
				switch 'ach' % achromatic grating
					case 'ach' % optimal locations for achromatic stimulus
						m.p.c.big = [-.008, -.023]; % good OS for ach., med. equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.7 Hz, achromatic
						%	m.p.c.dirPref = 54; m.p.c.freqSPref = 11.4; % 3.7 Hz, equiluminant
						%{
						m.p.c.big = [.023, .023]; % good OS for ach., medium equi., dO = 162
						m.p.c.dirPref = 18; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 4.9; % 2.7 Hz, equiluminant
						m.p.c.big = [-.008, -.031]; % good OS for ach., equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.9 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 9.8; % 4.1 Hz, equiluminant
						m.p.c.big = [.023, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = -144; m.p.c.freqSPref = 11.4; % 4.8 Hz, achromatic
						m.p.c.dirPref = -144; m.p.c.freqSPref = 9.8; % 3.8 Hz, equiluminant
						m.p.c.big = [-.008, .061]; % good OS for ach., equi., dO = 0
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = -18; m.p.c.freqSPref = 13.1; % 3.4 Hz, equiluminant
						%}
					case 'equi' % optimal locations for equiluminant stimulus
						%{
						m.p.c.big = [-.023, -.046]; % good OS for ach., equi., dO = 54
						m.p.c.dirPref = 180; m.p.c.freqSPref = 11.4; % 5.5 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 6.5; % 4.7 Hz, equiluminant
						m.p.c.big = [-.008, -.038]; % good OS for ach. and equi., dO = 54
						m.p.c.dirPref = 180; m.p.c.freqSPref = 11.4; % 5.2 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 6.5; % 4.4 Hz, equiluminant
						m.p.c.big = [.023, .008]; % good OS for ach., poor equi., dO = 180
						m.p.c.dirPref = -54; m.p.c.freqSPref = 9.8; % 5.1 Hz, achromatic
						m.p.c.dirPref = 126; m.p.c.freqSPref = 4.9; % 3.8 Hz, equiluminant
						%}
						m.p.c.big = [0, .023]; % good OS for ach., medium for equi., dO = 18
						m.p.c.dirPref = 162; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 6.5; % 4.5 Hz, equiluminant
						m.p.c.big = [-.069, -.008]; % poor OS for ach., good equi., dO = 36
						m.p.c.dirPref = 54; m.p.c.freqSPref = 6.5; % 12.9 Hz, achromatic
						m.p.c.dirPref = -90; m.p.c.freqSPref = 18; % 4.7 Hz, equiluminant
						%{
						m.p.c.big = [-.015, -.046]; % good OS for ach., poor equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.4 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 4.9; % 5.8 Hz, equiluminant
						m.p.c.big = [.023, .015]; % good OS for ach., poor equi., dO = 0
						m.p.c.dirPref = -36; m.p.c.freqSPref = 11.4; % 5.2 Hz, achromatic
						m.p.c.dirPref = 144; m.p.c.freqSPref = 4.9; % 5.5 Hz, equiluminant
						%}
				end
			case 52 % achromatic then equiluminant development; .975, 1.025
				switch 'equi' % achromatic grating
					case 'ach' % optimal locations for achromatic stimulus
						m.p.c.big = [-.008, -.023]; % good OS for ach., med. equi., dO = 54
						m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 5.7 Hz, achromatic
						m.p.c.dirPref = 54; m.p.c.freqSPref = 11.4; % 3.7 Hz, equiluminant
					case 'equi' % optimal locations for equiluminant stimulus
						m.p.c.big = [0, -.015]; % poor OS for ach., good equi., dO = 36
						m.p.c.dirPref = 144; m.p.c.freqSPref = 6.5; % 7.8 Hz, achromatic
						m.p.c.dirPref = -72; m.p.c.freqSPref = 9.8; % 8.6 Hz, equiluminant
						%{
						m.p.c.big = [.015, -.015]; % poor OS for ach., good equi., dO = 90
						m.p.c.dirPref = 180; m.p.c.freqSPref = 6.5; % 8.5 Hz, achromatic
						m.p.c.dirPref = 90; m.p.c.freqSPref = 9.8; % 7.8 Hz, equiluminant
						m.p.c.big = [.046, .031]; % poor OS for ach., good equi., dO = 72
						m.p.c.dirPref = 108; m.p.c.freqSPref = 6.5; % 8.4 Hz, achromatic
						m.p.c.dirPref = 180; m.p.c.freqSPref = 9.8; % 6.7 Hz, equiluminant
						%}
				end
			case 53 % ach. development, no off- or on-S inputs to 4CBeta
				m.p.c.big = [-.031, -.031]; % good OS
				m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 6.2 Hz
				m.p.c.a = [-.0305623, -.0305623]; % cell A for SfN 25 poster
			case 54 % achromatic then equiluminant development
				m.p.c.big = [.038, .046]; % good OS
				m.p.c.dirPref = 180; m.p.c.freqSPref = 8.2; % 6.5 Hz, achromatic
				m.p.c.b = m.p.c.big; % cell B for SfN 25 poster
				m.p.c.c = [-.0153, .0229]; % cell C for SfN 25 poster
			case 55 % ach. development, fixed rectification error in doMod
				m.p.c.big = [-.031, -.031]; % good OS
				m.p.c.dirPref = 0; m.p.c.freqSPref = 11.4; % 6.3 Hz
				m.p.c.a = m.p.c.big; % cell A in SfN 2025 poster
			case 56 % equiluminant development, fixed rectification error in doMod
				m.p.c.big = [.038, .046]; % good OS
				m.p.c.b = m.p.c.big; % cell B for SfN 25 poster
				m.p.c.c = [-.0153, .0229]; % cell C for SfN 25 poster
				m.p.c.d = [-.053, .038]; % PCA extreme
				m.p.c.dirPref = 162; m.p.c.freqSPref = 8.2; % 7.0 Hz
			case 57 % achromatic development, eccentricity = 10 deg
			case 58 % equiluminant development, eccentricity = 10 deg
			case 59 % achromatic development, eccentricity = 1 deg,
				% criterion-based v.f. width, ex. cell density, spatial frequencies
			case 60 % equiluminant development, eccentricity = 1 deg,
				% criterion-based v.f. width, ex. cell density, spatial frequencies
			case 61 % achromatic development; criterion-based v.f. width,
				% ex. cell density, spatial frequencies at arbitrary eccentricity
				m.p.c.a = [-.016, -.048]; % cell A, highest OSI
			case 62 % equiluminant development; criterion-based v.f. width,
				% ex. cell density, spatial frequencies at arbitrary eccentricity
				m.p.c.b = [-.048, .04]; % cell B, COI = .4, SOI = 1
				m.p.c.c = [-.016, -.064]; % cell C, COI = .85, SOI = 0
				m.p.c.d = [.048, .072]; % cell D, COI = .52, SOI = .63
			case 63 % achromatic development, eccentricity = 10 deg
			case 64 % equiluminant development
			case 65 % ach. dev., ecc. = 30 deg, no cropping in synaptic normalisation
			case 66 % equiluminant development
			case 67 % ach. dev., ecc. = 1 deg, about 1,000 excit'y cells at all ecc.
				m.p.c.a = [-.04, -.024]; % cell A, SOI = .95, OSI = .65
				m.p.c.dirPref = 18; m.p.c.freqSPref = 12; % 4.1 Hz
				m.p.c.b = [-.048, -.04]; % cell B, from Modulation 68
			case 68 % equiluminant development
				m.p.c.a = [-.04, -.024]; % cell A, from Modulation 67
				m.p.c.b = [-.048, -.04];
					% cell B, SOI = .91, COI = .54, OSI = .51, max = 5.6 Hz
				m.p.c.dirPref = -162; m.p.c.freqSPref = 10; % cell B, ach., 5.67 Hz
				m.p.c.c = [.024, .032]; % cell C, SOI = 0, COI = .84, OSI = 0
				m.p.c.d = [0, .04]; % cell D, midpoint of PCA
				%	m.p.c.e = [.047, .016]; % gang cell E, low SOI, COI = .56, SOI = .031
				m.p.c.e = [.047, .048]; % on-gang cell E, high SOI, COI = .53, SOI = .16
			case 69 % achromatic development, 10 deg
			case 70 % equiluminant development
			case 71 % achromatic development, 30 deg
			case 72 % equiluminant development
			case 73 % achromatic development, 1 deg, m.p.adapt = 0
				m.p.c.a = [-.04, -.024]; % cell A
			case 74 % equiluminant development
				m.p.c.a = [-.04, -.024]; % cell A, from Modulation 73
				m.p.c.b = [-.048, -.04];
					% cell B, SOI = .91, COI = .54, OSI = .51, max = 5.6 Hz *** check ***
				m.p.c.c = [.024, .032]; % cell C, SOI = 0, COI = .84, OSI = 0
		end

	else % set eccentricity and associated constants

		% Set eccentricity and data folder; set data file, d, to empty
		m.p.ecc = ecc; % eccentricity
		m.p.file = []; % no file
		[d, m] = m.readFile(m); % set the data folder

		% Set visual field width and stimulus preferences
		switch 'crit' % basis for settings
			case 'crit' % criterion-based setting
				% Field contains about 200 cones, for a convenient computation time.
				%	Set the width in setPrim > width.ecc.
				m.p.wid = readTab('width', ecc); % visual field width (deg)
			case 'guess' % guesswork
				switch ecc % field width increases with eccentricity
					case 0, m.p.wid = .1;
					case 1, m.p.wid = .2;
					case 3, m.p.wid = .5;
					case 10, m.p.wid = 1;
				end
		end
		m.p.c.big = [0, 0]; % default
		m.p.c.dirPref = m.p.dir; m.p.c.freqSPref = m.p.freqS; % default

	end
	
	%	Set stimulus constants
	m.p.c.cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1]; % contrast
	dir = linspace(-180, 180, 20 + 1); % stim. direction, closed interval (deg)
	m.p.c.dir = dir(1: end - 1)'; % stimulus direction, open interval (deg)

	%	Set spatial frequency constants
	switch 'range' % basis for settings
		case 'pref' % start with preferred spatial frequency as calc. from model
			switch m.p.ecc % ganglion cell preferred spatial frequency (cycles/deg)
				case 0, freqSPref = 11.0;
				case 1, freqSPref = 9.80;
				case 3, freqSPref = 9.18;
				case 10, freqSPref = 4.29;
			end
			m.p.freqSPref = freqSPref; % store
			f = [.5, 2] * freqSPref; %	preferred spatial frequency +-1 one octave
		case 'range' % range from 0 to spatial frequency cutoff
			freqMax = readTab('freqS', m.p.ecc); % spatial freq. cutoff (cycles/deg)
			f = [0, freqMax]; % spatial frequency range (cycles/deg)
	end
	m.p.c.freqSRange = f; % spatial frequency range (cycles/deg)
	m.p.c.freqS = linspace(f(1), f(2), 10)'; % spatial frequencies (cycles/deg)
