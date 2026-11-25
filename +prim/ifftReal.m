function r = ifftReal(r, dim) % calculate real part of inverse transform

%	Calculate inverse transform, ensuring a real output:
%	*	if warning is off, ignore the negative frequency components in transform
%	*	otherwise, take a standard inverse transform and report on error

	%	Initialise
	if nargin == 1 % the frequency dimension is unspecified
		dim = 1; % assume that the frequency dimension is 1
	end
	i = "runSup:LargeImaginary"; % warning ID
	s = warning("query", i); % check state of warning

	%	Inverse transform
	switch s.state % is warning state off or on?
		case "off" % assume that inverse transform is time-domain and therefore real
			r = ifft(r, [], dim, "symmetric"); % take inverse transform
		case "on" % take standard inverse transform, and check size of error
			r = ifft(r, [], dim); % take inverse transform
			mReal = max(abs(real(r(:)))); % maximum of real components
			mImag = max(abs(imag(r(:)))); % maximum of imaginary components
			ratio = mImag / mReal; % ratio of imaginary to real maxima
			if ratio > 1e-5 % issue warning
				w = "Inverse transform imaginary components are relatively large: %g";
				warning(i, w, ratio); % issue warning
				warning("off", i); % turn off further warnings with this ID
			end
			r = real(r); % make it real
	end