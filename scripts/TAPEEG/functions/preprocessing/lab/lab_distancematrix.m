% Calculate distmatrixpol for function lab_detect_bad (FASTER-routines)

function [distmatrixpol,distmatrixxyz,distmatrixproj] = lab_distancematrix(EEG,P)
eeg_chans = 1:size(EEG,1);
distmatrix = zeros(length(eeg_chans),length(eeg_chans));
distmatrixpol = [];
for chan2tst = eeg_chans;
	for q=eeg_chans
		distmatrixpol(chan2tst,q)=sqrt(((P.radius(chan2tst)^2)+(P.radius(q)^2))-(2*((P.radius(chan2tst))*...
			(P.radius(q))*cosd(P.theta(chan2tst) - P.theta(q)))));%calculates the distance between electrodes using polar format
	end
end

locs = P;
for u = eeg_chans
	if ~isempty(P.x(u))
		Xs(u) = P.x(u);
	else
		Xs(u) = 0;
	end
	if ~isempty(P.y(u))
		Ys(u) = P.y(u);
	else
		Ys(u) = 0;
		end
	if ~isempty(P.z(u))
		Zs(u) = P.z(u);
	else
		Zs(u) = 0;
	end
end
Xs = round2(Xs,6);
Ys = round2(Ys,6);
Zs = round2(Zs,6);

for u = eeg_chans
	for v=eeg_chans
		distmatrixxyz(u,v) = dist(Xs(u),Xs(v))+dist(Ys(u),Ys(v))+dist(Zs(u),Zs(v));
	end
end
D = max(max(distmatrixxyz));
distmatrixproj = (pi-2*(acos(distmatrixxyz./D))).*(D./2);
	function d = dist(in1,in2)
		d = sqrt(abs(in1.^2 - in2.^2));
	end

	function num = round2(num,decimal)
		num = num .* 10^decimal;
		num = round(num);
		num = num ./ 10^decimal;
	end
end