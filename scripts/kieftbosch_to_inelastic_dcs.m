clear;
% silicon
%M = load('~/cpo-share/inelastic_cs_silicon.txt');
%fid = fopen('~/scatter/inelastic-silicon.xml','w');
%rho_n = 4.99397e+28;
% sio2
%M = load('~/cpo-share/inelastic_cs_SiO2.txt');
%fid = fopen('~/scatter/inelastic-silicondioxide.xml','w');
%rho_n = 2.6541e+28;
% pmma
%M = load('~/cpo-share/inelastic_cs_PMMA.txt');
%fid = fopen('~/scatter/inelastic-pmma.xml','w');
%rho_n = 7.17e+27;
%
E = unique(M(:,1));
L = zeros(size(E));
scale = zeros(size(E));
DCS = zeros([size(E) 1000]);

x = [];
y = [];

fprintf(fid,'<cstable type="inelastic">\n');
for i = 1:numel(E)
	k = E(i);
	omega = M(M(:,1)==k,2);
	dcs = M(M(:,1)==k,3);
	mfp = M(M(:,1)==k,4);
	if mfp == inf
		continue
	end
	k = k*1.6e-19;
	omega = omega*1.6e-19;
	scaling = (1/mfp(end))/(rho_n*trapz(omega,dcs));
	dcs = dcs * scaling;
	%L(i) = 1./trapz(omega,dcs);
	L(i) = mfp(end);
	fprintf(fid,'<cross-section energy="%e*eV">\n',k/1.6e-19);
	for j = 1:numel(omega)
		fprintf(fid,'\t<insert omega0="%e*eV" dcs="%e*m^2" />\n',omega(j)/1.6e-19,dcs(j));
	end
	fprintf(fid,'</cross-section>\n');

	x = [x;k];
	y = [y;(trapz([0;omega],[0;dcs]))];
end
fprintf(fid,'</cstable>\n');
fclose(fid);

loglog(x/1.6e-19,y/1e-20);
