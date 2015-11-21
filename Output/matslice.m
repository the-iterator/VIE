% Matlab script to visualize total field with 3D slices
% To run:
% matlab2015a -nosplash -nodesktop -r matslice

load TotalField;

Dim = size(Etot);

NX = Dim(1);
NY = Dim(2);
NZ = Dim(3);

E1r = squeeze(real(Etot(:,:,:,1)));
E2r = squeeze(real(Etot(:,:,:,2)));
E3r = squeeze(real(Etot(:,:,:,3)));

Emag2 = abs(E1r).^2+abs(E2r).^2+abs(E3r).^2;

figure(1);
clf;

subplot(2,3,1)
h = slice(E1r,round(NX/2),[],[]);
h.EdgeColor = 'none';
hold on;

h = slice(E1r,[],round(NY/2),[]);
h.EdgeColor = 'none';

h = slice(E1r,[],[],round(NZ/2));
h.EdgeColor = 'none';

hold off;

axis equal tight;

subplot(2,3,2)
h = slice(E2r,round(NX/2),[],[]);
h.EdgeColor = 'none';
hold on;

h = slice(E2r,[],round(NY/2),[]);
h.EdgeColor = 'none';

h = slice(E2r,[],[],round(NZ/2));
h.EdgeColor = 'none';

hold off;

axis equal tight;

subplot(2,3,3)
h = slice(E3r,round(NX/2),[],[]);
h.EdgeColor = 'none';
hold on;

h = slice(E3r,[],round(NY/2),[]);
h.EdgeColor = 'none';

h = slice(E3r,[],[],round(NZ/2));
h.EdgeColor = 'none';

hold off;

axis equal tight;

subplot(2,3,5)
h = slice(Emag2,round(NX/2),[],[]);
h.EdgeColor = 'none';
hold on;

h = slice(Emag2,[],round(NY/2),[]);
h.EdgeColor = 'none';

h = slice(Emag2,[],[],round(NZ/2));
h.EdgeColor = 'none';

hold off;

axis equal tight;
