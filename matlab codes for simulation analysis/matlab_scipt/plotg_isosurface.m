function [ output_args ] = plotg_isosurface( G, maxR, maxZ, dr, dz )
% plotg_isosurface is a plot function that plots the isosurface of level
% set function G with size nz*nr;
G = G';
nr = maxR/dr; nz = maxZ/dz;
r = 0:1:(nr-1);
r = r*dr;
z = 1:nz;
z = (z-1)*dz;

x = linspace(-maxR,maxR,100);
y = x;
[X,Y] = meshgrid(x,y);
[THETA,RHO] = cart2pol(X,Y);

% [THETA,RHO,Z] = meshgrid(theta,r,z);
% [X,Y,Z] = pol2cart(THETA,RHO,Z); %now X,Y,Z are not in meshed form.
% Xcoord = reshape(X,size(X,1)*size(X,2)*size(X,3),1);
% Ycoord = reshape(Y,size(Y,1)*size(Y,2)*size(Y,3),1);
% Zcoord = reshape(Z,size(Z,1)*size(Z,2)*size(Z,3),1);

% Xq = linspace(-6,6,100);
% Yq = Xq;
% Zq = z;
% [Xq,Yq,Zq] = meshgrid(Xq,Yq,Zq);
% Xqcoord = reshape(Xq,size(Xq,1)*size(Xq,2)*size(Xq,3),1);
% Yqcoord = reshape(Yq,size(Yq,1)*size(Yq,2)*size(Yq,3),1);
% Zqcoord = reshape(Zq,size(Zq,1)*size(Zq,2)*size(Zq,3),1);
% 
% foo = PetscBinaryRead(['outputG_t_',num2str(time,'%.6f'),'.h5']);
% G = reshape(foo,nr,nz);
rg = r;
for j = 1:nz
    g1 = G(:,j);
    interpG(:,:,j) = interp1(rg,g1,RHO);
end

% G3 = repmat(G,1,1,length(theta));
% G3 = permute(G3,[2 3 1]);
% G3 = reshape(G3,size(G3,1)*size(G3,2)*size(G3,3),1);
% 
% Vq = interp3(Xcoord,Ycoord,Zcoord,G3,Xqcoord,Yqcoord,Zqcoord);
[X,Y,Z] = meshgrid(x,z,y);
interpG = permute(interpG, [3 1 2]);
% X = permute(X,[3 1 2]);
% Y = permute(Y,[3 1 2]);
% Z = permute(Z,[3 1 2]);
figure
fv = isosurface(X,Y,Z,interpG,0);
p1 = patch(fv,'FaceColor','red','EdgeColor','none');
isonormals(X,Y,Z,interpG,p1);
view(3)
camlight
camlight('headlight')
lighting gouraud
axis equal
% xlim([-10 10])
% ylim([-10 10])
% zlim([100 400])
end

