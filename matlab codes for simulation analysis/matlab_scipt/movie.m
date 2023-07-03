clear; close all;
addpath('../matlab_scipt/')
maxr = 10; dr = 0.1; maxz = 500; dz = 0.4; nr = maxr/dr; nz = maxz/dz;
moviename = ['movie.gif'];
t_be=100; t_end=1000; t_diff=1; noffiles = (t_end-t_be)/t_diff+1;
bx = 6; by = maxz;
for frame = 1:noffiles
    time = t_diff*(frame-1) + t_be;
    G = PetscBinaryRead(['outputG_t_',num2str(time,'%.6f')]);
    G = reshape(G,nr,nz)';
    plotg_isosurface(G,maxr,maxz,dr,dz);
    xlim([-bx bx])
    zlim([-bx bx])
    ylim([0 by])
    title(['t=',num2str(time,'%.3f')]);
    movieFrame(frame) = getframe(gcf);
    close
    im = frame2im(movieFrame(frame));
    [imind,cm] = rgb2ind(im,256);
    if frame==1
        imwrite(imind,cm,moviename,'gif','DelayTime',0,'LoopCount',inf);
    else
        imwrite(imind,cm,moviename,'gif','DelayTime',0,'WriteMode','append');
    end
end