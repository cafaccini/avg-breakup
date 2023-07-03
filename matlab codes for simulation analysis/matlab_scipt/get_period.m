function [ period ] = get_period( input_args )

C:\Users\cfaccini\Documents\cap_breakup\20200124\maxz_500.0_maxr_10.0_dz_0.40_vf_2.00000_tension_1500.00_r0_0.97_tlow_1400.0_thigh_1850.0_twidth1_970.0_thigh_width_30.0_twidth2_50.0_reinit_10\data

%path = './maxz_500.0_maxr_10.0_dz_0.40_vf_2.00000_tension_1500.00_r0_0.97_tlow_1400.0_thigh_1850.0_twidth1_970.0_thigh_width_30.0_twidth2_50.0_reinit_10/data';
%path = '../20200124/maxz_500.0_maxr_10.0_dz_0.40_vf_2.00000_tension_1500.00_r0_0.97_tlow_1400.0_thigh_1850.0_twidth1_970.0_thigh_width_30.0_twidth2_50.0_reinit_10/data';

addpath('../matlab_scipt/')
maxr = 10; dr = 0.1; nr = maxr/dr;
maxz = 500; dz = 0.4; nz = maxz/dz;
time = 25;
G = PetscBinaryRead([path, 'outputG_t_',num2str(time,'%.6f')]);
G = reshape(G,nr,nz)';

g = G(:,11);
count = 1;
for i = 2:length(g)-1
    if g(i)<g(i-1) && g(i)<g(i+1) && g(i) < -0.05
        peak(count) = i;
        if count > 1 && peak(count)-peak(count-1)<20
            peak(count-1) = (peak(count-1)+peak(count))/2;
            count = count - 1;
        end
        count = count + 1;
    end
    
end

if(length(peak) < 2)
    error('double check!!!!')
end

period = (peak(2) - peak(1)) * dz;
plotg_isocurve(G, maxr, maxz, dr, dz)
end

