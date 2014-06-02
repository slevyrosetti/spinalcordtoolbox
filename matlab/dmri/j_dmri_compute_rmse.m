% =========================================================================
% FUNCTION
% j_dmri_compute_rmse.m
%
% Compute RMSE between i and i+1 image in case of multiple averaging.
%
% INPUT
% diff				structure generated by j_dmri_process_dti.m
% 
% OUTPUT
% (-)
%
% COMMENTS
% Julien Cohen-Adad 2009-09-05
% =========================================================================
function j_dmri_compute_rmse(diff)

% parameters
file_data		= 'dti_FA.nii.gz';
% ref_nex			= diff.nex; % reference to compute the RMSE

% % load ref file (with max numer of averaging)
% nifti = load_nifti([diff.folder_average{ref_nex},filesep,file_data]);
% data_ref = nifti.vol;
% [nx ny nz] = size(data_ref);
% data_ref_2d = reshape(data_ref,nx*ny*nz,1);



% loop over nex
j_progress('Compute RMSE ............................................')
for i_nex = 1:diff.nex-1
	
	% load file i
	fname_data1 = [diff.folder_average{i_nex},filesep,file_data];
	nifti = load_nifti(fname_data1);
	data1 = nifti.vol;
	[nx ny nz] = size(data1);
	% reshape in 2D
	data1_2d = reshape(data1,nx*ny*nz,1);
	% find zero voxels
	mask=find(data1_2d);

	% load file i+1
	fname_data2 = [diff.folder_average{i_nex+1},filesep,file_data];
	nifti = load_nifti(fname_data2);
	data2 = nifti.vol;
	data2_2d = reshape(data2,nx*ny*nz,1);
	
	fname_rmse = [diff.folder_average{i_nex},filesep,'rmse.nii.gz'];

	% compute RMSE
	rmse = sqrt(sum((data2_2d(mask)-data1_2d(mask)).^2)/(nx*ny));
	
	% save structure
	diff.rmse(i_nex) = rmse;
	j_progress(i_nex/(diff.nex-1))
end

% save structure
save([diff.path,filesep,'diff'],'diff');

% display picture
h_fig = figure;
plot(diff.rmse,'o-','linewidth',2)
title('Root mean square error on the FA (reference is nav=10)')
xlabel('nav')
ylabel('RMSE')
% print figure
print(h_fig,'-dpng','-r150',[diff.path,filesep,'fig_rmse','.png']);
% close

