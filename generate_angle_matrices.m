function [thetas,phis] = generate_angle_matrices(analysis_settings,pp_coord,angle_out)
%generate_angle_matrices

%This function takes the process of generating the angle matrices out of
%the main program, saving time since the operation doesn't happen every run
%of the main program. 

%inputs: 
%analysis_settings structure, contains various settings
%pp_coord - the pupil plane coordinate list
%angle_out - corresponding angle list, i.e. asin mapping

%outputs: thetas and phis, i.e. the angle matrices..


%zjs 2018.1.22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_rows=analysis_settings.image_size(1,2);
no_columns=analysis_settings.image_size(1,1);
cir_x=analysis_settings.analysis_circle_x;
cir_y=analysis_settings.analysis_circle_y;
cir_rad=analysis_settings.analysis_circle_NA;
cir_rad_to_mask=analysis_settings.mask_r;

%pp_coord=xout;
%angle_out=yout;

for n=1:no_rows
    for m=1:no_columns
        y_fov=n-cir_y;
        x_fov=m-cir_x;
        %for theta determination, asin argument must be =<1
        argument=sqrt(x_fov^2+y_fov^2);
        if argument<=cir_rad_to_mask %to smaller mask circle.. 
            %This looks up for a given radius in pupil plane coordinates
            %the angle theta that it maps to. This is contained in the
            %input calibration pp_coord->angle_out.
            %Note: so this lininterp1 fn is much, about 50x faster            
            thetas(n,m)=lininterp1(pp_coord,angle_out,argument/(cir_rad*1.0));
            %problem here.. nice atan2 replaces a lot of conditional logic!
            phis(n,m)=atan2(y_fov,x_fov);            
        else
            %thetas and phis matrices, can't have a tan>1.. outside the
            %circle we set the values to nan.
            thetas(n,m)=nan;
            phis(n,m)=nan;
            %for mask too? 
            
        end
    end %no columns
end%no rows



%figure(102)
%imagesc(thetas)
%figure(103)
%imagesc(phis)

end

