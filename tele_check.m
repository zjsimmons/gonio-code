function [corrected_tele] = tele_check(tele_char_array)
%this checks for 'special' characters.. was screwing up some laser
%communication
%correction logic here
%10->5E 4A (094 074)
%13->5E 4D (094 077)
%94->5E 9E (094 158)

corrected_tele(1,:)=tele_char_array(1,:);

for n=2:length(tele_char_array)-1
    
    if strcmp(tele_char_array(n,:),'010')
        corrected_entry=['094';'074'];
      %  disp('special character error')
    elseif strcmp(tele_char_array(n,:),'013')
        corrected_entry=['094';'077'];
      %  disp('special character error')
    elseif strcmp(tele_char_array(n,:),'094')
        corrected_entry=['094';'158'];
      %  disp('special character error')
    else
        corrected_entry=tele_char_array(n,:);
    end
    
    corrected_tele=[corrected_tele;corrected_entry];
    
end

corrected_tele=[corrected_tele;tele_char_array(end,:)];





end

