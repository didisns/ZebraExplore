function [t, v, f_acq] = plotChannel(chStruct, ch)
%% chStruct contains the saved data. The structure has the 32
%% channels and each channel has value and time fieldnames
%% ch is the channel to extrapolate
name = fieldnames(chStruct);
v = chStruct.(name{ch+1}).value;
t_data =  chStruct.(name{ch+1}).time;
t_data (find(t_data(:,3) == 0),:) = [];
T_ACQ = 0;
t = [];
for i = 1:size(t_data,1)
    %%% da correggere
    %%t_data(i,2) = 1.9871e+04^-1;
%     t_acq = 2.0032e-04;% t_data(i,2);
    if (i == size(t_data,1))
        t_acq = T_ACQ/i;
    else
        t_acq =(t_data(i+1,1)-t_data(i,1))/ t_data(i,3);
%         if (t_acq > 0.5*t_data(i,2))  %% verifica che non ci sono salti di acquisizione
% %             keyboard
%             t_acq = t_data(i,2);
%         end
        T_ACQ = T_ACQ + t_acq;
    end
    t = [t, t_data(i,1)+([1:t_data(i,3)]-1)*t_acq];
end

f_acq = 1/t_acq;

v = v(1:length(t));
figure
plot(t,v);