 %%Time specifications:
   Fs = 1000;                   % samples per second
   dt = 1/Fs;                   % seconds per sample
   StopTime = 5;             % seconds
   t = (0:dt:StopTime-dt)';     % seconds
   %%Sine wave:

   frequencies = [0.1, 0.2, 0.4, 0.6, 0.8];
   f_list = zeros(length(t), length(frequencies));
   for i = 1:length(frequencies)
       f_list(:,i) = sin(2*pi*frequencies(i)*t);
   end

   figure;
   plot(t,f_list);
   
 
