function w=hann_window(N,Fs)
  % Hann window (symmetric if N is odd, periodic if N is even)
  rt=reduced_time_axis(N);
  rT=N;  % time multiplied by Fs
  w=0.5*(1+cos(2*pi*rt/rT));
end
