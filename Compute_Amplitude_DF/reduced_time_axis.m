function rt=reduced_time_axis(N)
  if isodd(N)
    H=(N-1)/2;
    rt=[-H:+H];
  else
    H=N/2-1;
    rt=[-H-1:H];
  end
end
