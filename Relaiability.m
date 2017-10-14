function e=Price(x)

global R;
global TaskSize;
global nTask;
r=R;
e=0;
n=TaskSize;
  for i=1:nTask  
     e=e+((n(i)*r(3,x(i))));
  end
end