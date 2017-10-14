function s=Makespan(ch)

global nResource;
global nTask;
global TaskSize;
global R;
m=ch;
n=TaskSize;
r=R;
m1=zeros(1,nResource);
ms=zeros(1,nResource);

for i=1:nResource
    for j=1:nTask
        if m(j)==i
            m1(i)=m1(i)+n(j);
        end
    end
    ms(i)=(m1(i)+r(3,i))./r(1,i);
end
s=max(ms);
end   