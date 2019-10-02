function y = trim_end(x)

lx=zeros(1,length(x));
idx=length(x);
for i=length(x):-1:1
    if x(i) ~= 0
        idx = i;
        break;
    end
end
lx(1:idx)=1;
lx((idx+1):end)=0;
y = x(lx==1);
