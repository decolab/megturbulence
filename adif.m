function c=adif(a,b)
index=find(abs(a-b)>pi);
c=abs(a-b);
c(index)=2*pi-abs(a(index)-b(index));

