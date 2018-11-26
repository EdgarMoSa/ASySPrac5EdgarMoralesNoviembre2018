function y=Dns(fun,a,b,n)
fa=subs(fun,a);
fb=subs(fun,b);
h=(b-a)/n;
s1=0;
s2=0;
for i=1:n-1
    r=a+i*h;
end
for i=1:2:n-1
    s1=s1+sub(fun,r);
end
for i=2:2:n-2
    s2=s2+sub(fun,r);
end
y=h*(fa+4*s1+2*s2+fb)/3;
    

