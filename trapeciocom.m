function integral=trapeciocom(f,a,b,n)
f0=subs(f,a);
fn=subs(f,b);
h=(b-a)/n;
fs=0;
for i=1:n-1
    x(i)=a+h*i;
    fs=fs+subs(f,x(i));
end
fs;
integral=0.5*h*(f0+2*fs+fn);