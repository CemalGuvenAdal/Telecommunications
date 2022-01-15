% EE431 Project 2
%1
%Creating random binary data
N=100000;
x=randi([0,1],1,N);
T=1
Ts=T/20;
t=[0:Ts:T-Ts]
p=sqrt(3)*t;
eksip=-1*p;
a=zeros(20*length(x),1);
for j=1:length(x)
    if x(j)==1
       a(20*(j-1)+1:j*20,1)=p;
    else
        a(20*(j-1)+1:j*20,1)=-p;
    end
end
s1=p;
b=zeros(N*20,1);
N0=logspace(-0.5,1,30);
errorarray=zeros(1,length(N0));
for k =1:length(N0)
% Noise adding
a=a+(N0(k)/2)*20*randn(length(a),1);
for j=1:N
    b(20*(j-1)+1:20*j,1)=(s1.').*(a(20*(j-1)+1:20*j,1));
end

for j=1:length(b)/20
    y(j)=trapz(1/N,b(20*(j-1)+1:20*j));
end
error=0;
y1=zeros(1,length(y));
for i=1:length(y)
    if y(1,i)>=0
        y1(1,i)=1;
    else
        y1(1,i)=0;
    end
end
        
var=sum(abs(y1-x))/length(x);
errorarray(1,k)=var;

end
figure();
semilogy(10*log10(1./N0),(errorarray))
grid on
title('Error vs SNR')
xlabel('SNR')
ylabel('Error')


% 
% % % part 2
N=100000;
x=randi([0,1],1,N);
T=1
Ts=T/10;
t=[0:Ts:T-Ts]
p=sqrt(3)*t;
eksip=-1*p;
a=zeros(20*length(x),1);
for j=1:length(x)
    if x(j)==1
       a(10*(j-1)+1:j*10,1)=p;
    else
        a(10*(j-1)+1:j*10,1)=-p;
    end
end
s1=p;
b=zeros(N*20,1);
N0=logspace(-0.5,1,30);
errorarray=zeros(1,length(N0));
for k =1:length(N0)
% Noise adding
a=a+(N0(k)/2)*20*randn(length(a),1);
for j=1:N
    b(10*(j-1)+1:10*j,1)=(s1.').*(a(10*(j-1)+1:10*j,1));
end

for j=1:length(b)/20
    y(j)=trapz(1/N,b(10*(j-1)+1:10*j));
end
error=0;
y1=zeros(1,length(y));
for i=1:length(y)
    if y(1,i)>=0
        y1(1,i)=1;
    else
        y1(1,i)=0;
    end
end
        
var=sum(abs(y1-x))/length(x);
errorarray(1,k)=var;

end
figure();
semilogy(10*log10(1./N0),(errorarray))
grid on
title('Error vs SNR')
xlabel('SNR')
ylabel('Error')
% part 3
N=100000;
x=randi([0,1],1,N);
T=1
Ts=T/20;
t=[0:Ts:T-Ts]
p=ones(1,20);
s1=sqrt(3)/2;
s0=-s1;
a=zeros(20*length(x),1);
for j=1:length(x)
    if x(j)==1
       a(20*(j-1)+1:j*20,1)=s1;
    else
        a(20*(j-1)+1:j*20,1)=s0;
    end
end

b=zeros(N*20,1);
N0=logspace(-0.5,1,30);
errorarray=zeros(1,length(N0));
for k =1:length(N0)
% Noise adding
a=a+(N0(k)/2)*20*randn(length(a),1);
for j=1:N
    b(20*(j-1)+1:20*j,1)=(p.').*(a(20*(j-1)+1:20*j,1));
end

for j=1:length(b)/20
    y(j)=trapz(1/N,b(20*(j-1)+1:20*j));
end
error=0;
y1=zeros(1,length(y));
for i=1:length(y)
    if y(1,i)>=0
        y1(1,i)=1;
    else
        y1(1,i)=0;
    end
end
        
var=sum(abs(y1-x))/length(x);
errorarray(1,k)=var;

end
figure();
semilogy(10*log10(1./N0),(errorarray))
grid on
title('Error vs SNR')
xlabel('SNR')
ylabel('Error')

% part 4
N=100000;
x=randi([0,1],1,N);
T=1
Ts=T/20;
t=[0:Ts:T-Ts]
p=sqrt(3)*t;
eksip=-1*p;
a=zeros(20*length(x),1);
for j=1:length(x)
    if x(j)==1
       a(20*(j-1)+1:j*20,1)=p;
    else
        a(20*(j-1)+1:j*20,1)=-p;
    end
end
s1=p;
b=zeros(N*20,1);
N0=logspace(-0.5,1,30);
errorarray=zeros(1,length(N0));
for k =1:length(N0)
% Noise adding
a=a+(N0(k)/2)*20*randn(length(a),1);
for j=1:N-1
    b(20*(j-1)+1:20*j,1)=(s1.').*(a(20*(j-1)+6:20*j+5,1));
end

for j=1:length(b)/20
    y(j)=trapz(1/N,b(20*(j-1)+1:20*j));
end
error=0;
y1=zeros(1,length(y));
for i=1:length(y)
    if y(1,i)>=0
        y1(1,i)=1;
    else
        y1(1,i)=0;
    end
end
        
var=sum(abs(y1-x))/length(x);
errorarray(1,k)=var;

end
figure();
semilogy(10*log10(1./N0),(errorarray))
grid on
title('Error vs SNR')
xlabel('SNR')
ylabel('Error')


% 
% % part 5
% 2
N0=logspace(-0.5,1,30);
S11=[-1,1];
signal=zeros(1,100000);
errore=0; 
varer=zeros(1,length(N0));
for  i= 1:100000
    a=randi([1,2]);
    signal(1,i)=S11(a);
end
for h =1:length(N0)
  

signal=signal+sqrt((N0(h)/2))*randn(1,100000);
signal2=zeros(1,100000);
for j= 1:100000
    [x,I]=min(abs(signal(1,j)-S11));
    signal2(1,j)=S11(I);
end
for i=1:100000
    if signal2(i)==signal(i)
        errore=errore;
    else
        errore=errore+1;
    end
end

varer(h)=errore/100000;   
end
figure()
semilogy(10*log10(2.5./N0),(varer))
title('2 PAM SNR VS ERROR');
xlabel('SNR')
ylabel('Error')

grid on 
% 4
N0=logspace(-0.5,1,30);
S11=[-3,-1,1,3];
signal=zeros(1,100000);
errore=0; 
varer=zeros(1,length(N0));
for  i= 1:100000
    a=randi([1,4]);
    signal(1,i)=S11(a);
end
for h =1:length(N0)
  

signal=signal+sqrt((N0(h)/2))*randn(1,100000);
signal2=zeros(1,100000);
for j= 1:100000
    [x,I]=min(abs(signal(1,j)-S11));
    signal2(1,j)=S11(I);
end
for i=1:100000
    if signal2(i)==signal(i)
        errore=errore;
    else
        errore=errore+1;
    end
end

varer(h)=errore/100000;   
end
figure()
semilogy(10*log10(2.5./N0),(varer))
grid on  
title('4 PAM SNR VS ERROR');
xlabel('SNR')
ylabel('Error')
% 8    
N0=logspace(-0.5,1,30);
S11=[-7,-5,-3,-1,1,3,5,7];
signal=zeros(1,100000);
errore=0; 
varer=zeros(1,length(N0));
for  i= 1:100000
    a=randi([1,8]);
    signal(1,i)=S11(a);
end
for h =1:length(N0)
  

signal=signal+sqrt((N0(h)/2))*randn(1,100000);
signal2=zeros(1,100000);
for j= 1:100000
    [x,I]=min(abs(signal(1,j)-S11));
    signal2(1,j)=S11(I);
end
for i=1:100000
    if signal2(i)==signal(i)
        errore=errore;
    else
        errore=errore+1;
    end
end

varer(h)=errore/100000;   
end
figure()
semilogy(10*log10(2.5./N0),(varer))
title('8 PAM SNR VS ERROR');
xlabel('SNR')
ylabel('Error')

grid on  
% 16
N0=logspace(-0.5,1,30);
S11=[-15,-13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15];
signal=zeros(1,100000);
errore=0; 
varer=zeros(1,length(N0));
for  i= 1:100000
    a=randi([1,16]);
    signal(1,i)=S11(a);
end
for h =1:length(N0)
  

signal=signal+sqrt((N0(h)/2))*randn(1,100000);
signal2=zeros(1,100000);
for j= 1:100000
    [x,I]=min(abs(signal(1,j)-S11));
    signal2(1,j)=S11(I);
end
for i=1:100000
    if signal2(i)==signal(i)
        errore=errore;
    else
        errore=errore+1;
    end
end

varer(h)=errore/100000;   
end
figure()
semilogy(10*log10(2.5./N0),(varer))
title('16 PAM SNR VS ERROR');
xlabel('SNR')
ylabel('Error')

grid on 
% 32
N0=logspace(-0.5,1,30);
S11=[-31,-29,-27,-25,-23,-21,-19,-17,-15,-13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31];
signal=zeros(1,100000);
errore=0; 
varer=zeros(1,length(N0));
for  i= 1:100000
    a=randi([1,32]);
    signal(1,i)=S11(a);
end
for h =1:length(N0)
  

signal=signal+sqrt((N0(h)/2))*randn(1,100000);
signal2=zeros(1,100000);
for j= 1:100000
    [x,I]=min(abs(signal(1,j)-S11));
    signal2(1,j)=S11(I);
end
for i=1:100000
    if signal2(i)==signal(i)
        errore=errore;
    else
        errore=errore+1;
    end
end

varer(h)=errore/100000;   
end
figure()
semilogy(10*log10(2.5./N0),(varer))
title('32 PAM SNR VS ERROR');
xlabel('SNR')
ylabel('Error')

grid on 
